!*******************************************************************************
!> author: Jacob Williams
!  date: October 29, 2016
!
!  Numerical differentiation module for computing the Jacobian matrix
!  (the derivative matrix of `m` functions w.r.t. `n` variables) using
!  finite differences.

    module numerical_differentiation_module

    use iso_fortran_env, only: wp => real64

    implicit none

    private

    real(wp),parameter :: zero = 0.0_wp

    type :: finite_diff_method

        !! defines the finite difference method
        !! used to compute the Jacobian.
        !!
        !! See: [[get_finite_difference_method]]
        !! for the different methods.

        private

        integer :: id = 0 !! unique ID for the method
        character(len=:),allocatable :: name  !! the name of the method
        integer :: class = 0 !! 2=backward diffs, 3=central diffs, etc...
        real(wp),dimension(:),allocatable :: dx_factors  !! multiplicative factors for `dx` perturbation
        real(wp),dimension(:),allocatable :: df_factors  !! multiplicative factors for accumulating function evaluations
        real(wp)                          :: df_den_factor = zero  !! denominator factor for finite difference equation (times `dx`)
    contains
        private
        procedure,public :: get_formula
    end type finite_diff_method
    interface finite_diff_method
        !! constructor
        module procedure initialize_finite_difference_method
    end interface

    type,public :: numdiff_type

        !! base type for sparsity and Jacobian computations.

        private

        integer :: n = 0 !! number of `x` variables
        integer :: m = 0 !! number of `f` functions

        real(wp),dimension(:),allocatable :: xlow  !! lower bounds on `x`
        real(wp),dimension(:),allocatable :: xhigh !! upper bounds on `x`

        integer :: chunk_size = 100  !! chuck size for allocating the arrays (>0)

        integer :: perturb_mode = 1  !! perturbation mode
                                     !! 1 - perturbation is `dx=dpert`
                                     !! 2 - perturbation is `dx=dpert*x`
                                     !! 3 - perturbation is `dx=dpert*(1+x)`
        real(wp),dimension(:),allocatable :: dpert !! perturbation vector for `x`

        logical :: sparsity_computed = .false. !! has the sparsity pattern already been computed?
        integer :: num_nonzero_elements = 0 !! number of nonzero elements in the jacobian
                                            !! (will be the dimension of `irow` and `icol`)
        integer,dimension(:),allocatable :: irow  !! sparsity pattern - rows of non-zero elements
        integer,dimension(:),allocatable :: icol  !! sparsity pattern - columns of non-zero elements

        type(finite_diff_method) :: meth    !! the finite difference method to use to
                                            !! compute the Jacobian

        ! these are required to be defined by the user:
        procedure(func),pointer    :: compute_function => null()
            !! the user-defined function

        procedure(spars_f),pointer :: compute_sparsity => null()
            !! for computing the sparsity pattern

        procedure(info_f),pointer :: info_function => null()
            !! an optional function the user can define
            !! which is called when each column of the jacobian is computed.
            !! It can be used to perform any setup operations.

    contains

        private

        procedure,public :: initialize => initialize_numdiff_type  !! initialize the class
        procedure,public :: compute_jacobian        !! main routine to compute the Jacobian
                                                    !! using the selected options. It
                                                    !! return the sparse (vector) form.
        procedure,public :: compute_jacobian_dense  !! return the dense `size(m,n)`
                                                    !! matrix form of the Jacobian.
        procedure,public :: destroy => destroy_numdiff_type  !! destroy the class
        procedure,public :: print_sparsity_pattern  !! print the sparsity pattern in vector form to a file
        procedure,public :: print_sparsity_matrix   !! print the sparsity pattern in matrix form to a file
        procedure,public :: set_sparsity_pattern    !! manually set the sparsity pattern

        ! internal routines:
        procedure :: destroy_sparsity            !! destroy the sparsity pattern
        procedure :: compute_perturbation_vector !! computes the variable perturbation factor
        procedure :: perturb_x_and_compute_f

    end type numdiff_type

    abstract interface
        subroutine func(me,x,f,indices_to_compute)
            !! The function (vector array of output functions `f`, computed
            !! from a vector of input variables `x`).
            !! This must be defined for all computations.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout) :: me
            real(wp),dimension(:),intent(in) :: x !! array of variables (size `n`)
            real(wp),dimension(:),intent(out) :: f !! array of functions (size `m`)
            integer,dimension(:),intent(in) :: indices_to_compute !! the elements of the
                                                                  !! function vector that need
                                                                  !! to be computed (the other
                                                                  !! are ignored)
        end subroutine func
        subroutine spars_f(me,x)
            !! The function to compute the sparsity pattern.
            !! It populates the `irow` and `icol` variables in the class.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout) :: me
            real(wp),dimension(:),intent(in) :: x !! vector of variables (size `n`)
        end subroutine spars_f
        subroutine info_f(me,column,i)
            !! User-defined info function (optional).
            !! Informs user what is being done during Jacobian computation.
            !! It can be used to perform any setup operations that need to
            !! done on the user's end.
            import :: numdiff_type
            implicit none
            class(numdiff_type),intent(inout)  :: me
            integer,intent(in) :: column !! the column being computed.
            integer,intent(in) :: i      !! perturbing this column for the `i`th time (1,2,...)
        end subroutine info_f
    end interface

    ! sparsity methods:
    public :: compute_sparsity_dense,compute_sparsity_random

    ! other:
    public :: get_finite_diff_formula

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Constructor for a [[finite_diff_method]].
!
!@note factors are input as integers for convienence, but are converted
!      to reals for the actual computations. (note: this means we can't
!      currently define methods that have non-integer factors).

    function initialize_finite_difference_method(id,name,class,dx_factors,&
                                                 df_factors,df_den_factor) result(me)

    implicit none

    type(finite_diff_method)        :: me
    integer,intent(in)              :: id            !! unique ID for the method
    character(len=*),intent(in)     :: name          !! the name of the method
    integer,intent(in)              :: class         !! 2=backward diffs, 3=central diffs, etc...
    integer,dimension(:),intent(in) :: dx_factors    !! multiplicative factors for dx perturbation
    integer,dimension(:),intent(in) :: df_factors    !! multiplicative factors for accumulating function evaluations
    integer,intent(in)              :: df_den_factor !! denominator factor for finite difference equation (times dx)

    if (size(dx_factors)/=size(df_factors)) then
        error stop 'Error: dx_factors and df_factors arrays must be the same size.'
    else

        me%id            = id
        me%name          = trim(name)
        me%class         = class
        me%dx_factors    = real(dx_factors,wp)
        me%df_factors    = real(df_factors,wp)
        me%df_den_factor = real(df_den_factor,wp)

    end if

    end function initialize_finite_difference_method
!*******************************************************************************

!*******************************************************************************
!>
!  Return a string with the finite difference formula.
!
!### Example
!  * For 3-point backward: `dfdx = (f(x-2h)-4f(x-h)+3f(x)) / (2h)`

    subroutine get_formula(me,formula)

    class(finite_diff_method),intent(in) :: me
    character(len=:),allocatable,intent(out) :: formula

    integer :: i !! counter
    integer :: istat !! write `iostat` flag
    character(len=10) :: x !! temp variable for integer to string conversion
    character(len=10) :: f !! temp variable for integer to string conversion

    if (allocated(me%dx_factors) .and. allocated(me%df_factors)) then

        formula = 'dfdx = ('

        do i = 1, size(me%dx_factors)

            if (int(me%df_factors(i))==1) then
                if (i==1) then
                    formula = formula//'f('
                else
                    formula = formula//'+f('
                end if
            elseif (int(me%df_factors(i))==-1) then
                formula = formula//'-f('
            else
                if (i==1) then
                    write(f,'(I10)',iostat=istat) int(me%df_factors(i))    ! integer to string
                else
                    write(f,'(SP,I10)',iostat=istat) int(me%df_factors(i)) ! integer to string (with sign)
                end if
                formula = formula//trim(adjustl(f))//'f('
            end if

            if (int(me%dx_factors(i))==0) then
                formula = formula//'x'
            elseif (int(me%dx_factors(i))==1) then
                formula = formula//'x+h'
            elseif (int(me%dx_factors(i))==-1) then
                formula = formula//'x-h'
            else
                write(x,'(SP,I10)',iostat=istat) int(me%dx_factors(i)) ! integer to string (with sign)
                formula = formula//'x'//trim(adjustl(x))//'h'
            end if

            formula = formula//')'

        end do

        write(f,'(I10)',iostat=istat) int(me%df_den_factor) ! integer to string
        if (int(me%df_den_factor)==1) then
            formula = formula//') / h'
        else
            formula = formula//') / ('//trim(adjustl(f))//'h)'
        end if

    else
        formula = ''
    end if

    end subroutine get_formula
!*******************************************************************************

!*******************************************************************************
!>
!  Return a string with the finite difference formula.
!  Input is the method `id` code.
!
!###See also:
!  * [[get_formula]]

    subroutine get_finite_diff_formula(id,formula)

    implicit none

    integer,intent(in) :: id  !! the id code for the method
    character(len=:),allocatable,intent(out) :: formula !! the formula string

    type(finite_diff_method) :: fd

    call get_finite_difference_method(id,fd)
    call get_formula(fd,formula)

    end subroutine get_finite_diff_formula
!*******************************************************************************

!*******************************************************************************
!>
!  Return a [[finite_diff_method]] given the `id` code.
!  (the `id` codes begin at 1, are sequential, and uniquely define the method).
!
!@note This is the only routine that has to be changed if a new
!      finite difference method is added.

    subroutine get_finite_difference_method(id,fd)

    implicit none

    integer,intent(in) :: id  !! the id code for the method
    type(finite_diff_method),intent(out) :: fd !! this method (can be used in [[compute_jacobian]])

    select case (id)
    case(1); fd = finite_diff_method(id,'2-point forward',  2,[1,0],[1,-1],1)      ! (f(x+h) - f(x)) / h
    case(2); fd = finite_diff_method(id,'2-point backward', 2,[0,-1],[1,-1],1)     ! (f(x) - f(x-h)) / h
    case(3); fd = finite_diff_method(id,'3-point forward',  3,[0,1,2],[-3,4,-1],2)
    case(4); fd = finite_diff_method(id,'3-point backward', 3,[-2,-1,0],[1,-4,3],2)
    case(5); fd = finite_diff_method(id,'3-point central',  3,[1,-1],[1,-1],2)     ! (f(x+h) - f(x-h)) / (2h)

    case default
        error stop 'Error: invalid id.'
    end select

    end subroutine get_finite_difference_method
!*******************************************************************************

!*******************************************************************************
!>
!  Initialize a [[numdiff_type]] class. This must be called first.

    subroutine initialize_numdiff_type(me,n,m,xlow,xhigh,perturb_mode,dpert,&
                        problem_func,sparsity_func,jacobian_method,info,chunk_size)

    implicit none

    class(numdiff_type),intent(out)     :: me
    integer,intent(in)                  :: n               !! number of `x` variables
    integer,intent(in)                  :: m               !! number of `f` functions
    real(wp),dimension(n),intent(in)    :: xlow            !! lower bounds on `x`
    real(wp),dimension(n),intent(in)    :: xhigh           !! upper bounds on `x`
    integer,intent(in)                  :: perturb_mode    !! perturbation mode (1,2,3)
    real(wp),dimension(n),intent(in)    :: dpert           !! perturbation vector for `x`
    procedure(func)                     :: problem_func    !!
    procedure(spars_f)                  :: sparsity_func   !!
    integer,intent(in)                  :: jacobian_method !! `id` code for the finite difference method
                                                           !! see [[get_finite_difference_method]]
    procedure(info_f),optional          :: info            !! a function the user can define
                                                           !! which is called when each column
                                                           !! of the jacobian is computed.
                                                           !! It can be used to perform any
                                                           !! setup operations.
    integer,intent(in),optional         :: chunk_size      !! chunk size for allocating the arrays
                                                           !! (must be >0) [default is 100]

    ! functions:
    me%compute_function => problem_func
    me%compute_sparsity => sparsity_func

    ! method:
    call get_finite_difference_method(jacobian_method,me%meth)

    ! size of the problem:
    me%n = n
    me%m = m

    ! input variable bounds:
    if (allocated(me%xlow)) deallocate(me%xlow)
    if (allocated(me%xhigh)) deallocate(me%xhigh)
    allocate(me%xlow(n))
    allocate(me%xhigh(n))
    me%xlow = xlow
    me%xhigh = xhigh

    ! perturbation options:
    me%perturb_mode = perturb_mode
    if (allocated(me%dpert)) deallocate(me%dpert)
    allocate(me%dpert(n))
    me%dpert = dpert

    ! optional:
    if (present(info))       me%info_function => info
    if (present(chunk_size)) me%chunk_size = chunk_size

    end subroutine initialize_numdiff_type
!*******************************************************************************

!*******************************************************************************
!>
!  destroy the [[numdiff_type]] class.

    subroutine destroy_numdiff_type(me)

    implicit none

    class(numdiff_type),intent(out) :: me

    end subroutine destroy_numdiff_type
!*******************************************************************************

!*******************************************************************************
!>
!  destroy the sparsity pattern in the class.

    subroutine destroy_sparsity(me)

    implicit none

    class(numdiff_type),intent(inout) :: me

    me%sparsity_computed = .false.
    me%num_nonzero_elements = 0
    if (allocated(me%irow)) deallocate(me%irow)
    if (allocated(me%icol)) deallocate(me%icol)

    end subroutine destroy_sparsity
!*******************************************************************************

!*******************************************************************************
!>
!  To specify the sparsity pattern directly if it is already known.

    subroutine set_sparsity_pattern(me,irow,icol)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,dimension(:),intent(in) :: irow
    integer,dimension(:),intent(in) :: icol

    call me%destroy_sparsity()

    if (size(irow)/=size(icol) .or. any(irow>me%m) .or. any(icol>me%n)) then
        error stop 'Error: invalid inputs to set_sparsity_pattern'
    else
        me%sparsity_computed = .true.
        me%num_nonzero_elements = size(irow)
        me%irow = irow
        me%icol = icol
    end if

    end subroutine set_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  assume all elements of Jacobian are non-zero.

    subroutine compute_sparsity_dense(me,x)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x !! vector of variables (size `n`)

    integer :: i !! counter
    integer :: r !! row counter
    integer :: c !! column counter

    call me%destroy_sparsity()

    me%num_nonzero_elements = me%m * me%n
    allocate(me%irow(me%num_nonzero_elements))
    allocate(me%icol(me%num_nonzero_elements))

    ! create the dense matrix:
    i = 0
    do c = 1, me%n
        do r = 1, me%m
            i = i + 1
            me%irow(i) = r
            me%icol(i) = c
        end do
    end do

    me%sparsity_computed = .true.

    end subroutine compute_sparsity_dense
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the sparsity pattern by computing the function at three
!  "random" points in the [xlow,xhigh] interval and checking if the
!  function values are the same.
!
!@note The input `x` is not used here.
!
!@note Could also allow the three coefficients to be user inputs.

    subroutine compute_sparsity_random(me,x)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x !! vector of variables (size `n`)

    integer :: i !! column counter
    integer :: j !! row counter
    integer :: n_icol  !! `icol` size counter
    integer :: n_irow  !! `irow` size counter
    integer,dimension(me%m) :: idx !! indices to compute [1,2,...,m]
    real(wp),dimension(me%n) :: x1 !! perturbed variable vector
    real(wp),dimension(me%n) :: x2 !! perturbed variable vector
    real(wp),dimension(me%n) :: x3 !! perturbed variable vector
    real(wp),dimension(me%m) :: f1 !! function evaluation
    real(wp),dimension(me%m) :: f2 !! function evaluation
    real(wp),dimension(me%m) :: f3 !! function evaluation

    real(wp),dimension(3),parameter :: coeffs = [0.251234567_wp,&
                                                 0.512345678_wp,&
                                                 0.751234567_wp]
        !! Pick three points roughly equally spaced.
        !! (add some noise in attempt to avoid freak zeros)
        !!````
        !! xlow---|----|--x--|---xhigh
        !!        1    2     3
        !!````

    ! initialize:
    call me%destroy_sparsity()

    ! we will compute all the functions:
    idx = [(i,i=1,me%m)]

    n_icol = 0  ! initialize vector size counters
    n_irow = 0

    ! define a nominal point roughly in the middle:
    x2 = me%xlow + (me%xhigh-me%xlow)*coeffs(2)
    call me%compute_function(x2,f2,idx)

    do i = 1, me%n  ! columns

        ! restore nominal:
        x1 = x2
        x3 = x2

        x1(i) = me%xlow(i) + (me%xhigh(i)-me%xlow(i))*coeffs(1)
        x3(i) = me%xlow(i) + (me%xhigh(i)-me%xlow(i))*coeffs(3)

        call me%compute_function(x1,f1,idx)
        call me%compute_function(x3,f3,idx)

        do j = 1, me%m ! each function (rows of Jacobian)
            if (f1(j)/=f2(j) .or. f3(j)/=f2(j)) then
                call expand_vector(me%icol,n_icol,me%chunk_size,val=i)
                call expand_vector(me%irow,n_irow,me%chunk_size,val=j)
            end if
        end do
        ! resize to correct size:
        call expand_vector(me%icol,n_icol,me%chunk_size,finished=.true.)
        call expand_vector(me%irow,n_irow,me%chunk_size,finished=.true.)

    end do

    ! finished:
    me%sparsity_computed = .true.
    me%num_nonzero_elements = size(me%irow)

contains

    pure subroutine expand_vector(vec,n,chunk_size,val,finished)

    !! add elements to the vector in chunks.

    implicit none

    integer,dimension(:),allocatable,intent(inout) :: vec !! the vector to add element to
    integer,intent(inout) :: n     !! counter for last element added to `vec`.
                                   !! must be initialized to `size(vec)`
                                   !! (or 0 if not allocated) before first call
    integer,intent(in) :: chunk_size  !! allocate `vec` in blocks of this size (>0)
    integer,intent(in),optional :: val !! the value to add to `vec`
    logical,intent(in),optional :: finished !! set to true to return `vec`
                                            !! as its correct size (`n`)

    integer,dimension(:),allocatable :: tmp  !! temporary array

    if (present(val)) then
        if (allocated(vec)) then
            if (n==size(vec)) then
                ! have to add another chunk:
                allocate(tmp(size(vec)+chunk_size))
                tmp(1:size(vec)) = vec
                call move_alloc(tmp,vec)
            end if
            n = n + 1
        else
            ! the first element:
            allocate(vec(chunk_size))
            n = 1
        end if
        vec(n) = val
    end if

    if (present(finished)) then
        if (finished) then
            ! set vec to actual size (n):
            if (allocated(tmp)) deallocate(tmp)
            allocate(tmp(n))
            tmp = vec(1:n)
            call move_alloc(tmp,vec)
        end if
    end if

    end subroutine expand_vector

    end subroutine compute_sparsity_random
!*******************************************************************************

!*******************************************************************************
!>
!  just a wrapper for `compute_jacobian`, that returns a dense matrix.

    subroutine compute_jacobian_dense(me,x,jac)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x !! vector of variables (size `n`)
    real(wp),dimension(:,:),allocatable,intent(out) :: jac !! the jacobian matrix

    real(wp),dimension(:),allocatable :: jac_vec  !! sparse jacobian representation
    integer :: i !! counter

    ! compute sparse form of jacobian:
    call me%compute_jacobian(x,jac_vec)

    ! size output matrix:
    allocate(jac(me%m,me%n))

    ! convert to dense form:
    jac = zero
    do i = 1, me%num_nonzero_elements
        jac(me%irow(i),me%icol(i)) = jac_vec(i)
    end do

    end subroutine compute_jacobian_dense
!*******************************************************************************

!*******************************************************************************
!>
!  Perturb the specified optimization variable, and compute the function.
!  This routine is designed so that `df` is accumulated as each function
!  evaluation is done, to avoid having to allocate more temporary storage.

    subroutine perturb_x_and_compute_f(me,x,dx_factor,dx,&
                                       df_factor,column,idx,df,df_den_factor)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x         !! nominal variable vector
    real(wp),intent(in)               :: dx_factor !! factor to multiply `dx`
    real(wp),dimension(:),intent(in)  :: dx        !! the perturbation value for this column
    real(wp),intent(in)               :: df_factor !! factor to multiply function value
    integer,intent(in)                :: column    !! the variable to perturb
    integer,dimension(:),intent(in)   :: idx       !! the elements in this
                                                   !! column of the Jacobian
                                                   !! to compute (passed to function)
    real(wp),dimension(me%m),intent(inout) :: df   !! the accumulated function value
                                                   !! note: for the first call, this
                                                   !! should be set to zero
    real(wp),intent(in),optional :: df_den_factor  !! if present, `df` is divided by
                                                   !! `df_den_factor*dx(column)`

    real(wp),dimension(me%n) :: xp  !! the perturbed variable vector
    real(wp),dimension(me%m) :: f   !! function evaluation

    xp = x
    if (dx_factor/=zero) xp(column) = xp(column) + dx_factor * dx(column)
    call me%compute_function(xp,f,idx)
    df(idx) = df(idx) + df_factor * f(idx)
    if (present(df_den_factor)) df(idx) = df(idx) / (df_den_factor*dx(column))

    end subroutine perturb_x_and_compute_f
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian using forward differences.

    subroutine compute_jacobian(me,x,jac)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x  !! vector of variables (size `n`)
    real(wp),dimension(:),allocatable,intent(out) :: jac  !! sparse jacobian vector

    real(wp),dimension(:),allocatable :: dfdx    !! the non-zero elements of a
                                                 !! column of the Jacobian matrix
    logical,dimension(me%m) :: funcs_to_compute  !! the elements in a given
                                                 !! column of the Jacobian
                                                 !! that are non-zero
    integer :: i  !! column counter
    real(wp),dimension(me%n) :: dx !! absolute perturbation (>0) for each variable
    integer,dimension(:),allocatable :: nonzero_elements_in_col !! the indices of the
                                                                !! nonzero Jacobian
                                                                !! elements in a column
    integer,dimension(:),allocatable :: indices  !! index vector
                                                 !! `[1,2,...,num_nonzero_elements]`
                                                 !! for putting `dfdx` into `jac`
    integer :: j !! function evaluation counter
    real(wp),dimension(me%m) :: df  !! accumulated function

    ! if we don't have a sparsity pattern yet then compute it:
    if (.not. me%sparsity_computed) call me%compute_sparsity(x)

    ! initialize:
    allocate(jac(me%num_nonzero_elements))
    jac = zero
    indices = [(i,i=1,me%num_nonzero_elements)] !NOTE could save this in the class so we don't have to keep allocating it

    ! compute dx vector:
    call me%compute_perturbation_vector(x,dx)

    ! compute Jacobian matrix column-by-column:
    do i=1,me%n

        ! determine functions to compute for this column:
        nonzero_elements_in_col = pack(me%irow,mask=me%icol==i)
        if (size(nonzero_elements_in_col)/=0) then ! there are functions to compute

            ! compute this column of the Jacobian:
            df = zero
            do j = 1, size(me%meth%dx_factors)-1
                if (associated(me%info_function)) call me%info_function(i,j)
                call me%perturb_x_and_compute_f(x,me%meth%dx_factors(j),&
                                                dx,me%meth%df_factors(j),&
                                                i,nonzero_elements_in_col,df)
            end do
            ! the last one has the denominator:
            if (associated(me%info_function)) call me%info_function(i,j)
            call me%perturb_x_and_compute_f(x,me%meth%dx_factors(j),&
                                            dx,me%meth%df_factors(j),&
                                            i,nonzero_elements_in_col,df,&
                                            me%meth%df_den_factor)

            ! put result into the output vector:
            jac(pack(indices,mask=me%icol==i)) = df(nonzero_elements_in_col)

        end if

    end do

    end subroutine compute_jacobian
!*******************************************************************************

!*******************************************************************************
!>
!  Compute `dx`, the perturbation vector for `x` used
!  when computing the gradients.

    subroutine compute_perturbation_vector(me,x,dx)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(me%n),intent(in)  :: x  !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(out) :: dx !! absolute perturbation (>0) for each variable

    real(wp),parameter :: eps = epsilon(1.0_wp) !! the smallest allowed absolute step

    select case (me%perturb_mode)
    case(1)
        dx = abs(me%dpert)
    case(2)
        dx = abs(me%dpert * x)
    case(3)
        dx = abs(me%dpert) * (1.0_wp + abs(x))
    case default
        error stop 'Error: invalid value for perturb_mode (must be 1, 2, or 3)'
    end select

    ! make sure none are too small:
    where (dx<eps) dx = eps

    end subroutine compute_perturbation_vector
!*******************************************************************************

!*******************************************************************************
!>
!  Print the sparsity pattern in vector form (`irow`, `icol`).

    subroutine print_sparsity_pattern(me,iunit)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in) :: iunit !! file unit to write to.
                                !! (assumed to be already opened)

    if (allocated(me%irow) .and. allocated(me%icol)) then
        write(iunit,'(A,1X,*(I3,","))') 'irow=',me%irow
        write(iunit,'(A,1X,*(I3,","))') 'icol=',me%icol
    else
        error stop 'Error: sparsity pattern not available.'
    end if

    end subroutine print_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  Print the sparsity pattern in matrix form.

    subroutine print_sparsity_matrix(me,iunit)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in) :: iunit !! file unit to write to.
                                !! (assumed to be already opened)

    integer :: r   !! row counter
    character(len=1),dimension(me%n) :: row  !! a row of the sparsity matrix

    if (allocated(me%irow) .and. allocated(me%icol)) then
        ! print by row:
        do r = 1,me%m
            row = '0'
            row(pack(me%icol,mask=me%irow==r)) = 'X'
            write(iunit,'(*(A1))') row
        end do
    else
        error stop 'Error: sparsity pattern not available.'
    end if

    end subroutine print_sparsity_matrix
!*******************************************************************************

!*******************************************************************************
    end module numerical_differentiation_module
!*******************************************************************************
