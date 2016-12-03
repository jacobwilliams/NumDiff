!*******************************************************************************
!> author: Jacob Williams
!  date: October 29, 2016
!
!  Numerical differentiation module for computing the Jacobian matrix
!  (the derivative matrix of `m` functions w.r.t. `n` variables) using
!  finite differences.

    module numerical_differentiation_module

    use kinds_module
    use utilities_module
    use dsm_module,      only: dsm
    use iso_fortran_env, only: error_unit
    use diff_module,     only: diff_func
    use cache_module,    only: function_cache

    implicit none

    private

    real(wp),parameter :: zero = 0.0_wp

    type,public :: finite_diff_method

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
        procedure,public :: print => print_finite_difference_method
    end type finite_diff_method
    interface finite_diff_method
        !! constructor
        module procedure initialize_finite_difference_method
    end interface

    type,public :: meth_array
        !! to store an array of [[finite_diff_method]] types
        !! this is used when the `mode=2` option is used
        !! in [[numdiff_type]]
        private
        type(finite_diff_method),dimension(:),allocatable :: meth
    end type meth_array

    type,public :: sparsity_pattern

        !! A sparsity pattern

        private

        logical :: sparsity_computed = .false. !! has the sparsity pattern already been computed?
        integer :: num_nonzero_elements = 0 !! number of nonzero elements in the jacobian
                                            !! (will be the dimension of `irow` and `icol`)
        integer,dimension(:),allocatable :: irow  !! sparsity pattern - rows of non-zero elements
        integer,dimension(:),allocatable :: icol  !! sparsity pattern - columns of non-zero elements

        integer,dimension(:),allocatable :: indices  !! index vector
                                                     !! `[1,2,...,num_nonzero_elements]`
                                                     !! for putting `df` into `jac`

        integer :: maxgrp = 0 !! the number of groups in the partition
                              !! of the columns of `a`.
        integer,dimension(:),allocatable :: ngrp !! specifies the partition of the columns of `a`.
                                                 !! column `jcol` belongs to group `ngrp(jcol)`.
                                                 !! `size(n)`

        contains
        private
        procedure,public :: destroy => destroy_sparsity
        procedure,public :: print => print_sparsity
        procedure,public :: columns_in_partition_group
    end type sparsity_pattern

    type,public :: numdiff_type

        !! base type for sparsity and Jacobian computations.

        private

        integer :: n = 0 !! number of `x` variables
        integer :: m = 0 !! number of `f` functions

        real(wp),dimension(:),allocatable :: xlow  !! lower bounds on `x`
        real(wp),dimension(:),allocatable :: xhigh !! upper bounds on `x`

        integer :: chunk_size = 100  !! chuck size for allocating the arrays (>0)

        integer :: perturb_mode = 1  !! perturbation mode:
                                     !! **1** - perturbation is `dx=dpert`,
                                     !! **2** - perturbation is `dx=dpert*x`,
                                     !! **3** - perturbation is `dx=dpert*(1+x)`
        real(wp),dimension(:),allocatable :: dpert !! perturbation vector for `x`

        logical :: partition_sparsity_pattern = .false.  !! to partition the sparsity pattern using [[dsm]]
        type(sparsity_pattern) :: sparsity  !! the sparsity pattern

        integer :: mode = 1 !! **1** = use `meth` (specified methods),
                            !! **2** = use `class` (specified class, method is selected on-the-fly).
        type(finite_diff_method),dimension(:),allocatable :: meth   !! the finite difference method to use
                                                                    !! compute the `n`th column of the Jacobian
                                                                    !! `size(n)`.  Either this or `class` is used
        integer,dimension(:),allocatable :: class  !! the class of method to use to
                                                   !! compute the `n`th column of the Jacobian
                                                   !! `size(n)`. Either this or `meth` is used
        type(meth_array),dimension(:),allocatable :: class_meths !! array of methods for the specified classes.
                                                                 !! used with `class` when `mode=2`


        ! parameters when using diff:
        real(wp) :: eps = 1.0e-9_wp !! tolerance parameter for [[diff]]
        real(wp) :: acc = 0.0_wp    !! tolerance parameter for [[diff]]

        type(function_cache) :: cache    !! if using the function cache

        logical :: use_diff = .false. !! if we are using the Neville's process method,
                                      !! rather than finite differences

        procedure(func),pointer :: problem_func => null()
            !! the user-defined function

        ! these are required to be defined by the user:
        procedure(func),pointer :: compute_function => null()
            !! compute the user-defined function
            !! this can point to the `problem_func` or, if using
            !! the cache, it points to [[compute_function_with_cache]].

        procedure(spars_f),pointer :: compute_sparsity => null()
            !! for computing the sparsity pattern

        procedure(info_f),pointer :: info_function => null()
            !! an optional function the user can define
            !! which is called when each column of the jacobian is computed.
            !! It can be used to perform any setup operations.

        procedure(jacobian_f),pointer :: jacobian_function => null()
            !! the low-level function called by [[compute_jacobian]]
            !! that actually computes the jacobian.

    contains

        private

        procedure,public :: initialize => initialize_numdiff !! initialize the class
        procedure,public :: diff_initialize => initialize_numdiff_for_diff !! initialize the class

        procedure,public :: compute_jacobian              !! main routine to compute the Jacobian
                                                          !! using the selected options. It
                                                          !! returns the sparse (vector) form.
        procedure,public :: compute_jacobian_dense        !! return the dense `size(m,n)`
                                                          !! matrix form of the Jacobian.
        procedure,public :: compute_jacobian_times_vector !! returns the product of the Jacobian
                                                          !! matrix and an input vector
        procedure,public :: destroy => destroy_numdiff_type  !! destroy the class
        procedure,public :: print_sparsity_pattern  !! print the sparsity pattern in vector form to a file
        procedure,public :: print_sparsity_matrix   !! print the sparsity pattern in matrix form to a file
        procedure,public :: set_sparsity_pattern    !! manually set the sparsity pattern
        procedure,public :: select_finite_diff_method  !! select a method in a specified class so
                                                       !! that the variable bounds are not violated
                                                       !! when by the perturbations.
        procedure,public :: select_finite_diff_method_for_partition_group  !! version of [[select_finite_diff_method]]
                                                                           !! for partitioned sparsity pattern.

        ! internal routines:
        procedure :: destroy_sparsity_pattern      !! destroy the sparsity pattern
        procedure :: compute_perturbation_vector   !! computes the variable perturbation factor
        procedure :: perturb_x_and_compute_f
        procedure :: perturb_x_and_compute_f_partitioned

    end type numdiff_type

    abstract interface
        subroutine func(me,x,f,funcs_to_compute)
            !! The function (vector array of output functions `f`, computed
            !! from a vector of input variables `x`).
            !! This must be defined for all computations.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout) :: me
            real(wp),dimension(:),intent(in) :: x !! array of variables (size `n`)
            real(wp),dimension(:),intent(out) :: f !! array of functions (size `m`)
            integer,dimension(:),intent(in) :: funcs_to_compute !! the elements of the
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
            integer,dimension(:),intent(in) :: column !! the columns being computed.
            integer,intent(in) :: i  !! perturbing these columns for the `i`th time (1,2,...)
        end subroutine info_f
        subroutine jacobian_f(me,x,dx,jac)
            !! Actual function for computing the Jacobian
            !! called by [[compute_jacobian]].
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout)   :: me
            real(wp),dimension(:),intent(in)    :: x    !! vector of variables (size `n`)
            real(wp),dimension(me%n),intent(in) :: dx   !! absolute perturbation (>0) for each variable
            real(wp),dimension(:),intent(out)   :: jac  !! sparse jacobian vector (size `num_nonzero_elements`)
        end subroutine jacobian_f

    end interface

    ! other:
    public :: get_finite_diff_formula
    public :: get_all_methods_in_class

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
!  Print the contents of a [[finite_diff_method]]. Used for debugging.

    subroutine print_finite_difference_method(me,iunit)

    implicit none

    class(finite_diff_method),intent(in) :: me
    integer,intent(in)                   :: iunit  !! file unit for printing
                                                   !! (assumed to be opened)

    write(iunit,'(A,1X,I5)')        'id            :', me%id
    write(iunit,'(A,1X,A)')         'name          :', me%name
    write(iunit,'(A,1X,I5)')        'class         :', me%class
    write(iunit,'(A,1X,*(I5,","))') 'dx_factors    :', int(me%dx_factors)
    write(iunit,'(A,1X,*(I5,","))') 'df_factors    :', int(me%df_factors)
    write(iunit,'(A,1X,I5)')        'df_den_factor :', int(me%df_den_factor)

    end subroutine print_finite_difference_method
!*******************************************************************************

!*******************************************************************************
!>
!  Wrapper for computing the function, using the cache.

    subroutine compute_function_with_cache(me,x,f,funcs_to_compute)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x !! array of variables (size `n`)
    real(wp),dimension(:),intent(out) :: f !! array of functions (size `m`)
    integer,dimension(:),intent(in) :: funcs_to_compute !! the elements of the
                                                        !! function vector that need
                                                        !! to be computed (the other
                                                        !! are ignored)

    integer :: i !! index in the cache
    logical,dimension(size(funcs_to_compute)) :: ffound  !! functions found in the cache
    logical :: xfound  !! if `x` was found in the cache

    call me%cache%get(x,funcs_to_compute,i,f,xfound,ffound)

    if (xfound .and. any(ffound)) then

        ! at least one of the functions was found
        if (all(ffound)) return ! all were found

        ! compute the ones that weren't found,
        ! and add them to the cache:
        call me%problem_func(x,f,pack(funcs_to_compute,mask=(.not. ffound)))
        call me%cache%put(i,x,f,pack(funcs_to_compute,mask=(.not. ffound)))

    else

        ! compute the function and add it to the cache:
        call me%problem_func(x,f,funcs_to_compute)
        call me%cache%put(i,x,f,funcs_to_compute)

    end if

    end subroutine compute_function_with_cache
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
    logical :: found

    call get_finite_difference_method(id,fd,found)
    call get_formula(fd,formula)

    end subroutine get_finite_diff_formula
!*******************************************************************************

!*******************************************************************************
!>
!  Return a [[finite_diff_method]] given the `id` code.
!  (the `id` codes begin at 1, are sequential, and uniquely define the method).
!
!  The available methods are:
!
!   * \( (f(x+h)-f(x)) / h                                       \)
!   * \( (f(x)-f(x-h)) / h                                       \)
!   * \( (f(x+h)-f(x-h)) / (2h)                                  \)
!   * \( (-3f(x)+4f(x+h)-f(x+2h)) / (2h)                         \)
!   * \( (f(x-2h)-4f(x-h)+3f(x)) / (2h)                          \)
!   * \( (-2f(x-h)-3f(x)+6f(x+h)-f(x+2h)) / (6h)                 \)
!   * \( (f(x-2h)-6f(x-h)+3f(x)+2f(x+h)) / (6h)                  \)
!   * \( (-11f(x)+18f(x+h)-9f(x+2h)+2f(x+3h)) / (6h)             \)
!   * \( (-2f(x-3h)+9f(x-2h)-18f(x-h)+11f(x)) / (6h)             \)
!   * \( (f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h)) / (12h)               \)
!   * \( (-3f(x-h)-10f(x)+18f(x+h)-6f(x+2h)+f(x+3h)) / (12h)     \)
!   * \( (-f(x-3h)+6f(x-2h)-18f(x-h)+10f(x)+3f(x+h)) / (12h)     \)
!   * \( (-25f(x)+48f(x+h)-36f(x+2h)+16f(x+3h)-3f(x+4h)) / (12h) \)
!   * \( (3f(x-4h)-16f(x-3h)+36f(x-2h)-48f(x-h)+25f(x)) / (12h)  \)
!
!  Where \(f(x)\) is the user-defined function of \(x\)
!  and \(h\) is a "small" perturbation.
!
!@note This is the only routine that has to be changed if a new
!      finite difference method is added.
!
!@note The order within a class is assumed to be the order that we would perfer
!      to use them (e.g., central diffs are first, etc.) This is used in
!      the [[select_finite_diff_method]] routine.

    subroutine get_finite_difference_method(id,fd,found)

    implicit none

    integer,intent(in)                   :: id     !! the id code for the method
    type(finite_diff_method),intent(out) :: fd     !! this method (can be used in [[compute_jacobian]])
    logical,intent(out)                  :: found  !! true if it was found

    found = .true.

    select case (id)
    case(1)
        ! (f(x+h)-f(x)) / h
        fd = finite_diff_method(id,'2-point forward',    2,[1,0],[1,-1],1)
    case(2)
        ! (f(x)-f(x-h)) / h
        fd = finite_diff_method(id,'2-point backward',   2,[0,-1],[1,-1],1)
    case(3)
        ! (f(x+h)-f(x-h)) / (2h)
        fd = finite_diff_method(id,'3-point central',    3,[1,-1],[1,-1],2)
    case(4)
        ! (-3f(x)+4f(x+h)-f(x+2h)) / (2h)
        fd = finite_diff_method(id,'3-point forward',    3,[0,1,2],[-3,4,-1],2)
    case(5)
        ! (f(x-2h)-4f(x-h)+3f(x)) / (2h)
        fd = finite_diff_method(id,'3-point backward',   3,[-2,-1,0],[1,-4,3],2)
    case(6)
        ! (-2f(x-h)-3f(x)+6f(x+h)-f(x+2h)) / (6h)
        fd = finite_diff_method(id,'4-point forward 1',  4,[-1,0,1,2],[-2,-3,6,-1],6)
    case(7)
        ! (f(x-2h)-6f(x-h)+3f(x)+2f(x+h)) / (6h)
        fd = finite_diff_method(id,'4-point backward 1', 4,[-2,-1,0,1],[1,-6,3,2],6)
    case(8)
        ! (-11f(x)+18f(x+h)-9f(x+2h)+2f(x+3h)) / (6h)
        fd = finite_diff_method(id,'4-point forward 2',  4,[0,1,2,3],[-11,18,-9,2],6)
    case(9)
        ! (-2f(x-3h)+9f(x-2h)-18f(x-h)+11f(x)) / (6h)
        fd = finite_diff_method(id,'4-point backward 2', 4,[-3,-2,-1,0],[-2,9,-18,11],6)
    case(10)
        ! (f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h)) / (12h)
        fd = finite_diff_method(id,'5-point central',    5,[-2,-1,1,2],[1,-8,8,-1],12)
    case(11)
        ! (-3f(x-h)-10f(x)+18f(x+h)-6f(x+2h)+f(x+3h)) / (12h)
        fd = finite_diff_method(id,'5-point forward 2',  5,[-1,0,1,2,3],[-3,-10,18,-6,1],12)
    case(12)
        ! (-f(x-3h)+6f(x-2h)-18f(x-h)+10f(x)+3f(x+h)) / (12h)
        fd = finite_diff_method(id,'5-point backward 2', 5,[-3,-2,-1,0,1],[-1,6,-18,10,3],12)
    case(13)
        ! (-25f(x)+48f(x+h)-36f(x+2h)+16f(x+3h)-3f(x+4h)) / (12h)
        fd = finite_diff_method(id,'5-point forward 1',  5,[0,1,2,3,4],[-25,48,-36,16,-3],12)
    case(14)
        ! (3f(x-4h)-16f(x-3h)+36f(x-2h)-48f(x-h)+25f(x)) / (12h)
        fd = finite_diff_method(id,'5-point backward 1', 5,[-4,-3,-2,-1,0],[3,-16,36,-48,25],12)
    case default
        found = .false.
    end select

    end subroutine get_finite_difference_method
!*******************************************************************************

!*******************************************************************************
!>
!  Returns all the methods with the given `class`.

    function get_all_methods_in_class(class) result(list_of_methods)

    implicit none

    integer,intent(in) :: class  !! the `class` is the number of points in the
                                 !! finite different computation
    type(meth_array) :: list_of_methods  !! all the methods with the specified `class`

    type(finite_diff_method) :: fd  !! temp variable for getting a
                                    !! method from [[get_finite_difference_method]]
    integer :: id     !! method id counter
    logical :: found  !! status flag
    integer :: n      !! temp size variable
    type(finite_diff_method),dimension(:),allocatable :: tmp  !! for resizing `meth` array

    ! currently, the only way to do this is to call the
    ! get_finite_difference_method routine and see if there
    ! is one available.
    id = 0
    do
        id = id + 1
        call get_finite_difference_method(id,fd,found)
        if (found) then
            if (fd%class==class) then
                if (allocated(list_of_methods%meth)) then
                    ! add to the list
                    n = size(list_of_methods%meth)
                    allocate(tmp(n+1))
                    tmp(1:n) = list_of_methods%meth
                    tmp(n+1) = fd
                    call move_alloc(tmp,list_of_methods%meth)
                    ! ... this doesn't appear to work on Intel compiler ...
                    !list_of_methods%meth = [list_of_methods%meth,fd]  ! add to the list
                else
                    ! first element:
                    allocate(list_of_methods%meth(1))
                    list_of_methods%meth = fd
                end if
            elseif (fd%class>class) then ! we assume they are in increasing order
                exit
            end if
        else
            exit ! done
        end if
    end do

    end function get_all_methods_in_class
!*******************************************************************************

!*******************************************************************************
!>
!  Select a finite diff method of a given `class` so that the perturbations
!  of `x` will not violate the variable bounds.

    subroutine select_finite_diff_method(me,x,xlow,xhigh,dx,list_of_methods,fd,status_ok)

    implicit none

    class(numdiff_type),intent(inout)    :: me
    real(wp),intent(in)                  :: x               !! the variable value
    real(wp),intent(in)                  :: xlow            !! the variable lower bound
    real(wp),intent(in)                  :: xhigh           !! the variable upper bound
    real(wp),intent(in)                  :: dx              !! the perturbation value (>0)
    type(meth_array),intent(in)          :: list_of_methods !! list of available methods to choose from
    type(finite_diff_method),intent(out) :: fd              !! this method can be used
    logical,intent(out)                  :: status_ok       !! true if it really doesn't violate the bounds
                                                            !! (say, the bounds are very close to each other)
                                                            !! if `status_ok=False`, then the first method in
                                                            !! the given class is returned in `fd`.

    integer  :: i   !! counter
    integer  :: j   !! counter
    real(wp) :: xp  !! perturbed `x` value

    ! initialize:
    status_ok = .false.

    ! try all the methods in the class:
    do i = 1, size(list_of_methods%meth)
        status_ok = .true. ! will be set to false if any
                           ! perturbation violates the bounds
        ! check each of the perturbations:
        do j = 1, size(list_of_methods%meth(i)%dx_factors)
            xp = x + list_of_methods%meth(i)%dx_factors(j)*dx
            if (xp < xlow .or. xp > xhigh) then
                status_ok = .false.
                exit
            end if
        end do
        if (status_ok) then   ! this one is OK to use
            fd = list_of_methods%meth(i)
            exit
        end if
    end do

    if (.not. status_ok) then
        ! no method was found that doesn't violate the bounds,
        ! so just return the first one in the list.
        fd = list_of_methods%meth(1)
    end if

    end subroutine select_finite_diff_method
!*******************************************************************************

!*******************************************************************************
!>
!  Select a finite diff method of a given `class` so that the perturbations
!  of `x` will not violate the variable bounds for any variable in the group.
!
!  The `x` vector are only the variables in a group (not the full variable vector)

    subroutine select_finite_diff_method_for_partition_group(me,x,xlow,xhigh,dx,&
                                                             list_of_methods,fd,status_ok)

    implicit none

    class(numdiff_type),intent(inout)        :: me
    real(wp),dimension(:),intent(in)         :: x               !! the variable values
    real(wp),dimension(:),intent(in)         :: xlow            !! the variable lower bounds
    real(wp),dimension(:),intent(in)         :: xhigh           !! the variable upper bounds
    real(wp),dimension(:),intent(in)         :: dx              !! the perturbation values (>0)
    type(meth_array),intent(in)              :: list_of_methods !! list of available methods to choose from
    type(finite_diff_method),intent(out)     :: fd              !! this method can be used
    logical,intent(out)                      :: status_ok       !! true if it really doesn't violate the bounds
                                                                !! (say, the bounds are very close to each other)
                                                                !! if `status_ok=False`, then the first method in
                                                                !! the given class is returned in `fd`.

    integer  :: i   !! counter
    integer  :: j   !! counter
    real(wp),dimension(size(x)) :: xp  !! perturbed `x` values

    ! initialize:
    status_ok = .false.

    ! try all the methods in the class:
    do i = 1, size(list_of_methods%meth)
        status_ok = .true. ! will be set to false if any
                           ! perturbation violates the bounds
        ! check each of the perturbations:
        do j = 1, size(list_of_methods%meth(i)%dx_factors)
            xp = x + list_of_methods%meth(i)%dx_factors(j)*dx
            if (any(xp < xlow) .or. any(xp > xhigh)) then
                status_ok = .false.
                exit
            end if
        end do
        if (status_ok) then   ! this one is OK to use
            fd = list_of_methods%meth(i)
            exit
        end if
    end do

    if (.not. status_ok) then
        ! no method was found that doesn't violate the bounds,
        ! so just return the first one in the list.
        fd = list_of_methods%meth(1)
    end if

    end subroutine select_finite_diff_method_for_partition_group
!*******************************************************************************

!*******************************************************************************
!>
!  Alternate version of [[initialize_numdiff]] routine when
!  using [[diff]] to compute the Jacobian.

    subroutine initialize_numdiff_for_diff(me,n,m,xlow,xhigh,&
                                    problem_func,sparsity_mode,info,&
                                    chunk_size,eps,acc,cache_size)

    implicit none

    class(numdiff_type),intent(out)  :: me
    integer,intent(in)               :: n             !! number of `x` variables
    integer,intent(in)               :: m             !! number of `f` functions
    real(wp),dimension(n),intent(in) :: xlow          !! lower bounds on `x`
    real(wp),dimension(n),intent(in) :: xhigh         !! upper bounds on `x`
    procedure(func)                  :: problem_func  !! the user function that defines the problem
                                                      !! (returns `m` functions)
    integer,intent(in)               :: sparsity_mode !! the sparsity computation method:
                                                      !! **1** - assume dense,
                                                      !! **2** - three-point method,
                                                      !! **3** - will be specified by the user in
                                                      !! a subsequent call to [[set_sparsity_pattern]].
    procedure(info_f),optional       :: info          !! a function the user can define
                                                      !! which is called when each column
                                                      !! of the jacobian is computed.
                                                      !! It can be used to perform any
                                                      !! setup operations.
    integer,intent(in),optional      :: chunk_size    !! chunk size for allocating the arrays
                                                      !! (must be >0) [default is 100]
    real(wp),intent(in),optional     :: eps           !! tolerance parameter for [[diff]]
                                                      !! if not present, default is `1.0e-9_wp`
    real(wp),intent(in),optional     :: acc           !! tolerance parameter for [[diff]]
                                                      !! if not present, default is `0.0_wp`
    integer,intent(in),optional      :: cache_size    !! if present, this is the cache size
                                                      !! for the function cache
                                                      !! (default is not to enable cache)

    logical :: cache  !! if the cache is to be used

    me%use_diff = .true.

    ! cache:
    cache = present(cache_size)
    if (cache) cache = cache_size>0
    if (cache) then
        me%compute_function  => compute_function_with_cache
        call me%cache%initialize(cache_size,n,m)
    else
        me%compute_function  => problem_func
        call me%cache%destroy()
    end if

    ! functions:
    me%problem_func      => problem_func
    me%jacobian_function => compute_jacobian_with_diff

    ! set the sparsity function
    select case (sparsity_mode)
    case(1)  ! dense
        me%compute_sparsity => compute_sparsity_dense
    case(2)  ! three-point method
        me%compute_sparsity => compute_sparsity_random
    case(3)  ! user defined
        me%compute_sparsity => null()
    case default
        error stop 'Error: sparsity_mode must be 1, 2, or 3.'
    end select

    ! size of the problem:
    me%n = n
    me%m = m

    ! input variable bounds:
    if (any(xlow>=xhigh)) error stop 'Error: all xlow must be < xhigh'
    allocate(me%xlow(n))
    allocate(me%xhigh(n))
    me%xlow = xlow
    me%xhigh = xhigh

    ! optional:
    if (present(chunk_size)) me%chunk_size = abs(chunk_size)
    if (present(eps))        me%eps = eps
    if (present(acc))        me%acc = acc
    if (present(info))       me%info_function => info

    end subroutine initialize_numdiff_for_diff
!*******************************************************************************

!*******************************************************************************
!>
!  Initialize a [[numdiff_type]] class. This must be called first.
!
!@note Only one of the following inputs can be used: `jacobian_method`,
!      `jacobian_methods`, `class`, or `classes`.

    subroutine initialize_numdiff(me,n,m,xlow,xhigh,perturb_mode,dpert,&
                        problem_func,sparsity_mode,jacobian_method,jacobian_methods,&
                        class,classes,info,chunk_size,partition_sparsity_pattern,&
                        cache_size)

    implicit none

    class(numdiff_type),intent(out)     :: me
    integer,intent(in)                  :: n               !! number of `x` variables
    integer,intent(in)                  :: m               !! number of `f` functions
    real(wp),dimension(n),intent(in)    :: xlow            !! lower bounds on `x`
    real(wp),dimension(n),intent(in)    :: xhigh           !! upper bounds on `x`
    integer,intent(in)                  :: perturb_mode    !! perturbation mode:
                                                           !! **1** - perturbation is `dx=dpert`,
                                                           !! **2** - perturbation is `dx=dpert*x`,
                                                           !! **3** - perturbation is `dx=dpert*(1+x)`
    real(wp),dimension(n),intent(in)    :: dpert           !! perturbation vector for `x`
    procedure(func)                     :: problem_func    !! the user function that defines the problem
                                                           !! (returns `m` functions)
    integer,intent(in)                  :: sparsity_mode   !! the sparsity computation method:
                                                           !! **1** - assume dense,
                                                           !! **2** - three-point method,
                                                           !! **3** - will be specified by the user in
                                                           !! a subsequent call to [[set_sparsity_pattern]].
    integer,intent(in),optional         :: jacobian_method !! `id` code for the finite difference method
                                                           !! to use for all `n` variables.
                                                           !! see [[get_finite_difference_method]]
    integer,dimension(n),intent(in),optional :: jacobian_methods !! `id` codes for the finite difference method
                                                                 !! to use for each variable.
                                                                 !! see [[get_finite_difference_method]]
    integer,intent(in),optional :: class  !! method class for the finite difference method
                                          !! to use for all `n` variables.
                                          !! see [[get_finite_difference_method]]
    integer,dimension(n),intent(in),optional :: classes   !! method class for the finite difference methods
                                                          !! to use for each variable.
                                                          !! see [[get_finite_difference_method]]
    procedure(info_f),optional :: info  !! a function the user can define
                                        !! which is called when each column
                                        !! of the jacobian is computed.
                                        !! It can be used to perform any
                                        !! setup operations.
    integer,intent(in),optional :: chunk_size  !! chunk size for allocating the arrays
                                               !! (must be >0) [default is 100]
    logical,intent(in),optional :: partition_sparsity_pattern  !! if the sparisty pattern is to
                                                               !! be partitioned using [[DSM]]
                                                               !! [default is False]
    integer,intent(in),optional      :: cache_size    !! if present, this is the cache size
                                                      !! for the function cache
                                                      !! (default is not to enable cache)

    integer :: i !! counter
    logical :: found
    logical :: cache  !! if the cache is to be used

    me%use_diff = .false.

    ! cache:
    cache = present(cache_size)
    if (cache) cache = cache_size>0
    if (cache) then
        me%compute_function  => compute_function_with_cache
        call me%cache%initialize(cache_size,n,m)
    else
        me%compute_function  => problem_func
        call me%cache%destroy()
    end if

    ! functions:
    me%problem_func => problem_func

    if (present(partition_sparsity_pattern)) then
        me%partition_sparsity_pattern = partition_sparsity_pattern
    else
        me%partition_sparsity_pattern = .false.
    end if

    ! method:
    if (allocated(me%meth)) deallocate(me%meth)
    if (allocated(me%class)) deallocate(me%class)
    if (allocated(me%class_meths)) deallocate(me%class_meths)

    if (      present(jacobian_method) .and. .not. present(jacobian_methods) .and. &
        .not. present(class) .and. .not. present(classes)) then
        ! use the same for all variable
        me%mode = 1
        allocate(me%meth(n))
        do i=1,n
            call get_finite_difference_method(jacobian_method,me%meth(i),found)
            if (.not. found) error stop 'Error: invalid jacobian_method'
        end do
    elseif (.not. present(jacobian_method) .and. present(jacobian_methods) .and. &
            .not. present(class) .and. .not. present(classes)) then
        ! specify a separate method for each variable
        me%mode = 1
        allocate(me%meth(n))
        do i=1,n
            call get_finite_difference_method(jacobian_methods(i),me%meth(i),found)
            if (.not. found) error stop 'Error: invalid jacobian_methods'
        end do
        if (me%partition_sparsity_pattern) error stop 'Error: when using partitioned '//&
            'sparsity pattern, all columns must use the same finite diff method.'
    elseif (.not. present(jacobian_method) .and. .not. present(jacobian_methods) .and. &
                  present(class) .and. .not. present(classes)) then
        ! use the class for all variables
        me%mode = 2
        allocate(me%class(n))
        me%class = class
        allocate(me%class_meths(n))
        me%class_meths(1) = get_all_methods_in_class(class)
        if (n>1) me%class_meths(2:n) = me%class_meths(1)  ! just copy them over
    elseif (.not. present(jacobian_method) .and. .not. present(jacobian_methods) .and. &
            .not. present(class) .and. present(classes)) then
        ! specify a separate class for each variable
        me%mode = 2
        me%class = classes
        allocate(me%class_meths(n))
        do i=1,n
            me%class_meths(i) = get_all_methods_in_class(me%class(i))
        end do
        if (me%partition_sparsity_pattern) error stop 'Error: when using partitioned '//&
            'sparsity pattern, all columns must use the same finite diff method.'
    else
        error stop 'Error: must specify one of either jacobian_method, jacobian_methods, class, or classes.'
    end if

    ! size of the problem:
    me%n = n
    me%m = m

    ! input variable bounds:
    if (any(xlow>=xhigh)) error stop 'Error: all xlow must be < xhigh'
    if (allocated(me%xlow)) deallocate(me%xlow)
    if (allocated(me%xhigh)) deallocate(me%xhigh)
    allocate(me%xlow(n))
    allocate(me%xhigh(n))
    me%xlow = xlow
    me%xhigh = xhigh

    ! perturbation options:
    select case (perturb_mode)
    case(1:3)
        me%perturb_mode = perturb_mode
    case default
        error stop 'Error: perturb_mode must be 1, 2, or 3.'
    end select
    if (allocated(me%dpert)) deallocate(me%dpert)
    allocate(me%dpert(n))
    me%dpert = abs(dpert)

    ! optional:
    if (present(info))       me%info_function => info
    if (present(chunk_size)) me%chunk_size = abs(chunk_size)

    ! set the jacobian function, depending on the options:
    if (me%partition_sparsity_pattern) then
        me%jacobian_function => compute_jacobian_partitioned
    else
        me%jacobian_function => compute_jacobian_standard
    end if

    ! set the sparsity function
    select case (sparsity_mode)
    case(1)  ! dense
        me%compute_sparsity => compute_sparsity_dense
    case(2)  ! three-point method
        me%compute_sparsity => compute_sparsity_random
    case(3)  ! user defined
        me%compute_sparsity => null()
    case default
        error stop 'Error: sparsity_mode must be 1, 2, or 3.'
    end select

    end subroutine initialize_numdiff
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
!  destroy a [[sparsity_pattern]] type.

    subroutine destroy_sparsity(me)

    implicit none

    class(sparsity_pattern),intent(out) :: me

    end subroutine destroy_sparsity
!*******************************************************************************

!*******************************************************************************
!>
!  Returns the columns in a sparsity partition group.
!
!@note This is just a wrapper to get data from `ngrp`.

    subroutine columns_in_partition_group(me,igroup,n_cols,cols,nonzero_rows,indices)

    implicit none

    class(sparsity_pattern),intent(in)           :: me
    integer,intent(in)                           :: igroup       !! group number. Should be `>0` and `<=me%mxgrp`
    integer,intent(out)                          :: n_cols       !! number of columns in the `igroup` group.
    integer,dimension(:),allocatable,intent(out) :: cols         !! the column numbers in the `igroup` group.
                                                                 !! (if none, then it is not allocated)
    integer,dimension(:),allocatable,intent(out) :: nonzero_rows !! the row numbers of all the nonzero
                                                                 !! Jacobian elements in this group
    integer,dimension(:),allocatable,intent(out) :: indices      !! nonzero indices in `jac` for a group

    integer :: i  !! counter
    integer :: num_nonzero_elements_in_col          !! number of nonzero elements in a column
    integer :: num_nonzero_elements_in_group        !! number of nonzero elements in a group
    integer,dimension(:),allocatable :: col_indices !! nonzero indices in `jac` for a column

    if (me%maxgrp>0 .and. allocated(me%ngrp)) then

        n_cols = count(me%ngrp==igroup)
        if (n_cols>0) then
            allocate(cols(n_cols))
            cols = pack([(i,i=1,size(me%ngrp))],mask=me%ngrp==igroup)
        end if

        ! get all the non-zero elements in each column:
        num_nonzero_elements_in_group = 0  ! initialize
        do i = 1, n_cols
            num_nonzero_elements_in_col = count(me%icol==cols(i))
            if (num_nonzero_elements_in_col/=0) then ! there are functions to
                                                     ! compute in this column
                num_nonzero_elements_in_group = num_nonzero_elements_in_group + &
                                                num_nonzero_elements_in_col
                if (allocated(col_indices)) deallocate(col_indices)
                allocate(col_indices(num_nonzero_elements_in_col))
                col_indices = pack(me%indices,mask=me%icol==cols(i))
                if (allocated(nonzero_rows)) then
                    nonzero_rows = [nonzero_rows,me%irow(col_indices)]
                    indices = [indices,col_indices]
                else
                    nonzero_rows = me%irow(col_indices)
                    indices = col_indices
                end if
            end if
        end do

    else
        error stop 'Error: the partition has not been computed.'
    end if

    end subroutine columns_in_partition_group
!*******************************************************************************

!*******************************************************************************
!>
!  Destroy the sparsity pattern in the class.

    subroutine destroy_sparsity_pattern(me)

    implicit none

    class(numdiff_type),intent(inout) :: me

    call me%sparsity%destroy()

    end subroutine destroy_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  To specify the sparsity pattern directly if it is already known.

    subroutine set_sparsity_pattern(me,irow,icol)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,dimension(:),intent(in)   :: irow  !! sparsity pattern nonzero elements row indices
    integer,dimension(:),intent(in)   :: icol  !! sparsity pattern nonzero elements column indices

    integer                   :: Mingrp !! for call to [[dsm]]
    integer                   :: Info   !! for call to [[dsm]]
    integer,dimension(me%m+1) :: ipntr  !! for call to [[dsm]]
    integer,dimension(me%n+1) :: jpntr  !! for call to [[dsm]]
    integer                   :: i      !! counter

    call me%destroy_sparsity_pattern()

    if (size(irow)/=size(icol) .or. any(irow>me%m) .or. any(icol>me%n)) then
        error stop 'Error: invalid inputs to set_sparsity_pattern'
    else

        me%sparsity%sparsity_computed = .true.
        me%sparsity%num_nonzero_elements = size(irow)
        me%sparsity%irow = irow
        me%sparsity%icol = icol

        allocate(me%sparsity%indices(me%sparsity%num_nonzero_elements))
        me%sparsity%indices = [(i,i=1,me%sparsity%num_nonzero_elements)]

        if (me%partition_sparsity_pattern) then
            associate( s => me%sparsity )
                allocate(s%ngrp(me%n))
                call dsm(me%m,me%n,s%num_nonzero_elements,&
                         s%irow,s%icol,&
                         s%ngrp,s%maxgrp,&
                         mingrp,info,ipntr,jpntr)
                if (info/=1) error stop 'Error partitioning sparsity pattern.'
            end associate
        end if

    end if

    end subroutine set_sparsity_pattern
!*******************************************************************************

!*******************************************************************************
!>
!  assume all elements of Jacobian are non-zero.

    subroutine compute_sparsity_dense(me,x)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x   !! vector of variables (size `n`)

    integer :: i !! counter
    integer :: r !! row counter
    integer :: c !! column counter

    call me%destroy_sparsity_pattern()

    me%sparsity%num_nonzero_elements = me%m * me%n
    allocate(me%sparsity%irow(me%sparsity%num_nonzero_elements))
    allocate(me%sparsity%icol(me%sparsity%num_nonzero_elements))

    allocate(me%sparsity%indices(me%sparsity%num_nonzero_elements))
    me%sparsity%indices = [(i,i=1,me%sparsity%num_nonzero_elements)]

    ! create the dense matrix:
    i = 0
    do c = 1, me%n
        do r = 1, me%m
            i = i + 1
            me%sparsity%irow(i) = r
            me%sparsity%icol(i) = c
        end do
    end do

    ! ... no real need for this, since it can't  ...
    ! ... be partitioned (all elements are true) ...
    if (me%partition_sparsity_pattern) then
        ! generate a "dense" partition
        me%sparsity%maxgrp = me%n
        allocate(me%sparsity%ngrp(me%n))
        me%sparsity%ngrp = [(i, i=1,me%n)]
    end if

    me%sparsity%sparsity_computed = .true.

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
!
!@note Instead of three points, the number of points could be a user input.

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

    integer                   :: Mingrp  !! for call to [[dsm]]
    integer                   :: Info    !! for call to [[dsm]]
    integer,dimension(me%m+1) :: ipntr   !! for call to [[dsm]]
    integer,dimension(me%n+1) :: jpntr   !! for call to [[dsm]]

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
    call me%destroy_sparsity_pattern()

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
                call expand_vector(me%sparsity%icol,n_icol,me%chunk_size,val=i)
                call expand_vector(me%sparsity%irow,n_irow,me%chunk_size,val=j)
            end if
        end do
        ! resize to correct size:
        call expand_vector(me%sparsity%icol,n_icol,me%chunk_size,finished=.true.)
        call expand_vector(me%sparsity%irow,n_irow,me%chunk_size,finished=.true.)

    end do

    me%sparsity%num_nonzero_elements = size(me%sparsity%irow)

    allocate(me%sparsity%indices(me%sparsity%num_nonzero_elements))
    me%sparsity%indices = [(i,i=1,me%sparsity%num_nonzero_elements)]

    if (me%partition_sparsity_pattern) then
        associate( s => me%sparsity )
            allocate(s%ngrp(me%n))
            call dsm(me%m,me%n,s%num_nonzero_elements,&
                     s%irow,s%icol,&
                     s%ngrp,s%maxgrp,&
                     mingrp,info,ipntr,jpntr)
            if (info/=1) error stop 'Error partitioning sparsity pattern.'
        end associate
    end if

    ! finished:
    me%sparsity%sparsity_computed = .true.

    end subroutine compute_sparsity_random
!*******************************************************************************

!*******************************************************************************
!>
!  just a wrapper for [[compute_jacobian]], that returns a dense (`m x n`) matrix.

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
    do i = 1, me%sparsity%num_nonzero_elements
        jac(me%sparsity%irow(i),me%sparsity%icol(i)) = jac_vec(i)
    end do

    end subroutine compute_jacobian_dense
!*******************************************************************************

!*******************************************************************************
!>
!  Perturb the specified optimization variable, and compute the function.
!  This routine is designed so that `df` is accumulated as each function
!  evaluation is done, to avoid having to allocate more temporary storage.

    subroutine perturb_x_and_compute_f(me,x,dx_factor,dx,&
                                       df_factor,column,idx,df)

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

    real(wp),dimension(me%n) :: xp  !! the perturbed variable vector
    real(wp),dimension(me%m) :: f   !! function evaluation

    xp = x
    if (dx_factor/=zero) xp(column) = xp(column) + dx_factor * dx(column)
    call me%compute_function(xp,f,idx)
    df(idx) = df(idx) + df_factor * f(idx)

    end subroutine perturb_x_and_compute_f
!*******************************************************************************

!*******************************************************************************
!>
!  Returns the product `J*v`, where `J` is the `m x n` Jacobian matrix
!  and `v` is an `n x 1` vector.

    subroutine compute_jacobian_times_vector(me,x,v,z)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x    !! vector of variables (size `n`)
    real(wp),dimension(:),intent(in)  :: v    !! vector (size `n`)
    real(wp),dimension(:),intent(out) :: z    !! The product `J*v` (size `m`)

    real(wp),dimension(:),allocatable :: jac !! sparse jacobian vector
    integer :: i !! counter
    integer :: r !! row number in full jacobian
    integer :: c !! column number in full jacobian

    ! first compute the jacobian in sparse vector form:
    call me%compute_jacobian(x,jac)

    ! initialize output vector:
    z = 0.0_wp

    if (allocated(jac)) then

        !   J    v   z
        !  ---   -   -
        !  X0X   X   X
        !  0X0 * X = X
        !  00X   X   X
        !  00X       X

        ! multiplication by input v:
        do i = 1, me%sparsity%num_nonzero_elements
            r = me%sparsity%irow(i)
            c = me%sparsity%icol(i)
            z(r) = z(r) + jac(i)*v(c)
        end do

        deallocate(jac)

    end if

    end subroutine compute_jacobian_times_vector
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian.

    subroutine compute_jacobian(me,x,jac)

    implicit none

    class(numdiff_type),intent(inout)             :: me
    real(wp),dimension(:),intent(in)              :: x    !! vector of variables (size `n`)
    real(wp),dimension(:),allocatable,intent(out) :: jac  !! sparse jacobian vector

    real(wp),dimension(me%n) :: dx  !! absolute perturbation (>0) for each variable

    ! if we don't have a sparsity pattern yet then compute it:
    ! [also computes the indices vector]
    if (.not. me%sparsity%sparsity_computed) call me%compute_sparsity(x)

    ! size the jacobian vector:
    allocate(jac(me%sparsity%num_nonzero_elements))

    ! compute dx vector:
    !if (.not. associated(me%jacobian_function,compute_jacobian_with_diff)) then  ! this doesn't work with Intel compiler ...
    if (.not. me%use_diff) then
        ! only need this for the finite difference methods (not diff)
        call me%compute_perturbation_vector(x,dx)
    end if

    ! compute the jacobian:
    if (associated(me%jacobian_function)) then
        call me%jacobian_function(x,dx,jac)
    else
        error stop 'Error: jacobian_function has not been associated.'
    end if

    end subroutine compute_jacobian
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian using finite differences.
!  (one column at a time)

    subroutine compute_jacobian_standard(me,x,dx,jac)

    implicit none

    class(numdiff_type),intent(inout)   :: me
    real(wp),dimension(:),intent(in)    :: x    !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(in) :: dx   !! absolute perturbation (>0)
                                                !! for each variable
    real(wp),dimension(:),intent(out)   :: jac  !! sparse jacobian vector (size
                                                !! `num_nonzero_elements`)

    integer,dimension(:),allocatable :: nonzero_elements_in_col  !! the indices of the
                                                                 !! nonzero Jacobian
                                                                 !! elements in a column
    integer                   :: i   !! column counter
    integer                   :: j   !! function evaluation counter
    real(wp),dimension(me%m)  :: df  !! accumulated function
    type(finite_diff_method)  :: fd  !! a finite different method (when
                                     !! specifying class rather than the method)
    logical :: status_ok   !! error flag
    integer :: num_nonzero_elements_in_col  !! number of nonzero elements in a column

    ! initialize:
    jac = zero

    ! compute Jacobian matrix column-by-column:
    do i=1,me%n

        ! determine functions to compute for this column:
        num_nonzero_elements_in_col = count(me%sparsity%icol==i)
        if (num_nonzero_elements_in_col/=0) then ! there are functions to compute

            nonzero_elements_in_col = pack(me%sparsity%irow,mask=me%sparsity%icol==i)

            select case (me%mode)
            case(1) ! use the specified methods

                ! compute this column of the Jacobian:
                df = zero
                do j = 1, size(me%meth(i)%dx_factors)
                    if (associated(me%info_function)) call me%info_function([i],j)
                    call me%perturb_x_and_compute_f(x,me%meth(i)%dx_factors(j),&
                                                    dx,me%meth(i)%df_factors(j),&
                                                    i,nonzero_elements_in_col,df)
                end do

                df(nonzero_elements_in_col) = df(nonzero_elements_in_col) / &
                                              (me%meth(i)%df_den_factor*dx(i))

            case(2) ! select the method from the class so as not to violate the bounds

                call me%select_finite_diff_method(x(i),me%xlow(i),me%xhigh(i),&
                                                  dx(i),me%class_meths(i),fd,status_ok)
                if (.not. status_ok) write(error_unit,'(A,1X,I5)') &
                    'Error: variable bounds violated for column: ',i

                ! compute this column of the Jacobian:
                df = zero
                do j = 1, size(fd%dx_factors)
                    if (associated(me%info_function)) call me%info_function([i],j)
                    call me%perturb_x_and_compute_f(x,fd%dx_factors(j),&
                                                    dx,fd%df_factors(j),&
                                                    i,nonzero_elements_in_col,df)
                end do
                df(nonzero_elements_in_col) = df(nonzero_elements_in_col) / &
                                              (fd%df_den_factor*dx(i))

            case default
                error stop 'Error: invalid mode'
            end select

            ! put result into the output vector:
            jac(pack(me%sparsity%indices,mask=me%sparsity%icol==i)) = &
                df(nonzero_elements_in_col)

        end if

    end do

    end subroutine compute_jacobian_standard
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian one element at a time using the Neville's process
!  algorithm [[diff]]. This takes a very large number of function evaluations,
!  but should give a very accurate answer.

    subroutine compute_jacobian_with_diff(me,x,dx,jac)

    implicit none

    class(numdiff_type),intent(inout)   :: me
    real(wp),dimension(:),intent(in)    :: x    !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(in) :: dx   !! absolute perturbation (>0)
                                                !! for each variable
    real(wp),dimension(:),intent(out)   :: jac  !! sparse jacobian vector (size
                                                !! `num_nonzero_elements`)

    integer,parameter  :: iord  = 1  !! tells [[diff]] to compute first derivative

    integer                  :: i      !! counter for nonzero elements of jacobian
    type(diff_func)          :: d      !! the [[diff]] class to use
    real(wp)                 :: x0     !! value of ith variable for [[diff]]
    real(wp)                 :: xmin   !! ith variable lower bound
    real(wp)                 :: xmax   !! ith variable upper bound
    real(wp)                 :: deriv  !! derivative value df/dx from [[diff]]
    real(wp)                 :: error  !! estimated error from [[diff]]
    integer                  :: ifail  !! [[diff]] error flag
    real(wp),dimension(me%n) :: xp     !! used by [[dfunc]].
    real(wp),dimension(me%m) :: fvec   !! used by [[dfunc]].
    integer                  :: ir     !! row index of next non-zero element
    integer                  :: ic     !! column index of next non-zero element

    integer :: ic_prev  !! previous column perturbed
    integer :: icount   !! count of number of times a column has been perturbed
    logical :: use_info !! if we are reporting to the user

    ! initialize:
    jac = zero
    use_info = associated(me%info_function)
    if (use_info) ic_prev = -1

    ! set the function for diff:
    call d%set_function(dfunc)

    ! each element is computed one by one
    do i = 1, me%sparsity%num_nonzero_elements

        ir   = me%sparsity%irow(i)
        ic   = me%sparsity%icol(i)
        x0   = x(ic)
        xmin = me%xlow(ic)
        xmax = me%xhigh(ic)

        ! if reporting to the user, have to keep track
        ! of which column is being perturbed, and how
        ! many times:
        if (use_info) then
            if (ic/=ic_prev) icount = 0
        end if

        call d%compute_derivative(iord,x0,xmin,xmax,me%eps,me%acc,deriv,error,ifail)

        if (ifail==0 .or. ifail==1) then
            jac(i) = deriv
        else
            error stop 'Error computing derivative with DIFF.'
        end if

        if (use_info) ic_prev = ic

    end do

    contains

        function dfunc(this,xval) result(fx)

        !! interface to user function for [[diff]] routine.

        implicit none

        class(diff_func),intent(inout)  :: this
        real(wp),intent(in)             :: xval  !! input variable (`ic` variable)
        real(wp)                        :: fx    !! derivative of `ir` function
                                                 !! w.r.t. `xval` variable

        if (use_info) then
            icount = icount + 1
            call me%info_function([ic],icount)
        end if

        xp = x
        xp(ic) = xval
        call me%compute_function(xp,fvec,funcs_to_compute=[ir])

        fx = fvec(ir)

        end function dfunc

    end subroutine compute_jacobian_with_diff
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian using finite differences,
!  (using the partitioned sparsity pattern to compute multiple columns
!  at a time).

    subroutine compute_jacobian_partitioned(me,x,dx,jac)

    implicit none

    class(numdiff_type),intent(inout)    :: me
    real(wp),dimension(:),intent(in)     :: x    !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(in)  :: dx   !! absolute perturbation (>0) for each variable
    real(wp),dimension(:),intent(out)    :: jac  !! sparse jacobian vector

    integer                          :: i             !! column counter
    integer                          :: j             !! function evaluation counter
    integer                          :: igroup        !! group number counter
    integer                          :: n_cols        !! number of columns in a group
    integer,dimension(:),allocatable :: cols          !! array of column indices in a group
    integer,dimension(:),allocatable :: nonzero_rows  !! the indices of the nonzero Jacobian
                                                      !! elementes (row numbers) in a group
    integer,dimension(:),allocatable :: indices       !! nonzero indices in `jac` for a group
    integer,dimension(:),allocatable :: col_indices   !! nonzero indices in `jac` for a column
    real(wp),dimension(me%m)         :: df            !! accumulated function
    type(finite_diff_method)         :: fd            !! a finite different method (when
                                                      !! specifying class rather than the method)
    logical                          :: status_ok     !! error flag

    integer :: num_nonzero_elements_in_col

    ! initialize:
    jac = zero

    ! compute by group:
    do igroup = 1, me%sparsity%maxgrp

        ! get the columns in this group:
        call me%sparsity%columns_in_partition_group(igroup,n_cols,cols,&
                                                    nonzero_rows,indices)
        if (n_cols>0) then

            if (allocated(nonzero_rows)) then

                select case (me%mode)
                case(1) ! use the specified methods

                    ! note: all the methods must be the same within a group

                    ! compute the columns of the Jacobian in this group:
                    df = zero
                    do j = 1, size(me%meth(1)%dx_factors)
                         if (associated(me%info_function)) call me%info_function(cols,j)
                         call me%perturb_x_and_compute_f_partitioned(x,me%meth(1)%dx_factors(j),&
                                                         dx,me%meth(1)%df_factors(j),&
                                                         cols,nonzero_rows,df)
                    end do
                    ! divide by the denominator, which can be different for each column:
                    do i = 1, n_cols
                        num_nonzero_elements_in_col = count(me%sparsity%icol==cols(i))
                        if (allocated(col_indices)) deallocate(col_indices)
                        allocate(col_indices(num_nonzero_elements_in_col))
                        col_indices = pack(me%sparsity%indices,mask=me%sparsity%icol==cols(i))
                        df(me%sparsity%irow(col_indices)) = df(me%sparsity%irow(col_indices)) / &
                                                            (me%meth(1)%df_den_factor*dx(cols(i)))
                    end do

                case(2) ! select the method from the class so as not to violate
                        ! the bounds on *any* of the variables in the group

                    ! note: all the classes must be the same within a group

                    call me%select_finite_diff_method_for_partition_group( &
                                    x(cols),me%xlow(cols),me%xhigh(cols),&
                                    dx(cols),me%class_meths(1),fd,status_ok)

                    if (.not. status_ok) then
                        write(error_unit,'(A,1X,I5,1X,A,1X,*(I5,1X))') &
                              'Error: variable bounds violated for group: ',&
                              igroup,'. columns: ',cols
                    end if

                    ! compute the columns of the Jacobian in this group:
                    df = zero
                    do j = 1, size(fd%dx_factors)
                        if (associated(me%info_function)) call me%info_function(cols,j)
                        call me%perturb_x_and_compute_f_partitioned(x,fd%dx_factors(j),&
                                                        dx,fd%df_factors(j),&
                                                        cols,nonzero_rows,df)
                    end do
                    ! divide by the denominator, which can be different for each column:
                    do i = 1, n_cols
                        num_nonzero_elements_in_col = count(me%sparsity%icol==cols(i))
                        if (allocated(col_indices)) deallocate(col_indices)
                        allocate(col_indices(num_nonzero_elements_in_col))
                        col_indices = pack(me%sparsity%indices,mask=me%sparsity%icol==cols(i))
                        df(me%sparsity%irow(col_indices)) = df(me%sparsity%irow(col_indices)) / &
                                                            (fd%df_den_factor*dx(cols(i)))
                    end do

                case default
                    error stop 'Error: invalid mode'
                end select

                ! put result into the output vector:
                jac(indices) = df(nonzero_rows)

            end if

        end if

    end do

    end subroutine compute_jacobian_partitioned
!*******************************************************************************

!*******************************************************************************
!>
!  Perturb the specified optimization variable, and compute the function.
!  This routine is designed so that `df` is accumulated as each function
!  evaluation is done, to avoid having to allocate more temporary storage.

    subroutine perturb_x_and_compute_f_partitioned(me,x,dx_factor,dx,&
                                       df_factor,columns,idx,df)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x         !! nominal variable vector
    real(wp),intent(in)               :: dx_factor !! factor to multiply `dx`
    real(wp),dimension(:),intent(in)  :: dx        !! the perturbation value for this column
    real(wp),intent(in)               :: df_factor !! factor to multiply function value
    integer,dimension(:),intent(in)   :: columns   !! the variables to perturb
    integer,dimension(:),intent(in)   :: idx       !! the elements in this
                                                   !! column of the Jacobian
                                                   !! to compute (passed to function)
    real(wp),dimension(me%m),intent(inout) :: df   !! the accumulated function value
                                                   !! note: for the first call, this
                                                   !! should be set to zero

    real(wp),dimension(me%n) :: xp  !! the perturbed variable vector
    real(wp),dimension(me%m) :: f   !! function evaluation

    xp = x
    if (dx_factor/=zero) xp(columns) = xp(columns) + dx_factor * dx(columns)
    call me%compute_function(xp,f,idx)
    df(idx) = df(idx) + df_factor * f(idx)

    end subroutine perturb_x_and_compute_f_partitioned
!*******************************************************************************

!*******************************************************************************
!>
!  Compute `dx`, the perturbation vector for `x` used
!  when computing the gradients.

    subroutine compute_perturbation_vector(me,x,dx)

    implicit none

    class(numdiff_type),intent(inout)    :: me
    real(wp),dimension(me%n),intent(in)  :: x  !! vector of variables (size `n`)
    real(wp),dimension(me%n),intent(out) :: dx !! absolute perturbation (>0)
                                               !! for each variable

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
!  Print the sparsity pattern.

    subroutine print_sparsity(me,n,m,iunit,dense)

    implicit none

    class(sparsity_pattern),intent(inout) :: me
    integer,intent(in) :: n  !! number of variables (columns of jacobian)
    integer,intent(in) :: m  !! number of functions (rows of jacobian)
    integer,intent(in) :: iunit !! file unit to write to.
                                !! (assumed to be already opened)
    logical,intent(in),optional :: dense  !! if present and true, the matrix form
                                          !! of the sparsity pattern is printed
                                          !! (default is vector form)

    logical :: print_matrix !! if the matrix form is to be printed
    integer :: r   !! row counter
    character(len=1),dimension(n) :: row  !! a row of the sparsity matrix

    if (allocated(me%irow) .and. allocated(me%icol)) then

        if (present(dense)) then
            print_matrix = dense
        else
            print_matrix = .false.  ! default
        end if

        write(iunit,'(A)') '---Sparsity pattern---'
        if (print_matrix) then
            do r = 1,m    ! print by row
                row = '0'
                row(pack(me%icol,mask=me%irow==r)) = 'X'
                write(iunit,'(*(A1))') row
            end do
        else
            write(iunit,'(A,1X,*(I3,","))') 'irow:',me%irow
            write(iunit,'(A,1X,*(I3,","))') 'icol:',me%icol
        end if
        write(iunit,'(A)') ''

        if (allocated(me%ngrp)) then
            write(iunit,'(A)') '---Sparsity partition---'
            write(iunit,'(A,1x,I5)')       'Number of groups:',me%maxgrp
            write(iunit,'(A,1x,*(I5,1X))') 'Group array:     ',me%ngrp
            write(iunit,'(A)') ''
        end if

    else
        error stop 'Error: sparsity pattern not available.'
    end if

    end subroutine print_sparsity
!*******************************************************************************

!*******************************************************************************
!>
!  Print the sparsity pattern in vector form (`irow`, `icol`).

    subroutine print_sparsity_pattern(me,iunit)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in) :: iunit !! file unit to write to.
                                !! (assumed to be already opened)

    call me%sparsity%print(me%n,me%m,iunit,dense=.false.)

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

    call me%sparsity%print(me%n,me%m,iunit,dense=.true.)

    end subroutine print_sparsity_matrix
!*******************************************************************************

!*******************************************************************************
    end module numerical_differentiation_module
!*******************************************************************************
