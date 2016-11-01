!*******************************************************************************
!> author: Jacob Williams
!  date: October 29, 2016
!
!  Numerical differentiation module.

    module numerical_differentiation_module

    use iso_fortran_env, only: real64

    implicit none

    private

    integer,parameter,public :: wp = real64  !! default real kind

    type,abstract,public :: numdiff_type

        !! base type for sparsity and Jacobian computations.

        private

        integer :: n = 0 !! number of `x` variables
        integer :: m = 0 !! number of `f` functions

        real(wp),dimension(:),allocatable :: xlow  !! lower bounds on `x`
        real(wp),dimension(:),allocatable :: xhigh !! upper bounds on `x`

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

        procedure(info_f),pointer :: info_function => null()  !! a function the user can define
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

        ! these are required to be defined when the class is defined:
        procedure(func),deferred,public    :: compute_function  !! the user-defined function
        procedure(spars_f),deferred,public :: compute_sparsity  !! for computing the sparsity pattern
        procedure(jac_f),deferred,public   :: jacobian_method   !! for computing the Jacobian matrix

        ! internal routines:
        procedure :: destroy_sparsity            !! destroy the sparsity pattern
        procedure :: compute_perturbation_vector !! computes the variable perturbation factor
        procedure :: print_sparsity_pattern      !! print the sparsity pattern to a file
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
        subroutine jac_f(me,x,dx,column,indices_to_compute,dfdx)
            !! The function to compute a column of jacobian.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout)  :: me
            real(wp),dimension(:),intent(in)   :: x       !! vector of variables (size `n`)
            real(wp),dimension(:),intent(in)   :: dx      !! absolute perturbation (>0)
                                                          !! for each variable
            integer,intent(in)                 :: column  !! column number to compute
            integer,dimension(:),intent(in)    :: indices_to_compute !! the indices in the
                                                                     !! column of the Jacobian
                                                                     !! to compute
            real(wp),dimension(:),intent(out)  :: dfdx  !! the non-zero elements
                                                        !! in the Jacobian column
                                                        !! length is `size(indices_to_compute)`
        end subroutine jac_f
        subroutine info_f(me,column)
            !! User-defined info function (optional)
            import :: numdiff_type
            implicit none
            class(numdiff_type),intent(inout)  :: me
            integer,intent(in) :: column  !! the column about to be computed
        end subroutine info_f
    end interface

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the class. This must be called first.

    subroutine initialize_numdiff_type(me,n,m,xlow,xhigh,perturb_mode,dpert,info)

    implicit none

    class(numdiff_type),intent(inout)   :: me
    integer,intent(in)                  :: n            !! number of `x` variables
    integer,intent(in)                  :: m            !! number of `f` functions
    real(wp),dimension(n),intent(in)    :: xlow         !! lower bounds on `x`
    real(wp),dimension(n),intent(in)    :: xhigh        !! upper bounds on `x`
    integer,intent(in)                  :: perturb_mode !! perturbation mode (1,2,3)
    real(wp),dimension(n),intent(in)    :: dpert        !! perturbation vector for `x`
    procedure(info_f),optional          :: info         !! a function the user can define
                                                        !! which is called when each column of the jacobian is computed.
                                                        !! It can be used to perform any setup operations.

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
    if (present(info)) then
        me%info_function => info
    else
        me%info_function => null()
    end if

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
!  just a wrapper for `compute_jacobian`, that return a dense matrix.

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
    jac = 0.0_wp
    do i = 1, me%num_nonzero_elements
        jac(me%irow(i),me%icol(i)) = jac_vec(i)
    end do

    end subroutine compute_jacobian_dense
!*******************************************************************************

!*******************************************************************************
!>
!  Compute a column of the Jacobian matrix using basic forward differences.

    subroutine forward_diff(me,x,dx,column,idx,dfdx)

    implicit none

    class(numdiff_type),intent(inout)  :: me
    real(wp),dimension(:),intent(in)   :: x           !! vector of variables (size `n`)
    real(wp),dimension(:),intent(in)   :: dx          !! absolute perturbation (>0)
                                                      !! for each variable
    integer,intent(in)                 :: column      !! column number to compute
    integer,dimension(:),intent(in)    :: idx         !! the elements in the
                                                      !! column of the Jacobian
                                                      !! to compute
    real(wp),dimension(size(idx)),intent(out) :: dfdx !! the non-zero elements in the
                                                      !! column of the Jacobian matrix.

    real(wp),dimension(me%n) :: xp  !! perturbed `x` vector
    real(wp),dimension(me%m) :: f0  !! function evaluation `f(x)`
    real(wp),dimension(me%m) :: f1  !! function evaluation `f(x+dx)`

    ! function evaluations:
    xp = x
    call me%compute_function(xp,f0,idx)

    xp(column) = x(column) + dx(column)
    call me%compute_function(xp,f1,idx)

    ! forward difference:
    dfdx = (f1(idx) - f0(idx)) / dx(column)

    end subroutine forward_diff
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
    integer :: nf !! number of functions to compute in a column
    real(wp),dimension(me%n) :: dx !! absolute perturbation (>0) for each variable
    integer,dimension(:),allocatable :: nonzero_elements_in_col !! the indices of the
                                                                !! nonzero Jacobian
                                                                !! elements in a column
    integer,dimension(:),allocatable :: indices  !! index vector
                                                 !! `[1,2,...,num_nonzero_elements]`
                                                 !! for putting `dfdx` into `jac`

    ! if we don't have a sparsity pattern yet then compute it:
    if (.not. me%sparsity_computed) call me%compute_sparsity(x)

    ! initialize:
    allocate(jac(me%num_nonzero_elements))
    jac = 0.0_wp
    indices = [(i,i=1,me%num_nonzero_elements)]

    ! compute dx vector:
    call me%compute_perturbation_vector(x,dx)

    ! compute Jacobian matrix column-by-column:
    do i=1,me%n

        ! determine functions to compute for this column:
        nonzero_elements_in_col = pack(me%irow,mask=me%icol==i)
        nf = size(nonzero_elements_in_col)  ! number of functions to compute
        if (nf/=0) then

            ! possibly inform caller what is happening:
            if (associated(me%info_function)) call me%info_function(i)

            ! size the output array:
            allocate(dfdx(nf))

            ! compute this column of the Jacobian:
            call me%jacobian_method(x,dx,i,nonzero_elements_in_col,dfdx)

            ! put result into the output vector:
            jac(pack(indices,mask=me%icol==i)) = dfdx

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

    ! make sure none are zero:
    where (dx<=0.0_wp) dx = epsilon(1.0_wp)

    end subroutine compute_perturbation_vector
!*******************************************************************************

!*******************************************************************************
!>
!  Print the sparsity pattern in matrix form.

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
    end module numerical_differentiation_module
!*******************************************************************************
