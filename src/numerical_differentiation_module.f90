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
        procedure :: compute_sparsity_column     !! gets a columns from the sparsity pattern
        procedure :: pack_jacobian_column        !! pack a column of the jacobian into the sparse vector
        procedure :: compute_perturbation_vector !! computes the variable perturbation factor

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
            logical,dimension(:),intent(in) :: funcs_to_compute !! the elements of `f` that
                                                            !! need to be computed (size `m`).
                                                            !! This is computed from the sparsity
                                                            !! pattern.
        end subroutine func
        subroutine spars_f(me,x)
            !! The function to compute the sparsity pattern.
            !! It populates the `irow` and `icol` variables in the class.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout) :: me
            real(wp),dimension(:),intent(in) :: x !! vector of variables (size `n`)
        end subroutine spars_f
        subroutine jac_f(me,x,dx,column,funcs_to_compute,dfdx)
            !! The function to compute a column of jacobian.
            import :: numdiff_type,wp
            implicit none
            class(numdiff_type),intent(inout)  :: me
            real(wp),dimension(:),intent(in)   :: x       !! vector of variables (size `n`)
            real(wp),dimension(:),intent(in)   :: dx      !! absolute perturbation (>0) for each variable
            integer,intent(in)                 :: column  !! column number to compute
            logical,dimension(me%m),intent(in) :: funcs_to_compute   !! the elements in the
                                                                     !! column of the Jacobian
                                                                     !! to compute
            real(wp),dimension(me%m),intent(out)  :: dfdx  !! the column of the full Jacobian
        end subroutine jac_f
    end interface

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the class. This must be called first.

    subroutine initialize_numdiff_type(me,n,m,xlow,xhigh,perturb_mode,dpert)

    implicit none

    class(numdiff_type),intent(inout)   :: me
    integer,intent(in)                  :: n            !! number of `x` variables
    integer,intent(in)                  :: m            !! number of `f` functions
    real(wp),dimension(n),intent(in)    :: xlow         !! lower bounds on `x`
    real(wp),dimension(n),intent(in)    :: xhigh        !! upper bounds on `x`
    integer,intent(in)                  :: perturb_mode !! perturbation mode (1,2,3)
    real(wp),dimension(n),intent(in)    :: dpert        !! perturbation vector for `x`

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
!  Given the column number (variable number), compute the
!  functions that need to be computed, based on the sparsity pattern.

    subroutine compute_sparsity_column(me,ivar,funcs_to_compute)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in) :: ivar !! column number to compute
    logical,dimension(me%m),intent(out) :: funcs_to_compute !! the functions that need
                                                            !! to be computed (non-zero)
                                                            !! elements of the sparsity pattern.

    integer :: i !! counter

    ! initialize:
    funcs_to_compute = .false.

    ! locate the elements for this column:
    do i = 1, me%num_nonzero_elements
        if (me%icol(i)==ivar) then
            funcs_to_compute(me%irow(i)) = .true.
        elseif (me%icol(i)>ivar) then
            ! all the elements in this columns have
            ! been found (data is stored by column)
            exit
        end if
    end do

    end subroutine compute_sparsity_column
!*******************************************************************************

!*******************************************************************************
!>
!  Compute a column of the Jacobian matrix using basic forward differences.

    subroutine forward_diff(me,x,dx,column,funcs_to_compute,dfdx)

    implicit none

    class(numdiff_type),intent(inout)  :: me
    real(wp),dimension(:),intent(in)   :: x       !! vector of variables (size `n`)
    real(wp),dimension(:),intent(in)   :: dx      !! absolute perturbation (>0) for each variable
    integer,intent(in)                 :: column  !! column number to compute
    logical,dimension(me%m),intent(in) :: funcs_to_compute   !! the elements in the
                                                             !! column of the Jacobian
                                                             !! to compute
    real(wp),dimension(me%m),intent(out)  :: dfdx  !! the column of the full Jacobian

    real(wp),dimension(me%n) :: xp  !! perturbed `x` vector
    real(wp),dimension(me%m) :: f0  !! function evaluation `f(x)`
    real(wp),dimension(me%m) :: f1  !! function evaluation `f(x+dx)`

    xp = x
    call me%compute_function(xp,f0,funcs_to_compute)

    xp(column) = x(column) + dx(column)
    call me%compute_function(xp,f1,funcs_to_compute)

    ! note: we do the subtraction for all elements,
    ! although we only need to do them for the non-zero
    ! elements .... maybe change this ...

    dfdx = (f1 - f0) / dx(column)

    end subroutine forward_diff
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the Jacobian using forward differences.

    subroutine compute_jacobian(me,x,j)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x  !! vector of variables (size `n`)
    real(wp),dimension(:),allocatable,intent(out) :: j  !! sparse jacobian vector

    real(wp),dimension(me%m) :: dfdx !! a column of the full Jacobian
    logical,dimension(me%m) :: funcs_to_compute  !! the elements in a given
                                                 !! column of the Jacobian
                                                 !! that are non-zero
    integer :: i  !! column counter
    integer :: nf !! number of functions to compute in a column
    real(wp),dimension(me%n) :: dx !! absolute perturbation (>0) for each variable

    ! if we don't have a sparsity pattern yet then compute it:
    if (.not. me%sparsity_computed) call me%compute_sparsity(x)

    ! initialize:
    allocate(j(me%num_nonzero_elements))
    j = 0.0_wp

    ! compute dx vector:
    call me%compute_perturbation_vector(x,dx)

    ! compute Jacobian matrix column-by-column:
    do i=1,me%n

        ! determine functions to compute for this column:
        call me%compute_sparsity_column(i,funcs_to_compute)
        nf = count(funcs_to_compute)  ! number of functions to compute
        if (nf/=0) then

            ! compute this column of the Jacobian:
            call me%jacobian_method(x,dx,i,funcs_to_compute,dfdx)

            ! put result into the output vector:
            call me%pack_jacobian_column(i,dfdx,nf,j)

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
!  Take a column of the Jacobian, and put it
!  in the sparse `j` jacobian vector

    subroutine pack_jacobian_column(me,column_number,dfdx_column,num_nonzero_funcs,j)

    implicit none

    class(numdiff_type),intent(inout) :: me
    integer,intent(in)                :: column_number     !! the columns number
    real(wp),dimension(:)             :: dfdx_column       !! a column of the full Jacobian
    integer,intent(in)                :: num_nonzero_funcs !! number of non-zero Jacobian
                                                           !! elements in this column
    real(wp),dimension(:),intent(inout) :: j               !! the full sparse jacobian vector

    integer :: k !! counter
    integer :: n_accumulated !! for determing when all elements of
                             !! a Jacobian column have been
                             !! accumulated into `j`

    ! put result into the output vector:
    n_accumulated = 0
    do k=1,me%num_nonzero_elements
        if (me%icol(k)==column_number) then
            n_accumulated = n_accumulated + 1
            j(k) = dfdx_column(me%irow(k))
            if (n_accumulated==num_nonzero_funcs) exit !done with this column
        end if
    end do

    end subroutine pack_jacobian_column
!*******************************************************************************

!*******************************************************************************
    end module numerical_differentiation_module
!*******************************************************************************
