!*******************************************************************************
!>
!  Test that the estimated sparsity pattern does not depend on the
!  scale of the functions (the sparsity tolerances are relative).

    program sparsity_scaling_test

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_kinds_module, only: wp

    implicit none

    integer,parameter :: n = 3 !! number of variables
    integer,parameter :: m = 2 !! number of functions
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),dimension(n),parameter :: dpert = 1.0e-5_wp
    real(wp),dimension(n),parameter :: x = [1.0_wp,2.0_wp,3.0_wp]
    integer,dimension(*),parameter :: sparsity_modes = [2,4]
    real(wp),dimension(*),parameter :: scales = [1.0e-20_wp, 1.0_wp, 1.0e20_wp]

    real(wp) :: scale  !! the function scale factor (used in [[func]])
    integer :: n_failed = 0  !! number of failed checks
    integer :: i, j
    type(numdiff_type) :: prob
    integer,dimension(:),allocatable :: irow, icol
    character(len=100) :: label

    do i = 1, size(sparsity_modes)
        do j = 1, size(scales)

            scale = scales(j)
            write(label,'(A,I0,A,ES8.1)') 'sparsity_mode=',sparsity_modes(i),', scale=',scale

            call prob%destroy()
            call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,&
                                 problem_func=func,sparsity_mode=sparsity_modes(i),&
                                 jacobian_method=3)
            call prob%compute_sparsity_pattern(x,irow,icol)

            if (prob%failed()) then
                call fail(trim(label)//': exception raised')
            else if (.not. allocated(irow) .or. .not. allocated(icol)) then
                call fail(trim(label)//': no sparsity pattern')
            else if (size(irow)/=3 .or. size(icol)/=3) then
                call fail(trim(label)//': wrong number of elements')
            else if (any(irow/=[1,2,2]) .or. any(icol/=[1,2,3])) then
                call fail(trim(label)//': wrong sparsity pattern')
            end if

        end do
    end do

    if (n_failed==0) then
        write(output_unit,'(A)') 'sparsity_scaling_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'sparsity_scaling_test: ', n_failed, ' check(s) failed'
        error stop 1
    end if

contains

    subroutine fail(msg)
        !! record a failed check
        character(len=*),intent(in) :: msg
        n_failed = n_failed + 1
        write(error_unit,'(A)') 'FAILED: '//msg
    end subroutine fail

    subroutine func(me,x,f,funcs_to_compute)
        !! test function: f1 depends on x1, f2 depends on x2 and x3.
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        f(1) = scale * x(1)**2
        f(2) = scale * x(2)*x(3)
    end subroutine func

    end program sparsity_scaling_test
!*******************************************************************************
