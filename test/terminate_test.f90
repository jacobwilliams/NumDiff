!*******************************************************************************
!>
!  Test that [[terminate]] stops the computation immediately, when it is
!  called from the user function or from the info function, at any point
!  in the computation (including during the sparsity computation), for
!  both the finite difference methods and the [[diff]] method.

    program terminate_test

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_kinds_module, only: wp

    implicit none

    integer,parameter :: n = 3 !! number of variables
    integer,parameter :: m = 2 !! number of functions
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),dimension(n),parameter :: dpert = 1.0e-6_wp
    real(wp),dimension(n),parameter :: x0 = [1.0_wp, 2.0_wp, 3.0_wp]
    integer,parameter :: n_configs = 6 !! number of configurations to test

    integer :: n_failed = 0     !! number of failed checks
    integer :: nf = 0           !! number of function evaluations
    integer :: ni = 0           !! number of info function calls
    integer :: stop_f = 0       !! terminate in the function on this evaluation (if >0)
    integer :: stop_i = 0       !! terminate in the info function on this call (if >0)
    integer :: nf_at_stop = 0   !! value of `nf` when terminate was called
    integer :: iconfig

    do iconfig = 1, n_configs
        call test_config(iconfig)
    end do

    if (n_failed==0) then
        write(output_unit,'(A)') 'terminate_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'terminate_test: ', n_failed, ' check(s) failed'
        error stop 1
    end if

contains

    subroutine check(ok, msg)
        !! record the result of a check
        logical,intent(in) :: ok
        character(len=*),intent(in) :: msg
        if (.not. ok) then
            n_failed = n_failed + 1
            write(error_unit,'(A)') 'FAILED: '//msg
        end if
    end subroutine check

    subroutine setup(prob, iconfig, label)
        !! initialize the problem for a configuration
        type(numdiff_type),intent(out) :: prob
        integer,intent(in) :: iconfig
        character(len=:),allocatable,intent(out) :: label
        select case (iconfig)
        case(1)
            label = 'diff'
            call prob%diff_initialize(n,m,xlow,xhigh,func,sparsity_mode=1,info=info)
        case(2)
            label = 'diff with cache'
            call prob%diff_initialize(n,m,xlow,xhigh,func,sparsity_mode=1,info=info,cache_size=100)
        case(3)
            label = 'central diffs, dense'
            call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1,&
                                 jacobian_method=3,info=info)
        case(4)
            label = 'forward diffs, sparsity_mode=2, cache'
            call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=2,&
                                 jacobian_method=1,info=info,cache_size=100)
        case(5)
            label = 'class 3, sparsity_mode=4, partitioned'
            call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=4,&
                                 class=3,info=info,partition_sparsity_pattern=.true.)
        case(6)
            label = 'diff, sparsity_mode=4'
            call prob%diff_initialize(n,m,xlow,xhigh,func,sparsity_mode=4,info=info,&
                                      dpert_for_sparsity=dpert,sparsity_perturb_mode=1)
        end select
    end subroutine setup

    subroutine test_config(iconfig)
        !! terminate at every function evaluation and every info call
        integer,intent(in) :: iconfig
        type(numdiff_type) :: prob
        character(len=:),allocatable :: label
        real(wp),dimension(:),allocatable :: jac
        integer :: nf_total, ni_total, k, istat
        character(len=:),allocatable :: error_msg
        logical :: ok_f, ok_i

        ! normal run:
        stop_f = 0; stop_i = 0; nf = 0; ni = 0
        call setup(prob, iconfig, label)
        call prob%compute_jacobian(x0,jac)
        call check(.not. prob%failed(), label//': normal run')
        nf_total = nf
        ni_total = ni

        ! terminate in the function:
        ok_f = .true.
        do k = 1, nf_total
            stop_f = k; stop_i = 0; nf = 0; ni = 0; nf_at_stop = -1
            call setup(prob, iconfig, label)
            call prob%compute_jacobian(x0,jac)
            call prob%get_error_status(istat=istat,error_msg=error_msg)
            if (.not. prob%failed() .or. istat/=-1 .or. nf/=nf_at_stop) then
                ok_f = .false.
                write(error_unit,'(A,I0,A,I0,A,I0)') '  '//label//': terminated at evaluation ', k, &
                                                   ', istat=', istat, ', evaluations=', nf
                exit
            end if
        end do
        call check(ok_f, label//': terminate in the function stops all evaluations')

        ! terminate in the info function:
        ok_i = ni_total>0
        do k = 1, ni_total
            stop_f = 0; stop_i = k; nf = 0; ni = 0; nf_at_stop = -1
            call setup(prob, iconfig, label)
            call prob%compute_jacobian(x0,jac)
            call prob%get_error_status(istat=istat)
            if (.not. prob%failed() .or. istat/=-1 .or. nf/=nf_at_stop) then
                ok_i = .false.
                write(error_unit,'(A,I0,A,I0,A,I0,A,I0,A)') '  '//label//': terminated at info call ', k, &
                                        ', istat=', istat, ', evaluations=', nf, ' (expected ', nf_at_stop, ')'
                exit
            end if
        end do
        call check(ok_i, label//': terminate in the info function stops all evaluations')

        call check(error_msg=='Terminated by the user', label//': error message')
    end subroutine test_config

    subroutine func(me,x,f,funcs_to_compute)
        !! test function
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        nf = nf + 1
        f(1) = x(1)**2 + 3.0_wp*x(2)
        f(2) = x(2)*x(3)
        if (nf==stop_f) then
            nf_at_stop = nf
            call me%terminate()
        end if
    end subroutine func

    subroutine info(me,column,i,x)
        !! info function
        class(numdiff_type),intent(inout) :: me
        integer,dimension(:),intent(in) :: column
        integer,intent(in) :: i
        real(wp),dimension(:),intent(in) :: x
        ni = ni + 1
        if (ni==stop_i) then
            nf_at_stop = nf
            call me%terminate()
        end if
    end subroutine info

    end program terminate_test
!*******************************************************************************
