!*******************************************************************************
!>
!  Tests for the [[numdiff_type]] options that are not covered elsewhere:
!
!  * per-variable methods (`jacobian_methods`) and classes (`classes`)
!  * perturbation modes 2 and 3, and [[set_dpert]]
!  * selecting methods within the variable bounds (including the case
!    where no method in the class fits), with and without partitioning
!  * separate sparsity bounds, `dpert_for_sparsity` and `sparsity_perturb_mode`
!  * the [[diff]] method when `x` is outside the bounds
!  * functions that do not depend on `x` (empty sparsity pattern)
!
!  The test function records the perturbation actually applied to each
!  variable, and whether any evaluation was outside the bounds.

    module options_test_functions

    use numerical_differentiation_module, only: numdiff_type
    use numdiff_kinds_module, only: wp

    implicit none

    integer,parameter :: n = 3 !! number of variables
    integer,parameter :: m = 3 !! number of functions

    real(wp),dimension(n) :: x_nominal = 0.0_wp !! the nominal point (for recording the steps)
    real(wp),dimension(n) :: step = 0.0_wp      !! largest perturbation of each variable from `x_nominal`
    real(wp),dimension(n) :: x_min = 0.0_wp     !! smallest value of each variable evaluated
    real(wp),dimension(n) :: x_max = 0.0_wp     !! largest value of each variable evaluated
    logical :: constant = .false.               !! if true, the function does not depend on `x`

    contains

    subroutine reset(x)
        !! reset the recorded values
        real(wp),dimension(n),intent(in) :: x
        x_nominal = x
        step  = 0.0_wp
        x_min = huge(1.0_wp)
        x_max = -huge(1.0_wp)
    end subroutine reset

    subroutine func(me,x,f,funcs_to_compute)
        !! the test function:
        !!```
        !!  f1 = x1**2 + x2
        !!  f2 = sin(x2)*x3
        !!  f3 = exp(x3) + x1*x3
        !!```
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        step  = max(step, abs(x - x_nominal))
        x_min = min(x_min, x)
        x_max = max(x_max, x)
        ! only the requested functions are computed (the others are zero):
        f = 0.0_wp
        if (constant) then
            if (any(funcs_to_compute==1)) f(1) = 1.0_wp
            if (any(funcs_to_compute==2)) f(2) = 2.0_wp
            if (any(funcs_to_compute==3)) f(3) = 3.0_wp
        else
            if (any(funcs_to_compute==1)) f(1) = x(1)**2 + x(2)
            if (any(funcs_to_compute==2)) f(2) = sin(x(2))*x(3)
            if (any(funcs_to_compute==3)) f(3) = exp(x(3)) + x(1)*x(3)
        end if
    end subroutine func

    pure function analytic_jacobian(x) result(jac)
        !! the exact Jacobian of [[func]]
        real(wp),dimension(n),intent(in) :: x
        real(wp),dimension(m,n) :: jac
        jac = 0.0_wp
        jac(1,1) = 2.0_wp*x(1)
        jac(1,2) = 1.0_wp
        jac(2,2) = cos(x(2))*x(3)
        jac(2,3) = sin(x(2))
        jac(3,1) = x(3)
        jac(3,3) = exp(x(3)) + x(1)
    end function analytic_jacobian

    end module options_test_functions
!*******************************************************************************

!*******************************************************************************
    program options_test

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_kinds_module, only: wp
    use options_test_functions

    implicit none

    real(wp),parameter :: eps = epsilon(1.0_wp)
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),dimension(n),parameter :: x0 = [1.0_wp, 2.0_wp, 0.5_wp]

    ! step sizes for each column: forward, central, 5-point central
    real(wp),dimension(n),parameter :: dpert_mixed = [sqrt(eps), eps**(1.0_wp/3.0_wp), eps**(1.0_wp/5.0_wp)]
    ! tolerances for each column:
    real(wp),dimension(n),parameter :: tol_mixed = 1000.0_wp * [sqrt(eps), eps**(2.0_wp/3.0_wp), eps**(4.0_wp/5.0_wp)]
    ! for central differences:
    real(wp),parameter :: h_central = eps**(1.0_wp/3.0_wp)
    real(wp),parameter :: tol_central = 1000.0_wp * eps**(2.0_wp/3.0_wp)

    integer :: n_failed = 0  !! number of failed checks

    call test_jacobian_methods_and_classes()
    call test_perturb_modes()
    call test_set_dpert()
    call test_bounds()
    call test_sparsity_options()
    call test_diff_outside_bounds()
    call test_constant_function()
    call test_cache()

    if (n_failed==0) then
        write(output_unit,'(A)') 'options_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'options_test: ', n_failed, ' check(s) failed'
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

    subroutine check_jacobian(prob, x, tol, label, jac)
        !! compute the dense jacobian and compare it to the analytic one.
        !! `tol` is the tolerance for each column.
        type(numdiff_type),intent(inout) :: prob
        real(wp),dimension(n),intent(in) :: x
        real(wp),dimension(n),intent(in) :: tol
        character(len=*),intent(in) :: label
        real(wp),dimension(:,:),allocatable,intent(out),optional :: jac
        real(wp),dimension(:,:),allocatable :: j
        real(wp),dimension(m,n) :: jexact
        integer :: c
        character(len=:),allocatable :: error_msg
        call reset(x)
        call prob%compute_jacobian_dense(x,j)
        if (prob%failed()) then
            call prob%get_error_status(error_msg=error_msg)
            call check(.false., label//': exception: '//error_msg)
            return
        end if
        jexact = analytic_jacobian(x)
        do c = 1, n
            call check(all(abs(j(:,c)-jexact(:,c)) <= tol(c)*max(1.0_wp,abs(jexact(:,c)))), &
                       label//': column '//achar(iachar('0')+c))
        end do
        if (present(jac)) jac = j
    end subroutine check_jacobian

    subroutine test_jacobian_methods_and_classes()
        !! per-variable methods and classes
        type(numdiff_type) :: prob
        real(wp),dimension(:,:),allocatable :: jac_methods, jac_classes
        integer :: istat

        ! forward, central, and 5-point central:
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert_mixed,func,sparsity_mode=1,&
                             jacobian_methods=[1,3,10])
        call check_jacobian(prob, x0, tol_mixed, 'jacobian_methods', jac_methods)

        ! the first method in each of these classes (when within the bounds)
        ! is the same as above, so the result must be identical:
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert_mixed,func,sparsity_mode=1,&
                             classes=[2,3,5])
        call check_jacobian(prob, x0, tol_mixed, 'classes', jac_classes)
        if (allocated(jac_methods) .and. allocated(jac_classes)) &
            call check(all(jac_methods==jac_classes), 'classes: same as the equivalent jacobian_methods')

        ! these can't be used with a partitioned sparsity pattern:
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert_mixed,func,sparsity_mode=1,&
                             jacobian_methods=[1,3,10],partition_sparsity_pattern=.true.)
        call prob%get_error_status(istat=istat)
        call check(istat==10, 'jacobian_methods with partitioning: exception 10')
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert_mixed,func,sparsity_mode=1,&
                             classes=[2,3,5],partition_sparsity_pattern=.true.)
        call prob%get_error_status(istat=istat)
        call check(istat==11, 'classes with partitioning: exception 11')
    end subroutine test_jacobian_methods_and_classes

    subroutine test_perturb_modes()
        !! the perturbation applied for each perturb_mode
        type(numdiff_type) :: prob
        real(wp),dimension(n),parameter :: x = [0.0_wp, 2.0_wp, -0.5_wp]  ! includes a zero
        real(wp),dimension(n) :: dpert, dx
        integer :: mode

        dpert = h_central
        do mode = 1, 3
            call prob%destroy()
            call prob%initialize(n,m,xlow,xhigh,mode,dpert,func,sparsity_mode=1,jacobian_method=3)
            call check_jacobian(prob, x, [tol_central,tol_central,tol_central], &
                                'perturb_mode='//achar(iachar('0')+mode))
            select case (mode)
            case(1); dx = dpert
            case(2); dx = abs(dpert*x)
                     where (dx < eps) dx = dpert  ! x=0 falls back to dpert
            case(3); dx = dpert*(1.0_wp + abs(x))
            end select
            call check(all(abs(step-dx) <= 10.0_wp*eps*max(1.0_wp,abs(x))), &
                       'perturb_mode='//achar(iachar('0')+mode)//': perturbation')
        end do
    end subroutine test_perturb_modes

    subroutine test_set_dpert()
        !! change the perturbation after initializing
        type(numdiff_type) :: prob
        real(wp),dimension(n),parameter :: dpert2 = [2.0e-5_wp, 3.0e-5_wp, 4.0e-5_wp]
        integer :: istat
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,[1.0e-5_wp,1.0e-5_wp,1.0e-5_wp],func,&
                             sparsity_mode=1,jacobian_method=3)
        call prob%set_dpert(-dpert2)  ! the sign is ignored
        call check_jacobian(prob, x0, [1.0e-6_wp,1.0e-6_wp,1.0e-6_wp], 'set_dpert')
        call check(all(abs(step-dpert2) <= 10.0_wp*eps*max(1.0_wp,abs(x0))), 'set_dpert: perturbation')
        call prob%set_dpert([1.0e-5_wp])
        call prob%get_error_status(istat=istat)
        call check(istat==29, 'set_dpert: wrong size: exception 29')
    end subroutine test_set_dpert

    subroutine test_bounds()
        !! selecting a method within the variable bounds (class mode)
        type(numdiff_type) :: prob
        real(wp),dimension(n) :: dpert, x, lo, hi, tol
        logical :: partition
        integer :: i

        dpert = h_central
        tol = tol_central
        do i = 1, 2
            partition = i==2

            ! x1 at the upper bound, x3 at the lower bound: a one-sided
            ! method must be selected so the bounds are not violated:
            x  = [xhigh(1), x0(2), xlow(3)]
            call prob%destroy()
            call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1,class=3,&
                                 partition_sparsity_pattern=partition)
            call check_jacobian(prob, x, tol, 'x on the bounds, partition='//merge('T','F',partition))
            call check(all(x_min>=xlow) .and. all(x_max<=xhigh), &
                       'x on the bounds: no evaluations outside the bounds, partition='//merge('T','F',partition))

            ! the interval for x2 is narrower than the perturbation, so
            ! no method in the class fits. The first method is used
            ! (with a warning), so the bounds are exceeded, but the
            ! result is still correct (the function is defined there):
            lo = xlow
            hi = xhigh
            lo(2) = x0(2) - 0.5_wp*dpert(2)
            hi(2) = x0(2) + 0.5_wp*dpert(2)
            call prob%destroy()
            call prob%initialize(n,m,lo,hi,1,dpert,func,sparsity_mode=1,class=3,&
                                 partition_sparsity_pattern=partition)
            write(error_unit,'(A)') '[the following bounds warning is expected]'
            call check_jacobian(prob, x0, tol, 'no method fits, partition='//merge('T','F',partition))
            call check(x_min(2)<lo(2) .and. x_max(2)>hi(2), &
                       'no method fits: the first method is used, partition='//merge('T','F',partition))
        end do
    end subroutine test_bounds

    subroutine test_sparsity_options()
        !! separate bounds and perturbations for computing the sparsity
        type(numdiff_type) :: prob
        real(wp),dimension(n),parameter :: xlow_s  = [0.5_wp, 1.5_wp, 0.25_wp]
        real(wp),dimension(n),parameter :: xhigh_s = [1.5_wp, 2.5_wp, 0.75_wp]
        real(wp),dimension(n),parameter :: dpert_s = [1.0e-4_wp, 2.0e-4_wp, 3.0e-4_wp]
        integer,dimension(:),allocatable :: irow, icol
        integer :: istat

        ! sparsity_mode=2 samples within the sparsity bounds:
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,[h_central,h_central,h_central],func,&
                             sparsity_mode=2,jacobian_method=3,&
                             xlow_for_sparsity=xlow_s,xhigh_for_sparsity=xhigh_s)
        call reset(x0)
        call prob%compute_sparsity_pattern(x0,irow,icol)
        call check(.not. prob%failed(), 'sparsity bounds (mode 2): no exception')
        call check(all(x_min>=xlow_s) .and. all(x_max<=xhigh_s), 'sparsity bounds (mode 2): samples within bounds')
        call check_pattern(irow, icol, 'sparsity bounds (mode 2)')

        ! sparsity_mode=4 with a separate perturbation for the sparsity:
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,[h_central,h_central,h_central],func,&
                             sparsity_mode=4,jacobian_method=3,&
                             xlow_for_sparsity=xlow_s,xhigh_for_sparsity=xhigh_s,&
                             dpert_for_sparsity=dpert_s,sparsity_perturb_mode=1)
        call reset(x0)
        call prob%compute_sparsity_pattern(x0,irow,icol)
        call check(.not. prob%failed(), 'sparsity dpert (mode 4): no exception')
        call check(all(x_min>=xlow_s) .and. all(x_max<=xhigh_s), 'sparsity dpert (mode 4): samples within bounds')
        call check_pattern(irow, icol, 'sparsity dpert (mode 4)')
        ! each variable is perturbed by dpert_s from each sample point, and the
        ! sample points are strictly inside the bounds, so:
        call check(all(x_max - x_min > dpert_s), 'sparsity dpert (mode 4): perturbation used')

        ! invalid inputs:
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,[h_central,h_central,h_central],func,&
                             sparsity_mode=2,jacobian_method=3,&
                             xlow_for_sparsity=xhigh_s,xhigh_for_sparsity=xlow_s)
        call prob%get_error_status(istat=istat)
        call check(istat==7, 'sparsity bounds reversed: exception 7')
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,[h_central,h_central,h_central],func,&
                             sparsity_mode=4,jacobian_method=3,sparsity_perturb_mode=4)
        call prob%get_error_status(istat=istat)
        call check(istat==14, 'invalid sparsity_perturb_mode: exception 14')
    end subroutine test_sparsity_options

    subroutine check_pattern(irow, icol, label)
        !! check the sparsity pattern of [[func]] (in column order)
        integer,dimension(:),allocatable,intent(in) :: irow, icol
        character(len=*),intent(in) :: label
        if (.not. allocated(irow) .or. .not. allocated(icol)) then
            call check(.false., label//': pattern computed')
        else if (size(irow)/=6 .or. size(icol)/=6) then
            call check(.false., label//': pattern size')
        else
            call check(all(irow==[1,3,1,2,2,3]) .and. all(icol==[1,1,2,2,3,3]), label//': pattern')
        end if
    end subroutine check_pattern

    subroutine test_diff_outside_bounds()
        !! the diff method requires x to be within the bounds
        type(numdiff_type) :: prob
        real(wp),dimension(:),allocatable :: jac
        integer :: istat
        call prob%destroy()
        call prob%diff_initialize(n,m,xlow,xhigh,func,sparsity_mode=1)
        call check_jacobian(prob, x0, [1.0e-8_wp,1.0e-8_wp,1.0e-8_wp], 'diff')
        call prob%compute_jacobian([xhigh(1)+1.0_wp, x0(2), x0(3)],jac)
        call prob%get_error_status(istat=istat)
        call check(istat==24, 'diff with x outside the bounds: exception 24')
    end subroutine test_diff_outside_bounds

    subroutine test_constant_function()
        !! a function that does not depend on x gives a zero jacobian,
        !! with and without partitioning, for each sparsity mode.
        type(numdiff_type) :: prob
        real(wp),dimension(:,:),allocatable :: jac
        integer :: mode, i
        logical :: partition
        character(len=60) :: label
        constant = .true.
        do i = 1, 2
            partition = i==2
            do mode = 1, 4
                write(label,'(A,I0,A,L1)') 'constant function, sparsity_mode=', mode, ', partition=', partition
                call prob%destroy()
                call prob%initialize(n,m,xlow,xhigh,1,[h_central,h_central,h_central],func,&
                                     sparsity_mode=mode,jacobian_method=3,&
                                     partition_sparsity_pattern=partition)
                if (mode==3) call prob%set_sparsity_pattern([integer ::],[integer ::])
                call prob%compute_jacobian_dense(x0,jac)
                call check(.not. prob%failed(), trim(label)//': no exception')
                if (allocated(jac)) call check(all(jac==0.0_wp), trim(label)//': zero jacobian')
            end do
        end do
        constant = .false.
    end subroutine test_constant_function

    subroutine test_cache()
        !! the function cache must not change the results.
        !!
        !! At the second point, x1 and x2 are so large that the perturbation
        !! does not change them, so columns 1 and 2 are both evaluated at
        !! the nominal point, with different rows (column 1: rows 1,3;
        !! column 2: rows 1,2). This gives a partial cache hit.
        type(numdiff_type) :: prob, prob_cached
        real(wp),dimension(:),allocatable :: jac, jac_cached
        real(wp),dimension(n,2) :: x
        integer :: i
        x(:,1) = x0
        x(:,2) = [1.0e20_wp, 1.0e20_wp, 1.0_wp]
        do i = 1, 2
            call prob%destroy()
            call prob%initialize(n,m,-2.0e20_wp*[1,1,1],2.0e20_wp*[1,1,1],1,[h_central,h_central,h_central],&
                                 func,sparsity_mode=3,jacobian_method=3)
            call prob%set_sparsity_pattern([1,3,1,2,2,3],[1,1,2,2,3,3])
            call prob_cached%destroy()
            call prob_cached%initialize(n,m,-2.0e20_wp*[1,1,1],2.0e20_wp*[1,1,1],1,[h_central,h_central,h_central],&
                                        func,sparsity_mode=3,jacobian_method=3,cache_size=100)
            call prob_cached%set_sparsity_pattern([1,3,1,2,2,3],[1,1,2,2,3,3])
            call prob%compute_jacobian(x(:,i),jac)
            call prob_cached%compute_jacobian(x(:,i),jac_cached)
            call check(.not. prob%failed() .and. .not. prob_cached%failed(), 'cache: no exception')
            if (allocated(jac) .and. allocated(jac_cached)) &
                call check(all(jac==jac_cached), 'cache: same result as without the cache, point '//achar(iachar('0')+i))
        end do
    end subroutine test_cache

    end program options_test
!*******************************************************************************
