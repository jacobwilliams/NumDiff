!*******************************************************************************
!>
!  Tests for [[diff_module]] (Neville's process for the first, second,
!  and third derivatives of a function of one variable).
!
!  The main checks are:
!
!  * the estimated error bound is not exceeded by the actual error,
!  * `ifail` is consistent with the requested tolerance,
!  * invalid inputs and too-small intervals are detected,
!  * a user termination (at any function evaluation) returns `ifail=-1`
!    with no further function evaluations.
!
!  The cases are also chosen to exercise the code paths of [[diff]]:
!  symmetric and one-sided formulas (both directions), all three orders,
!  absolute, relative, and minimal-error tolerances, `x0=0`, very narrow
!  and very wide intervals, and the overflow guards.

    module diff_test_functions

    use numdiff_kinds_module, only: wp
    use diff_module, only: diff_func

    implicit none

    integer :: ncalls  = 0 !! number of function evaluations
    integer :: stop_at = 0 !! if >0, call `terminate` on this function evaluation
    integer :: fn      = 1 !! which test function to use (see [[f]])

    ! test function IDs:
    integer,parameter :: f_exp      = 1 !! `exp(x)`
    integer,parameter :: f_sin      = 2 !! `sin(x)`
    integer,parameter :: f_const    = 3 !! `7`
    integer,parameter :: f_quantize = 4 !! `anint(10*x)/10` (a step function)
    integer,parameter :: f_square   = 5 !! `x**2`
    integer,parameter :: f_linear   = 6 !! `x`
    integer,parameter :: f_big      = 7 !! `big*x` (derivative near overflow)
    integer,parameter :: f_wiggle1  = 8 !! `exp(x) + 1e-6*sin(1e3*x)`
    integer,parameter :: f_wiggle2  = 9 !! `exp(x) + 1e-4*sin(1e2*x)`
    integer,parameter :: f_kink     = 10 !! `min(x, x_kink)`

    real(wp),parameter :: big    = huge(1.0_wp) / 20.0_wp !! slope for `f_big`
    real(wp),parameter :: x_kink = 1.0_wp + 45.0_wp*epsilon(1.0_wp) !! kink location for `f_kink`

    contains

    function f(me,x) result(fx)
        !! the test functions
        class(diff_func),intent(inout) :: me
        real(wp),intent(in) :: x
        real(wp) :: fx
        ncalls = ncalls + 1
        select case (fn)
        case(f_exp);      fx = exp(x)
        case(f_sin);      fx = sin(x)
        case(f_const);    fx = 7.0_wp
        case(f_quantize); fx = anint(x*10.0_wp)/10.0_wp
        case(f_square);   fx = x**2
        case(f_linear);   fx = x
        case(f_big);      fx = big*x
        case(f_wiggle1);  fx = exp(x) + 1.0e-6_wp*sin(1.0e3_wp*x)
        case(f_wiggle2);  fx = exp(x) + 1.0e-4_wp*sin(1.0e2_wp*x)
        case(f_kink);     fx = min(x, x_kink)
        case default;     error stop 'invalid test function'
        end select
        if (stop_at>0 .and. ncalls>=stop_at) call me%terminate()
    end function f

    pure function exact(ifn,iord,x) result(d)
        !! the exact derivative of order `iord` of test function `ifn` at `x`
        integer,intent(in) :: ifn
        integer,intent(in) :: iord
        real(wp),intent(in) :: x
        real(wp) :: d
        select case (ifn)
        case(f_exp)
            d = exp(x)
        case(f_sin)
            select case (iord)
            case(1); d = cos(x)
            case(2); d = -sin(x)
            case default; d = -cos(x)
            end select
        case(f_square)
            select case (iord)
            case(1); d = 2.0_wp*x
            case(2); d = 2.0_wp
            case default; d = 0.0_wp
            end select
        case(f_linear, f_big)
            d = 0.0_wp
            if (iord==1) d = 1.0_wp
            if (iord==1 .and. ifn==f_big) d = big
        case(f_wiggle1)
            d = exp(x) + wiggle(1.0e-6_wp, 1.0e3_wp)
        case(f_wiggle2)
            d = exp(x) + wiggle(1.0e-4_wp, 1.0e2_wp)
        case default
            d = 0.0_wp
        end select
    contains
        pure real(wp) function wiggle(a,b)
            !! derivative of `a*sin(b*x)`
            real(wp),intent(in) :: a, b
            select case (iord)
            case(1); wiggle =  a*b*cos(b*x)
            case(2); wiggle = -a*b**2*sin(b*x)
            case default; wiggle = -a*b**3*cos(b*x)
            end select
        end function wiggle
    end function exact

    end module diff_test_functions
!*******************************************************************************

!*******************************************************************************
    program diff_test

    use iso_fortran_env, only: output_unit, error_unit
    use numdiff_kinds_module, only: wp
    use diff_module, only: diff_func
    use diff_test_functions

    implicit none

    type :: test_case
        !! a set of inputs for `diff`
        character(len=40) :: label = ''
        integer  :: fn   = f_exp
        real(wp) :: x0   = 0.0_wp
        real(wp) :: xmin = 0.0_wp
        real(wp) :: xmax = 0.0_wp
        real(wp) :: eps  = 1.0e-9_wp
        real(wp) :: accr = 0.0_wp
    end type test_case

    ! machine constants used by diff (for the precision-independent narrow intervals):
    integer,parameter  :: inf = -minexponent(1.0_wp) - 2
    integer,parameter  :: sup = maxexponent(1.0_wp) - 1
    real(wp),parameter :: twoinf = 2.0_wp**(-inf)
    real(wp),parameter :: min_width_x0_zero = 128.0_wp*twoinf !! smallest half-width allowed at `x0=0`
    real(wp),parameter :: minh_iord3 = 2.0_wp**(-min(inf,sup)/3) !! smallest step for `iord=3`

    integer :: n_failed = 0  !! number of failed checks
    type(diff_func) :: d
    type(test_case),dimension(:),allocatable :: smooth_cases, terminate_cases
    integer :: i, order

    call test_invalid_inputs()   ! before set_function is called

    call d%set_function(f)

    ! smooth functions: accuracy and error bound, for all orders
    smooth_cases = [ &
        test_case('interior, near xmin side', f_exp, 0.3_wp, -1.0_wp, 2.0_wp), &
        test_case('interior, near xmax side', f_exp, 1.5_wp, -1.0_wp, 2.0_wp), &
        test_case('x0 = xmin (one-sided)',    f_exp, 0.0_wp,  0.0_wp, 1.0_wp), &
        test_case('x0 = xmax (one-sided)',    f_exp, 1.0_wp,  0.0_wp, 1.0_wp), &
        test_case('x0 close to xmin',         f_exp, 1.0e-3_wp, 0.0_wp, 1.0_wp), &
        test_case('x0 = 0',                   f_exp, 0.0_wp, -1.0_wp, 1.0_wp), &
        test_case('x0 = 0, close to xmin',    f_exp, 0.0_wp, -1.0e-3_wp, 1.0_wp), &
        test_case('x0 = 0, close to xmax',    f_exp, 0.0_wp, -1.0_wp, 1.0e-3_wp), &
        test_case('relative eps',             f_exp, 0.3_wp, -1.0_wp, 2.0_wp, eps=-1.0e-9_wp), &
        test_case('minimal error (eps=0)',    f_exp, 0.3_wp, -1.0_wp, 2.0_wp, eps=0.0_wp), &
        test_case('tight eps',                f_exp, 0.3_wp, -1.0_wp, 2.0_wp, eps=1.0e-15_wp), &
        test_case('relative accr',            f_exp, 0.3_wp, -1.0_wp, 2.0_wp, accr=-1.0e-15_wp), &
        test_case('absolute accr',            f_exp, 0.3_wp, -1.0_wp, 2.0_wp, accr=1.0e-15_wp), &
        test_case('wide interval',            f_sin, 0.3_wp, -1.0e7_wp, 1.0e7_wp), &
        test_case('wide interval, one-sided', f_sin, 0.0_wp, 0.0_wp, 1.0e7_wp), &
        test_case('wide interval',            f_square, 0.3_wp, -1.0e7_wp, 1.0e7_wp), &
        test_case('large x0',                 f_square, 1.0e20_wp, 0.0_wp, 2.0e20_wp), &
        test_case('high-frequency term',      f_wiggle1, 0.3_wp, -1.0_wp, 2.0_wp), &
        test_case('high-frequency term',      f_wiggle2, 0.3_wp, -1.0_wp, 2.0_wp), &
        test_case('high-frequency term, eps=0', f_wiggle2, 0.3_wp, -1.0_wp, 2.0_wp, eps=0.0_wp), &
        test_case('narrow interval',          f_exp, 1.0_wp, 1.0_wp-450.0_wp*epsilon(1.0_wp), &
                                                             1.0_wp+450.0_wp*epsilon(1.0_wp)) ]
    ! the same cases for sin(x), which has different values for each order:
    do i = 1, 13
        smooth_cases = [smooth_cases, smooth_cases(i)]
        smooth_cases(size(smooth_cases))%fn = f_sin
    end do
    do order = 1, 3
        do i = 1, size(smooth_cases)
            call test_accuracy(smooth_cases(i), order)
        end do
    end do

    call test_special_cases()

    ! user termination at every function evaluation:
    terminate_cases = [ &
        smooth_cases(1:8), &          ! exp(x)
        smooth_cases(22:29), &        ! sin(x)
        test_case('wide interval', f_sin, 0.3_wp, -1.0e7_wp, 1.0e7_wp), &
        test_case('wide interval, one-sided', f_sin, 0.0_wp, 0.0_wp, 1.0e7_wp), &
        test_case('minimal error (eps=0)', f_exp, 0.3_wp, -1.0_wp, 2.0_wp, eps=0.0_wp), &
        test_case('kink near x0', f_kink, 1.0_wp, 1.0_wp-180.0_wp*epsilon(1.0_wp), &
                                                  1.0_wp+180.0_wp*epsilon(1.0_wp)), &
        test_case('x0 = 0, minimum interval', f_linear, 0.0_wp, -1.5_wp*min_width_x0_zero, &
                                                                1.5_wp*min_width_x0_zero) ]
    do order = 1, 3
        do i = 1, size(terminate_cases)
            call test_terminate(terminate_cases(i), order)
        end do
    end do

    if (n_failed==0) then
        write(output_unit,'(A)') 'diff_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'diff_test: ', n_failed, ' check(s) failed'
        error stop 1
    end if

contains

    subroutine check(ok, msg, tc, iord)
        !! record the result of a check
        logical,intent(in) :: ok
        character(len=*),intent(in) :: msg
        type(test_case),intent(in),optional :: tc
        integer,intent(in),optional :: iord
        if (.not. ok) then
            n_failed = n_failed + 1
            if (present(tc) .and. present(iord)) then
                write(error_unit,'(A,I0,A,I0,A)') 'FAILED: '//trim(tc%label)//' (fn=', tc%fn, &
                                                   ', iord=', iord, '): '//msg
            else
                write(error_unit,'(A)') 'FAILED: '//msg
            end if
        end if
    end subroutine check

    subroutine run(tc, iord, deriv, error, ifail)
        !! call `diff` for a test case
        type(test_case),intent(in) :: tc
        integer,intent(in) :: iord
        real(wp),intent(out) :: deriv, error
        integer,intent(out) :: ifail
        fn = tc%fn
        ncalls = 0
        call d%compute_derivative(iord,tc%x0,tc%xmin,tc%xmax,tc%eps,tc%accr,deriv,error,ifail)
    end subroutine run

    subroutine test_accuracy(tc, iord)
        !! the actual error is within the estimated error, and
        !! `ifail` is consistent with the requested tolerance.
        type(test_case),intent(in) :: tc
        integer,intent(in) :: iord
        real(wp) :: deriv, error, tol
        integer :: ifail
        call run(tc, iord, deriv, error, ifail)
        call check(ifail==0 .or. ifail==1, 'ifail is 0 or 1', tc, iord)
        if (ifail/=0 .and. ifail/=1) return
        call check(abs(deriv - exact(tc%fn,iord,tc%x0)) <= error, &
                   'actual error is within the estimated error', tc, iord)
        if (tc%eps < 0.0_wp) then
            tol = abs(tc%eps*deriv)
        else
            tol = tc%eps
        end if
        if (tc%eps==0.0_wp) then
            call check(ifail==0, 'ifail=0 when minimizing the error', tc, iord)
        else if (ifail==0) then
            call check(error <= tol, 'ifail=0 only if the tolerance is met', tc, iord)
        else
            call check(error > tol, 'ifail=1 only if the tolerance is not met', tc, iord)
        end if
    end subroutine test_accuracy

    subroutine test_invalid_inputs()
        !! invalid inputs return `ifail=2`
        real(wp) :: deriv, error
        integer :: ifail
        type(test_case) :: tc
        tc = test_case('ok', f_exp, 0.3_wp, -1.0_wp, 2.0_wp)
        call run(tc, 1, deriv, error, ifail)
        call check(ifail==2, 'function not set: ifail=2')
        call d%set_function(f)
        call run(tc, 0, deriv, error, ifail)
        call check(ifail==2, 'iord=0: ifail=2')
        call run(tc, 4, deriv, error, ifail)
        call check(ifail==2, 'iord=4: ifail=2')
        call run(test_case('', f_exp, 0.3_wp, 2.0_wp, -1.0_wp), 1, deriv, error, ifail)
        call check(ifail==2, 'xmin>xmax: ifail=2')
        call run(test_case('', f_exp, 0.3_wp, 1.0_wp, 1.0_wp), 1, deriv, error, ifail)
        call check(ifail==2, 'xmin=xmax: ifail=2')
        call run(test_case('', f_exp, -3.0_wp, -1.0_wp, 2.0_wp), 1, deriv, error, ifail)
        call check(ifail==2, 'x0<xmin: ifail=2')
        call run(test_case('', f_exp, 3.0_wp, -1.0_wp, 2.0_wp), 1, deriv, error, ifail)
        call check(ifail==2, 'x0>xmax: ifail=2')
        call check(ncalls==0, 'invalid inputs: no function evaluations')
    end subroutine test_invalid_inputs

    subroutine test_special_cases()
        !! constant functions, intervals that are too small, and the overflow guard
        real(wp) :: deriv, error
        integer :: ifail, iord
        type(test_case) :: tc

        ! constant function:
        do iord = 1, 3
            call run(test_case('', f_const, 0.3_wp, -1.0_wp, 2.0_wp), iord, deriv, error, ifail)
            call check(ifail==0 .and. deriv==0.0_wp .and. error==0.0_wp, 'constant function')
        end do

        ! the interval is too small for x0:
        call run(test_case('', f_exp, 1.0_wp, 1.0_wp-epsilon(1.0_wp), 1.0_wp+epsilon(1.0_wp)), &
                 1, deriv, error, ifail)
        call check(ifail==3, 'interval too small: ifail=3')
        call check(ncalls==0, 'interval too small: no function evaluations')

        ! the function only changes on a scale that is too large for the interval:
        call run(test_case('', f_quantize, 0.33_wp, -1.0_wp, 2.0_wp), 1, deriv, error, ifail)
        call check(ifail==3, 'step function: ifail=3')

        ! at x0=0, the interval is large enough for the first derivative,
        ! but too small for the third:
        tc = test_case('narrow interval at x0=0', f_linear, 0.0_wp, -10.0_wp*minh_iord3, &
                                                                     10.0_wp*minh_iord3)
        call run(tc, 1, deriv, error, ifail)
        call check(ifail==0, 'narrow interval at x0=0, iord=1: ifail=0')
        call check(abs(deriv-1.0_wp) <= error, 'narrow interval at x0=0, iord=1: within error')
        call run(tc, 3, deriv, error, ifail)
        call check(ifail==3, 'narrow interval at x0=0, iord=3: ifail=3')

        ! the smallest interval allowed at x0=0:
        tc = test_case('minimum interval at x0=0', f_linear, 0.0_wp, -1.5_wp*min_width_x0_zero, &
                                                                      1.5_wp*min_width_x0_zero)
        call run(tc, 1, deriv, error, ifail)
        call check(ifail==0 .or. ifail==1, 'minimum interval at x0=0: ifail is 0 or 1')
        call check(abs(deriv-1.0_wp) <= error, 'minimum interval at x0=0: within error')

        ! a kink just inside the interval (the function is not smooth,
        ! so only check that the result is valid):
        do iord = 1, 3
            call run(test_case('', f_kink, 1.0_wp, 1.0_wp-180.0_wp*epsilon(1.0_wp), &
                                                   1.0_wp+180.0_wp*epsilon(1.0_wp)), iord, deriv, error, ifail)
            call check(ifail==1 .or. ifail==3, 'kink: ifail is 1 or 3')
        end do

        ! derivative near overflow: the difference quotient is capped to
        ! avoid overflow, so the result is not accurate, but it is finite
        ! and the requested accuracy is reported as not met.
        ! [note: the estimated error does not account for the cap]
        do iord = 1, 3
            call run(test_case('', f_big, 0.5_wp, 0.0_wp, 1.0_wp), iord, deriv, error, ifail)
            call check(ifail==1, 'derivative near overflow: ifail=1')
            call check(abs(deriv)<=huge(1.0_wp) .and. deriv==deriv, 'derivative near overflow: finite')
        end do
    end subroutine test_special_cases

    subroutine test_terminate(tc, iord)
        !! terminating at each function evaluation returns `ifail=-1`
        !! with no further evaluations, and the next call is not affected.
        type(test_case),intent(in) :: tc
        integer,intent(in) :: iord
        real(wp) :: deriv, error, deriv0, error0
        integer :: ifail, ifail0, k, n
        logical :: ok
        stop_at = 0
        call run(tc, iord, deriv0, error0, ifail0)
        n = ncalls
        ok = .true.
        do k = 1, n
            stop_at = k
            call run(tc, iord, deriv, error, ifail)
            if (ifail/=-1 .or. ncalls/=k) then
                ok = .false.
                exit
            end if
        end do
        stop_at = 0
        call check(ok, 'terminate at every function evaluation', tc, iord)
        ! a subsequent call is not affected by the termination:
        call run(tc, iord, deriv, error, ifail)
        call check(ifail==ifail0 .and. ncalls==n .and. deriv==deriv0 .and. error==error0, &
                   'not affected by a previous termination', tc, iord)
    end subroutine test_terminate

    end program diff_test
!*******************************************************************************
