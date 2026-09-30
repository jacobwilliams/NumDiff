!*******************************************************************************
!>
!  Unit tests for [[divide_interval]], and for the `num_sparsity_points`
!  limits in [[numdiff_type]] that keep it in its valid range.

    program divide_interval_test

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_utilities_module, only: divide_interval
    use numdiff_kinds_module, only: wp

    implicit none

    integer :: n_failed = 0  !! number of failed checks
    integer :: n !! number of points
    real(wp),dimension(:),allocatable :: points
    real(wp),parameter :: tol = 10.0_wp * epsilon(1.0_wp)

    ! documented example:
    points = divide_interval(3)
    call check(size(points)==3, 'n=3: size')
    if (size(points)==3) then
        call check(all(abs(points - [0.25308641972530865_wp, &
                                     0.5061728394506173_wp,  &
                                     0.759259259175926_wp]) <= tol), 'n=3: documented values')
    end if

    ! every allowed number of points gives that many distinct interior points:
    do n = 1, max_num_sparsity_points
        points = divide_interval(n)
        call check(size(points)==n, 'size', n)
        if (size(points)/=n) cycle
        call check(all(points > tol .and. points < 1.0_wp - tol), 'points inside (0,1)', n)
        if (n>1) call check(all(points(2:n) > points(1:n-1)), 'strictly increasing', n)
    end do

    ! num_sparsity_points is validated when initializing:
    call check_num_sparsity_points(0, expect_ok=.false.)
    call check_num_sparsity_points(max_num_sparsity_points+1, expect_ok=.false.)
    call check_num_sparsity_points(1, expect_ok=.true.)
    call check_num_sparsity_points(max_num_sparsity_points, expect_ok=.true.)

    if (n_failed==0) then
        write(output_unit,'(A)') 'divide_interval_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'divide_interval_test: ', n_failed, ' check(s) failed'
        error stop 1
    end if

contains

    subroutine check(ok, name, num_points)
        !! record the result of a check
        logical,intent(in) :: ok
        character(len=*),intent(in) :: name
        integer,intent(in),optional :: num_points  !! number of points being tested
        if (.not. ok) then
            n_failed = n_failed + 1
            if (present(num_points)) then
                write(error_unit,'(A,I0,A)') 'FAILED: n=', num_points, ': '//name
            else
                write(error_unit,'(A)') 'FAILED: '//name
            end if
        end if
    end subroutine check

    subroutine check_num_sparsity_points(num_sparsity_points, expect_ok)
        !! initialize with `sparsity_mode=4` and compute the jacobian of
        !! a simple function, checking that invalid values are rejected
        !! and that valid values give the correct sparsity pattern.
        integer,intent(in) :: num_sparsity_points
        logical,intent(in) :: expect_ok
        integer,parameter :: nv = 3 !! number of variables
        integer,parameter :: m  = 2 !! number of functions
        real(wp),dimension(nv),parameter :: xlow  = -10.0_wp
        real(wp),dimension(nv),parameter :: xhigh = 10.0_wp
        real(wp),dimension(nv),parameter :: dpert = 1.0e-5_wp
        type(numdiff_type) :: prob
        real(wp),dimension(:),allocatable :: jac
        integer,dimension(:),allocatable :: irow, icol
        character(len=40) :: label
        write(label,'(A,I0)') 'num_sparsity_points=', num_sparsity_points
        call prob%initialize(nv,m,xlow,xhigh,perturb_mode=1,dpert=dpert,&
                             problem_func=func,sparsity_mode=4,jacobian_method=3,&
                             num_sparsity_points=num_sparsity_points)
        if (.not. expect_ok) then
            call check(prob%failed(), trim(label)//': rejected')
            return
        end if
        call check(.not. prob%failed(), trim(label)//': accepted')
        call prob%compute_jacobian([1.0_wp,2.0_wp,3.0_wp],jac)
        call check(.not. prob%failed(), trim(label)//': jacobian computed')
        if (prob%failed()) return
        ! f1 depends on x1, f2 depends on x2 and x3:
        call prob%get_sparsity_pattern(irow,icol)
        call check(all(irow==[1,2,2]) .and. all(icol==[1,2,3]), trim(label)//': sparsity pattern')
    end subroutine check_num_sparsity_points

    subroutine func(me,x,f,funcs_to_compute)
        !! test function
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        f(1) = x(1)**2
        f(2) = x(2)*x(3)
    end subroutine func

    end program divide_interval_test
!*******************************************************************************
