!*******************************************************************************
!>
!  Tests for [[compute_jacobian_times_vector]].
!
!  The test function has a known Jacobian with a constant element,
!  a zero column, and a zero row:
!```
!  f1 = x1**2 + 3*x2         J = [ 2*x1     3   0   0 ]
!  f2 = x2*x3                    [ 0       x3  x2   0 ]
!  f3 = sin(x1)                  [ cos(x1)  0   0   0 ]
!  f4 = 5                        [ 0        0   0   0 ]
!```

    program jacobian_times_vector_test

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_kinds_module, only: wp

    implicit none

    integer,parameter :: n = 4 !! number of variables
    integer,parameter :: m = 4 !! number of functions
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),parameter :: h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)  !! step for central differences
    real(wp),dimension(n),parameter :: dpert = h
    real(wp),parameter :: tol = 100.0_wp * epsilon(1.0_wp)**(2.0_wp/3.0_wp) !! central difference accuracy
    real(wp),dimension(n),parameter :: x = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
    real(wp),dimension(n),parameter :: v = [0.5_wp, -1.0_wp, 2.0_wp, 7.0_wp]

    integer :: n_failed = 0  !! number of failed checks
    type(numdiff_type) :: prob
    real(wp),dimension(m) :: z
    real(wp),dimension(:,:),allocatable :: jac

    ! dense sparsity pattern:
    call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,problem_func=func,&
                         sparsity_mode=1,jacobian_method=3)
    call compute_and_check('dense', matmul(analytic_jacobian(x),v))

    ! computed sparsity pattern (three-point method), not partitioned and partitioned:
    call prob%destroy()
    call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,problem_func=func,&
                         sparsity_mode=2,jacobian_method=3)
    call compute_and_check('sparsity_mode=2', matmul(analytic_jacobian(x),v))

    call prob%destroy()
    call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,problem_func=func,&
                         sparsity_mode=2,jacobian_method=3,partition_sparsity_pattern=.true.)
    call compute_and_check('sparsity_mode=2, partitioned', matmul(analytic_jacobian(x),v))

    ! computed sparsity pattern (jacobians at several points):
    call prob%destroy()
    call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,problem_func=func,&
                         sparsity_mode=4,jacobian_method=3)
    call compute_and_check('sparsity_mode=4', matmul(analytic_jacobian(x),v))

    ! user-specified pattern, with the constant element (1,2)=3 given
    ! as a linear element (it is not in the nonlinear pattern):
    call prob%destroy()
    call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,problem_func=func,&
                         sparsity_mode=3,jacobian_method=3)
    call prob%set_sparsity_pattern(irow=[1,2,2,3], icol=[1,2,3,1],&
                                   linear_irow=[1], linear_icol=[2], linear_vals=[3.0_wp])
    call compute_and_check('user pattern + linear element', matmul(analytic_jacobian(x),v))

    ! the dense form must give the same product:
    call prob%compute_jacobian_dense(x,jac)
    call check(.not. prob%failed(), 'user pattern + linear element: dense jacobian computed')
    if (allocated(jac)) then
        call prob%compute_jacobian_times_vector(x,v,z)
        call check(all(abs(z - matmul(jac,v)) <= tol*max(1.0_wp,abs(z))), &
                   'user pattern + linear element: consistent with dense jacobian')
    end if

    ! only linear elements (the nonlinear pattern is empty,
    ! so compute_jacobian returns no elements):
    call prob%destroy()
    call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,problem_func=func,&
                         sparsity_mode=3,jacobian_method=3)
    call prob%set_sparsity_pattern(irow=[integer ::], icol=[integer ::],&
                                   linear_irow=[1], linear_icol=[2], linear_vals=[3.0_wp])
    call compute_and_check('linear elements only', [3.0_wp*v(2), 0.0_wp, 0.0_wp, 0.0_wp])

    ! a zero vector gives exactly zero:
    call prob%destroy()
    call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,problem_func=func,&
                         sparsity_mode=1,jacobian_method=3)
    call prob%compute_jacobian_times_vector(x,[0.0_wp,0.0_wp,0.0_wp,0.0_wp],z)
    call check(.not. prob%failed(), 'zero vector: no exception')
    call check(all(z==0.0_wp), 'zero vector: z is exactly zero')

    if (n_failed==0) then
        write(output_unit,'(A)') 'jacobian_times_vector_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'jacobian_times_vector_test: ', n_failed, ' check(s) failed'
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

    subroutine compute_and_check(label, z_expected)
        !! compute `J*v` with the current settings in `prob`
        !! and compare it to the expected value.
        character(len=*),intent(in) :: label
        real(wp),dimension(m),intent(in) :: z_expected
        integer :: i
        character(len=:),allocatable :: error_msg
        call prob%compute_jacobian_times_vector(x,v,z)
        if (prob%failed()) then
            call prob%get_error_status(error_msg=error_msg)
            call check(.false., label//': exception raised: '//error_msg)
            return
        end if
        do i = 1, m
            if (abs(z(i)-z_expected(i)) > tol*max(1.0_wp,abs(z_expected(i)))) then
                call check(.false., label//': wrong value for row '//achar(iachar('0')+i))
                write(error_unit,'(A,2(1X,ES24.16))') '   got, expected:', z(i), z_expected(i)
            end if
        end do
    end subroutine compute_and_check

    pure function analytic_jacobian(x) result(jac)
        !! the exact Jacobian of [[func]]
        real(wp),dimension(n),intent(in) :: x
        real(wp),dimension(m,n) :: jac
        jac = 0.0_wp
        jac(1,1) = 2.0_wp*x(1)
        jac(1,2) = 3.0_wp
        jac(2,2) = x(3)
        jac(2,3) = x(2)
        jac(3,1) = cos(x(1))
    end function analytic_jacobian

    subroutine func(me,x,f,funcs_to_compute)
        !! test function
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        f(1) = x(1)**2 + 3.0_wp*x(2)
        f(2) = x(2)*x(3)
        f(3) = sin(x(1))
        f(4) = 5.0_wp
    end subroutine func

    end program jacobian_times_vector_test
!*******************************************************************************
