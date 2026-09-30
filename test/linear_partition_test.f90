!*******************************************************************************
!>
!  Test the partitioned Jacobian when a linear (constant) sparsity pattern
!  is specified. The partition must account for the linear elements: two
!  columns that share a row through a linear element cannot be perturbed
!  together, even though the linear elements themselves are not computed.
!
!  Test function (n=5, m=4):
!```
!  f1 = x1**2 + 3*x2                  nonlinear: (1,1)       linear: (1,2)=3
!  f2 = x2*x3                         nonlinear: (2,2),(2,3)
!  f3 = sin(x4) - 2*x1 + 0.5*x5       nonlinear: (3,4)       linear: (3,1)=-2, (3,5)=0.5
!  f4 = exp(x5) + 4*x3                nonlinear: (4,5)       linear: (4,3)=4
!```
!  Using only the nonlinear pattern, just columns 2 and 3 share a row,
!  so a partition of it would perturb (for example) x1 and x2 together,
!  which corrupts row 1.

    program linear_partition_test

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_kinds_module, only: wp

    implicit none

    integer,parameter :: n = 5 !! number of variables
    integer,parameter :: m = 4 !! number of functions
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),dimension(n),parameter :: x0 = [1.0_wp, 2.0_wp, 3.0_wp, 0.5_wp, 0.25_wp]
    real(wp),dimension(n),parameter :: v = [0.5_wp, -1.0_wp, 2.0_wp, 7.0_wp, -3.0_wp]

    ! the sparsity pattern:
    integer,dimension(*),parameter  :: irow        = [1,2,2,3,4]
    integer,dimension(*),parameter  :: icol        = [1,2,3,4,5]
    integer,dimension(*),parameter  :: linear_irow = [1,3,3,4]
    integer,dimension(*),parameter  :: linear_icol = [2,1,5,3]
    real(wp),dimension(*),parameter :: linear_vals = [3.0_wp,-2.0_wp,0.5_wp,4.0_wp]

    ! step sizes and tolerances for forward and central differences:
    real(wp),parameter :: h_forward = sqrt(epsilon(1.0_wp))
    real(wp),parameter :: h_central = epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    real(wp),parameter :: tol_forward = 100.0_wp * sqrt(epsilon(1.0_wp))
    real(wp),parameter :: tol_central = 100.0_wp * epsilon(1.0_wp)**(2.0_wp/3.0_wp)

    integer :: n_failed = 0  !! number of failed checks

    ! DSM partition, specific methods:
    call test_dsm_partition('forward diffs', h_forward, tol_forward, jacobian_method=1)
    call test_dsm_partition('central diffs', h_central, tol_central, jacobian_method=3)
    ! DSM partition, method class (selected to stay within the bounds):
    call test_dsm_partition('class 3', h_central, tol_central, class=3)

    call test_user_partition()

    if (n_failed==0) then
        write(output_unit,'(A)') 'linear_partition_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'linear_partition_test: ', n_failed, ' check(s) failed'
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

    subroutine initialize(prob, h, partition, jacobian_method, class)
        !! initialize the problem with the user-specified (mode 3) sparsity
        type(numdiff_type),intent(out) :: prob
        real(wp),intent(in) :: h
        logical,intent(in) :: partition
        integer,intent(in),optional :: jacobian_method
        integer,intent(in),optional :: class
        real(wp),dimension(n) :: dpert
        dpert = h
        call prob%initialize(n,m,xlow,xhigh,perturb_mode=1,dpert=dpert,problem_func=func,&
                             sparsity_mode=3,jacobian_method=jacobian_method,class=class,&
                             partition_sparsity_pattern=partition)
    end subroutine initialize

    subroutine check_jacobian(prob, label, tol)
        !! check the dense Jacobian and J*v against the analytic values
        type(numdiff_type),intent(inout) :: prob
        character(len=*),intent(in) :: label
        real(wp),intent(in) :: tol
        real(wp),dimension(:,:),allocatable :: jac
        real(wp),dimension(m,n) :: jac_exact
        real(wp),dimension(m) :: z
        integer :: r, c
        character(len=:),allocatable :: error_msg
        character(len=100) :: element
        jac_exact = analytic_jacobian(x0)
        call prob%compute_jacobian_dense(x0,jac)
        if (prob%failed()) then
            call prob%get_error_status(error_msg=error_msg)
            call check(.false., label//': exception: '//error_msg)
            return
        end if
        do c = 1, n
            do r = 1, m
                if (abs(jac(r,c)-jac_exact(r,c)) > tol*max(1.0_wp,abs(jac_exact(r,c)))) then
                    write(element,'(A,I0,A,I0,A,ES12.4,A,ES12.4)') ' element (',r,',',c,'): got ', &
                                                                     jac(r,c), ', expected ', jac_exact(r,c)
                    call check(.false., label//':'//trim(element))
                end if
            end do
        end do
        call prob%compute_jacobian_times_vector(x0,v,z)
        call check(all(abs(z - matmul(jac_exact,v)) <= tol*max(1.0_wp,abs(z))), label//': J*v')
    end subroutine check_jacobian

    subroutine test_dsm_partition(label, h, tol, jacobian_method, class)
        !! partition computed by DSM
        character(len=*),intent(in) :: label
        real(wp),intent(in) :: h, tol
        integer,intent(in),optional :: jacobian_method, class
        type(numdiff_type) :: prob, prob_unpartitioned
        real(wp),dimension(:),allocatable :: jac, jac_unpartitioned
        integer,dimension(:),allocatable :: ngrp, irow_out, icol_out
        integer :: maxgrp

        call initialize(prob, h, .true., jacobian_method, class)
        call prob%set_sparsity_pattern(irow,icol,linear_irow,linear_icol,linear_vals)
        call check(.not. prob%failed(), label//': set_sparsity_pattern')
        call check_jacobian(prob, label//' [partitioned]', tol)

        ! the partition must group some columns (otherwise this test is not
        ! testing anything), and must be consistent with the full pattern:
        call prob%get_sparsity_pattern(irow_out,icol_out,maxgrp=maxgrp,ngrp=ngrp)
        call check(maxgrp<n, label//': some columns are grouped')
        call check(consistent(ngrp), label//': partition respects the linear elements')

        ! the columns in a group do not affect each other's rows,
        ! so the result must be identical to the unpartitioned one:
        call initialize(prob_unpartitioned, h, .false., jacobian_method, class)
        call prob_unpartitioned%set_sparsity_pattern(irow,icol,linear_irow,linear_icol,linear_vals)
        call prob%compute_jacobian(x0,jac)
        call prob_unpartitioned%compute_jacobian(x0,jac_unpartitioned)
        call check(all(jac==jac_unpartitioned), label//': identical to the unpartitioned jacobian')
    end subroutine test_dsm_partition

    subroutine test_user_partition()
        !! partition specified by the user
        type(numdiff_type) :: prob
        integer :: istat

        ! consistent with the full pattern (columns 2 and 4 share no row,
        ! including through the linear elements):
        call initialize(prob, h_central, .true., jacobian_method=3)
        call prob%set_sparsity_pattern(irow,icol,linear_irow,linear_icol,linear_vals,&
                                       maxgrp=4, ngrp=[1,2,3,2,4])
        call check(.not. prob%failed(), 'consistent user partition: accepted')
        call check_jacobian(prob, 'consistent user partition', tol_central)

        ! only consistent with the nonlinear pattern (columns 1 and 2 share
        ! row 1 through the linear element (1,2)), so it must be rejected:
        call initialize(prob, h_central, .true., jacobian_method=3)
        call prob%set_sparsity_pattern(irow,icol,linear_irow,linear_icol,linear_vals,&
                                       maxgrp=2, ngrp=[1,1,2,1,1])
        call prob%get_error_status(istat=istat)
        call check(prob%failed() .and. istat==33, 'inconsistent user partition: rejected')

        ! without the linear pattern, the same partition is fine:
        call initialize(prob, h_central, .true., jacobian_method=3)
        call prob%set_sparsity_pattern(irow,icol,maxgrp=2,ngrp=[1,1,2,1,1])
        call check(.not. prob%failed(), 'nonlinear-only user partition: accepted')
    end subroutine test_user_partition

    pure logical function consistent(ngrp)
        !! independent check that no two columns in the same group
        !! share a row in the full (nonlinear + linear) pattern
        integer,dimension(:),intent(in) :: ngrp
        integer,dimension(:),allocatable :: r, c
        integer :: i, j
        r = [irow, linear_irow]
        c = [icol, linear_icol]
        consistent = .true.
        do i = 1, size(r)
            do j = 1, size(r)
                if (r(i)==r(j) .and. c(i)/=c(j) .and. ngrp(c(i))==ngrp(c(j))) consistent = .false.
            end do
        end do
    end function consistent

    pure function analytic_jacobian(x) result(jac)
        !! the exact Jacobian of [[func]]
        real(wp),dimension(n),intent(in) :: x
        real(wp),dimension(m,n) :: jac
        jac = 0.0_wp
        jac(1,1) = 2.0_wp*x(1)
        jac(1,2) = 3.0_wp
        jac(2,2) = x(3)
        jac(2,3) = x(2)
        jac(3,1) = -2.0_wp
        jac(3,4) = cos(x(4))
        jac(3,5) = 0.5_wp
        jac(4,3) = 4.0_wp
        jac(4,5) = exp(x(5))
    end function analytic_jacobian

    subroutine func(me,x,f,funcs_to_compute)
        !! test function
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        f(1) = x(1)**2 + 3.0_wp*x(2)
        f(2) = x(2)*x(3)
        f(3) = sin(x(4)) - 2.0_wp*x(1) + 0.5_wp*x(5)
        f(4) = exp(x(5)) + 4.0_wp*x(3)
    end subroutine func

    end program linear_partition_test
!*******************************************************************************
