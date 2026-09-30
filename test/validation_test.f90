!*******************************************************************************
!>
!  Tests for the input validation in [[numdiff_type]] (each invalid
!  input raises the expected exception), and for the output routines
!  (formulas, printing, and returning the sparsity pattern).

    program validation_test

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

    integer :: n_failed = 0  !! number of failed checks

    call test_initialize_errors()
    call test_other_errors()
    call test_error_status()
    call test_formulas()
    call test_sparsity_output()

    if (n_failed==0) then
        write(output_unit,'(A)') 'validation_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'validation_test: ', n_failed, ' check(s) failed'
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

    subroutine check_error(prob, istat_expected, label, msg_contains)
        !! check that an exception with the expected code was raised,
        !! and that the Jacobian is then not computed.
        type(numdiff_type),intent(inout) :: prob
        integer,intent(in) :: istat_expected
        character(len=*),intent(in) :: label
        character(len=*),intent(in),optional :: msg_contains !! text expected in the message
        integer :: istat
        character(len=:),allocatable :: error_msg
        real(wp),dimension(:),allocatable :: jac
        character(len=60) :: codes
        call prob%get_error_status(istat=istat, error_msg=error_msg)
        write(codes,'(A,I0,A,I0)') ' (got ', istat, ', expected ', istat_expected
        call check(prob%failed() .and. istat==istat_expected, label//trim(codes)//')')
        if (present(msg_contains)) &
            call check(index(error_msg,msg_contains)>0, label//': message: '//error_msg)
        ! nothing is computed after an exception:
        call prob%compute_jacobian(x0,jac)
        call check(.not. allocated(jac), label//': no jacobian after an exception')
        call prob%get_error_status(istat=istat)
        call check(istat==istat_expected, label//': exception is preserved')
    end subroutine check_error

    subroutine test_initialize_errors()
        !! invalid inputs to initialize and diff_initialize
        type(numdiff_type) :: prob
        real(wp),dimension(n) :: lo

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1,jacobian_method=45)
        call check_error(prob, 8, 'invalid jacobian_method')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1,jacobian_methods=[1,45,3])
        call check_error(prob, 9, 'invalid jacobian_methods')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1)
        call check_error(prob, 12, 'no method specified')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1,jacobian_method=1,class=3)
        call check_error(prob, 12, 'two method options specified')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,4,dpert,func,sparsity_mode=1,jacobian_method=1)
        call check_error(prob, 13, 'invalid perturb_mode')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=5,jacobian_method=1)
        call check_error(prob, 5, 'invalid sparsity_mode')

        lo = xlow
        lo(2) = xhigh(2)
        call prob%destroy()
        call prob%initialize(n,m,lo,xhigh,1,dpert,func,sparsity_mode=1,jacobian_method=1)
        call check_error(prob, 4, 'xlow >= xhigh', msg_contains='Error for optimization variable 2')

        call prob%destroy()
        call prob%diff_initialize(n,m,lo,xhigh,func,sparsity_mode=1)
        call check_error(prob, 4, 'diff: xlow >= xhigh')

        call prob%destroy()
        call prob%diff_initialize(n,m,xlow,xhigh,func,sparsity_mode=4)
        call check_error(prob, 2, 'diff: sparsity_mode=4 without dpert_for_sparsity')

        call prob%destroy()
        call prob%diff_initialize(n,m,xlow,xhigh,func,sparsity_mode=4,&
                                  dpert_for_sparsity=dpert,sparsity_perturb_mode=0)
        call check_error(prob, 1, 'diff: invalid sparsity_perturb_mode')
    end subroutine test_initialize_errors

    subroutine test_other_errors()
        !! invalid inputs to the other routines
        type(numdiff_type) :: prob

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1,jacobian_method=1)
        call prob%set_numdiff_bounds([0.0_wp],[1.0_wp])
        call check_error(prob, 3, 'set_numdiff_bounds: wrong size')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=3,jacobian_method=1)
        call prob%set_sparsity_pattern([1,2],[1])
        call check_error(prob, 15, 'set_sparsity_pattern: size mismatch')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=3,jacobian_method=1)
        call prob%set_sparsity_pattern([1,2],[1,n+1])
        call check_error(prob, 15, 'set_sparsity_pattern: column out of range')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=3,jacobian_method=1)
        call prob%set_sparsity_pattern([1],[1],linear_irow=[0],linear_icol=[1],linear_vals=[1.0_wp])
        call check_error(prob, 17, 'set_sparsity_pattern: invalid linear pattern')

        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=3,jacobian_method=1,&
                             partition_sparsity_pattern=.true.)
        call prob%set_sparsity_pattern([1],[1],maxgrp=1,ngrp=[1,2,1])
        call check_error(prob, 28, 'set_sparsity_pattern: ngrp > maxgrp')
    end subroutine test_other_errors

    subroutine test_error_status()
        !! no error, then a user termination
        type(numdiff_type) :: prob
        integer :: istat
        character(len=:),allocatable :: error_msg
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1,jacobian_method=1)
        call prob%get_error_status(istat=istat, error_msg=error_msg)
        call check(.not. prob%failed() .and. istat==0 .and. error_msg=='', 'no error: istat=0, empty message')
        call prob%terminate()
        call check_error(prob, -1, 'terminate', msg_contains='Terminated by the user')
    end subroutine test_error_status

    subroutine test_formulas()
        !! finite difference formulas, and printing a method
        character(len=:),allocatable :: formula, name
        type(finite_diff_method) :: fd, fd_empty
        type(numdiff_type) :: prob
        logical :: status_ok
        character(len=:),allocatable :: text

        call get_finite_diff_formula(1, formula, name)
        call check(formula=='dfdx = (f(x+h)-f(x)) / h' .and. name=='2-point forward 1', 'formula: id 1')
        call get_finite_diff_formula(3, formula, name)
        call check(formula=='dfdx = (f(x+h)-f(x-h)) / (2h)' .and. name=='3-point central', 'formula: id 3')
        call get_finite_diff_formula(10, formula)
        call check(formula=='dfdx = (f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h)) / (12h)', 'formula: id 10')
        call get_finite_diff_formula(45, formula, name)
        call check(formula=='' .and. name=='', 'formula: unknown id')
        call fd_empty%get_formula(formula)
        call check(formula=='', 'formula: uninitialized method')

        ! get a method, and print it:
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=1,jacobian_method=1)
        call prob%select_finite_diff_method(0.0_wp,-1.0_wp,1.0_wp,1.0e-3_wp,&
                                            get_all_methods_in_class(3),fd,status_ok)
        call check(status_ok, 'select_finite_diff_method: status_ok')
        text = print_to_string(fd)
        call check(index(text,'name          : 3-point central')>0, 'print method: name')
        call check(index(text,'dx_factors    :     1,   -1,')>0, 'print method: dx_factors')
        call check(index(text,'df_den_factor :     2')>0, 'print method: df_den_factor')
    end subroutine test_formulas

    function print_to_string(fd) result(text)
        !! print a method to a scratch file, and return the output
        type(finite_diff_method),intent(in) :: fd
        character(len=:),allocatable :: text
        integer :: iunit
        open(newunit=iunit, status='scratch', form='formatted')
        call fd%print(iunit)
        text = read_back(iunit)
    end function print_to_string

    function read_back(iunit) result(text)
        !! rewind a scratch file, read all the lines, and close it
        integer,intent(in) :: iunit
        character(len=:),allocatable :: text
        integer :: istat
        character(len=1000) :: line
        rewind(iunit)
        text = ''
        do
            read(iunit,'(A)',iostat=istat) line
            if (istat/=0) exit
            text = text//trim(line)//' | '
        end do
        close(iunit)
    end function read_back

    subroutine test_sparsity_output()
        !! returning and printing a sparsity pattern with linear elements and a partition
        type(numdiff_type) :: prob
        integer,dimension(:),allocatable :: irow, icol, linear_irow, linear_icol, ngrp
        real(wp),dimension(:),allocatable :: linear_vals
        integer :: maxgrp, iunit
        character(len=:),allocatable :: text

        ! f1 = x1**2 (nonlinear), f2 = x2*x3 (nonlinear), plus a linear element (1,3)=5
        call prob%destroy()
        call prob%initialize(n,m,xlow,xhigh,1,dpert,func,sparsity_mode=3,jacobian_method=3,&
                             partition_sparsity_pattern=.true.)
        call prob%set_sparsity_pattern([1,2,2],[1,2,3],linear_irow=[1],linear_icol=[3],linear_vals=[5.0_wp])
        call check(.not. prob%failed(), 'sparsity output: set_sparsity_pattern')

        call prob%get_sparsity_pattern(irow,icol,linear_irow,linear_icol,linear_vals,maxgrp,ngrp)
        call check(all(irow==[1,2,2]) .and. all(icol==[1,2,3]), 'get_sparsity_pattern: pattern')
        call check(all(linear_irow==[1]) .and. all(linear_icol==[3]) .and. all(linear_vals==[5.0_wp]), &
                   'get_sparsity_pattern: linear pattern')
        call check(maxgrp>=1 .and. size(ngrp)==n, 'get_sparsity_pattern: partition')

        ! vector form:
        open(newunit=iunit, status='scratch', form='formatted')
        call prob%print_sparsity_pattern(iunit)
        text = read_back(iunit)
        call check(index(text,'irow:   1,  2,  2,')>0 .and. index(text,'icol:   1,  2,  3,')>0, &
                   'print_sparsity_pattern: pattern')
        call check(index(text,'---Sparsity partition---')>0, 'print_sparsity_pattern: partition')
        call check(index(text,'---Linear sparsity pattern---')>0 .and. index(text,'vals:')>0, &
                   'print_sparsity_pattern: linear pattern')

        ! matrix form:
        open(newunit=iunit, status='scratch', form='formatted')
        call prob%print_sparsity_matrix(iunit)
        text = read_back(iunit)
        ! nonlinear rows: X00, 0XX. linear rows: 00X, 000
        call check(index(text,'---Sparsity pattern--- | X00 | 0XX |')>0, 'print_sparsity_matrix: pattern')
        call check(index(text,'---Linear sparsity pattern--- | 00X | 000 |')>0, 'print_sparsity_matrix: linear pattern')
    end subroutine test_sparsity_output

    subroutine func(me,x,f,funcs_to_compute)
        !! test function
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        f(1) = x(1)**2 + 5.0_wp*x(3)
        f(2) = x(2)*x(3)
    end subroutine func

    end program validation_test
!*******************************************************************************
