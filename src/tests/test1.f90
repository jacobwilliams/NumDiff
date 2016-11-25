!*******************************************************************************
!> author: Jacob Williams
!  date: October 31, 2016
!
!  Test1 for the numerical differentiation module.

    program test1

    use iso_fortran_env, only: wp => real64, output_unit
    use numerical_differentiation_module

    implicit none

    integer,parameter :: n = 6
    integer,parameter :: m = 4
    real(wp),dimension(n),parameter :: x     = 1.0_wp
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),dimension(n),parameter :: dpert = 1.0e-5_wp ! 1.0e-6_wp
    integer,parameter :: perturb_mode = 1

    type(numdiff_type)                :: my_prob
    integer                           :: i          !! counter
    real(wp),dimension(:),allocatable :: jac
    character(len=:),allocatable      :: formula
    type(finite_diff_method)          :: fd
    logical                           :: status_ok
    type(meth_array)                  :: meths
    integer                           :: func_evals  !! function evaluation counter

    do i=1,9  ! try different finite diff methods

        func_evals = 0
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_func=compute_sparsity_random,&
                                jacobian_method=i,&
                                partition_sparsity_pattern=.false.)  ! 1 = forward diffs

        call get_finite_diff_formula(i,formula)
        call my_prob%compute_jacobian(x,jac)

        if (i==1) then
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_pattern(output_unit)
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_matrix(output_unit)
            write(output_unit,'(A)') ''

            write(output_unit,'(A)') ''
            write(output_unit,'(A)') '-------------------'
            write(output_unit,'(A)') ' specify method [no partitioning]'
            write(output_unit,'(A)') '-------------------'
            write(output_unit,'(A)') ''
        end if

        write(output_unit,'(A)') formula
        write(output_unit,'(A)') ''
        write(output_unit,'(A,1X,*(F17.12,","))') 'jac =',jac
        write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
        write(output_unit,'(A)') ''

    end do

    do i=1,9  ! try different finite diff methods

        func_evals = 0
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_func=compute_sparsity_random,&
                                jacobian_method=i,&
                                partition_sparsity_pattern=.true.)  ! 1 = forward diffs

        call get_finite_diff_formula(i,formula)
        call my_prob%compute_jacobian(x,jac)

        if (i==1) then
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_pattern(output_unit)
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_matrix(output_unit)
            write(output_unit,'(A)') ''

            write(output_unit,'(A)') ''
            write(output_unit,'(A)') '-------------------'
            write(output_unit,'(A)') ' specify method [with partitioning]'
            write(output_unit,'(A)') '-------------------'
            write(output_unit,'(A)') ''
        end if

        write(output_unit,'(A)') formula
        write(output_unit,'(A)') ''
        write(output_unit,'(A,1X,*(F17.12,","))') 'jac =',jac
        write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
        write(output_unit,'(A)') ''

    end do

    meths = get_all_methods_in_class(3)

    write(output_unit,'(A)') ''
    write(output_unit,'(A)') '-------------------'
    write(output_unit,'(A)') 'select_finite_diff_method'
    write(output_unit,'(A)') '-------------------'
    write(output_unit,'(A)') ''
    call my_prob%select_finite_diff_method(0.0_wp,0.0_wp,1.0_wp,0.001_wp,meths,fd,status_ok)
    call fd%get_formula(formula)
    write(output_unit,'(A)') formula
    write(output_unit,'(A)') ''

    call my_prob%select_finite_diff_method(0.9999_wp,0.0_wp,1.0_wp,0.001_wp,meths,fd,status_ok)
    call fd%get_formula(formula)
    write(output_unit,'(A)') formula
    write(output_unit,'(A)') ''

    do i=2,4

        write(output_unit,'(A)') ''
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ' specify class [no partitioning]'
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ''

        func_evals = 0
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_func=compute_sparsity_random,&
                                class=i)

        call my_prob%compute_jacobian(x,jac)
        write(output_unit,'(A,1X,*(F17.12,","))') 'jac =',jac
        write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
        write(output_unit,'(A)') ''

    end do

    do i=2,4

        write(output_unit,'(A)') ''
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ' specify class [with partitioning]'
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ''

        func_evals = 0
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_func=compute_sparsity_random,&
                                class=i,partition_sparsity_pattern=.true.)

        call my_prob%compute_jacobian(x,jac)
        write(output_unit,'(A,1X,*(F17.12,","))') 'jac =',jac
        write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
        write(output_unit,'(A)') ''

    end do

contains

    subroutine my_func(me,x,f,indices_to_compute)

    !! Problem function

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x
    real(wp),dimension(:),intent(out) :: f
    integer,dimension(:),intent(in)   :: indices_to_compute

    if (any(indices_to_compute==1)) f(1) = x(1)*x(2) - x(3)**3
    if (any(indices_to_compute==2)) f(2) = x(3) - 1.0_wp
    if (any(indices_to_compute==3)) f(3) = x(4)*x(5)
    if (any(indices_to_compute==4)) f(4) = 2.0_wp*x(6)

    func_evals = func_evals + 1

    end subroutine my_func

!*******************************************************************************
    end program test1
!*******************************************************************************
