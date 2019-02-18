!*******************************************************************************
!> author: Jacob Williams
!  date: October 31, 2016
!
!  Test1 for the numerical differentiation module.

    program test1

    use iso_fortran_env, only: wp => real64, output_unit
    use numerical_differentiation_module

    implicit none

    integer,parameter :: n = 10
    integer,parameter :: m = 6
    real(wp),dimension(n),parameter :: x     = 1.0_wp
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),dimension(n),parameter :: dpert = 1.0e-5_wp
    integer,parameter :: perturb_mode = 1

    type(numdiff_type)                :: my_prob
    integer                           :: i          !! counter
    real(wp),dimension(:),allocatable :: jac
    character(len=:),allocatable      :: formula
    type(finite_diff_method)          :: fd
    logical                           :: status_ok
    type(meth_array)                  :: meths
    integer                           :: func_evals  !! function evaluation counter

    integer,parameter :: cache_size = 0 !! `0` indicates not to use cache
    !integer,parameter :: cache_size = 1000 ! use cache

    do i=1,44  ! try different finite diff methods

        call get_finite_diff_formula(i,formula)
        if (formula=='') then
            write(output_unit,*) 'id ', i , 'not available'
            cycle
        end if
        func_evals = 0
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_mode=2,&
                                jacobian_method=i,&
                                partition_sparsity_pattern=.false.,&
                                cache_size=cache_size)

        call my_prob%compute_jacobian(x,jac)

        if (i==1) then
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_pattern(output_unit)
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_matrix(output_unit)
            write(output_unit,'(A)') ''

            write(output_unit,'(A)') ''
            write(output_unit,'(A)') '-------------------'
            write(output_unit,'(A,I2)') ' specify method [no partitioning] : ',i
            write(output_unit,'(A)') '-------------------'
            write(output_unit,'(A)') ''
        end if

        write(output_unit,'(A)') formula
        write(output_unit,'(A)') ''
        write(output_unit,'(A,1X,*(F27.16,","))') 'jac =',jac
        write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
        write(output_unit,'(A)') ''

    end do

    do i=1,44  ! try different finite diff methods

        call get_finite_diff_formula(i,formula)
        if (formula=='') then
            write(output_unit,*) 'id ', i , 'not available'
            cycle
        end if

        func_evals = 0
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_mode=2,&
                                jacobian_method=i,&
                                partition_sparsity_pattern=.true.,&
                                cache_size=cache_size)

        call my_prob%compute_jacobian(x,jac)

        if (i==1) then
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_pattern(output_unit)
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_matrix(output_unit)
            write(output_unit,'(A)') ''

            write(output_unit,'(A)') ''
            write(output_unit,'(A)') '-------------------'
            write(output_unit,'(A,I2)') ' specify method [with partitioning] : ',i
            write(output_unit,'(A)') '-------------------'
            write(output_unit,'(A)') ''
        end if

        write(output_unit,'(A)') formula
        write(output_unit,'(A)') ''
        write(output_unit,'(A,1X,*(F27.16,","))') 'jac =',jac
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

    do i=2,9

        write(output_unit,'(A)') ''
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A,I2)') ' specify class [no partitioning] : ', i
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ''

        func_evals = 0
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_mode=2,&
                                class=i,&
                                cache_size=cache_size)

        call my_prob%compute_jacobian(x,jac)
        write(output_unit,'(A,1X,*(F27.16,","))') 'jac =',jac
        write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
        write(output_unit,'(A)') ''

    end do

    do i=2,9

        write(output_unit,'(A)') ''
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A,I2)') ' specify class [with partitioning] : ',i
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ''

        func_evals = 0
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_mode=2,&
                                class=i,partition_sparsity_pattern=.true.,&
                                cache_size=cache_size)

        call my_prob%compute_jacobian(x,jac)
        write(output_unit,'(A,1X,*(F27.16,","))') 'jac =',jac
        write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
        write(output_unit,'(A)') ''

    end do

    write(output_unit,'(A)') ''
    write(output_unit,'(A)') '-------------------'
    write(output_unit,'(A)') ' diff algorithm'
    write(output_unit,'(A)') '-------------------'
    write(output_unit,'(A)') ''

    func_evals = 0
    call my_prob%diff_initialize(n,m,xlow,xhigh,my_func,sparsity_mode=2,&
                                    cache_size=cache_size)

    call my_prob%compute_jacobian(x,jac)
    write(output_unit,'(A,1X,*(F27.16,","))') 'jac =',jac
    write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
    write(output_unit,'(A)') ''

contains

    subroutine my_func(me,x,f,funcs_to_compute)

    !! Problem function

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x
    real(wp),dimension(:),intent(out) :: f
    integer,dimension(:),intent(in)   :: funcs_to_compute

    if (any(funcs_to_compute==1)) f(1) = x(1)*x(2) - x(3)**3
    if (any(funcs_to_compute==2)) f(2) = x(3) - 1.0_wp
    if (any(funcs_to_compute==3)) f(3) = x(4)*x(5)
    if (any(funcs_to_compute==4)) f(4) = 2.0_wp*x(6) + sin(x(7))
    if (any(funcs_to_compute==5)) f(5) = cos(x(8)) + sqrt(x(9))
    if (any(funcs_to_compute==6)) f(6) = 1.0_wp / (1.0_wp + exp(x(10)))

    func_evals = func_evals + 1

    end subroutine my_func

!*******************************************************************************
    end program test1
!*******************************************************************************
