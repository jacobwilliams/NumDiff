!*******************************************************************************
!> author: Jacob Williams
!  date: October 31, 2016
!
!  Test1 for the numerical differentiation module.

    program test1

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_kinds_module, only: wp

    implicit none

    integer,parameter :: n = 10
    integer,parameter :: m = 6
    real(wp),dimension(n),parameter :: x     = 1.0_wp
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),dimension(n),parameter :: dpert = 1.0e-5_wp
    integer,parameter :: perturb_mode = 1
    integer,parameter :: cache_size = 0 !! `0` indicates not to use cache
    integer,parameter :: sparsity_mode = 4

    type(numdiff_type)                :: my_prob
    integer                           :: i          !! counter
    integer                           :: j          !! counter
    real(wp),dimension(:),allocatable :: jac
    character(len=:),allocatable      :: formula
    type(finite_diff_method)          :: fd
    logical                           :: status_ok
    type(meth_array)                  :: meths
    integer                           :: func_evals  !! function evaluation counter
    integer,dimension(:),allocatable  :: methods     !! array of method IDs
    character(len=:),allocatable      :: error_msg   !! error message string

    methods = [(i, i = 1, 44)]
    methods = [methods, 500,600,700,800]  ! these only have central diffs

    do j=1,size(methods)  ! try different finite diff methods

        i = methods(j)

        call get_finite_diff_formula(i,formula)
        if (formula=='') then
            write(output_unit,*) 'id ', i , 'not available'
            cycle
        end if
        func_evals = 0
        call my_prob%destroy()
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_mode=sparsity_mode,&
                                jacobian_method=i,&
                                partition_sparsity_pattern=.false.,&
                                cache_size=cache_size)
        if (my_prob%failed()) then
            call my_prob%get_error_status(error_msg=error_msg)
            write(error_unit,'(A)') error_msg
            stop
        end if

        call my_prob%compute_jacobian(x,jac)
        if (my_prob%failed()) then
            call my_prob%get_error_status(error_msg=error_msg)
            write(error_unit,'(A)') error_msg
            stop
        end if

        if (i==1) then
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_pattern(output_unit)
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_matrix(output_unit)
            write(output_unit,'(A)') ''

            write(output_unit,'(A)') ''
            write(output_unit,'(A)') '--------------------------------------'
            write(output_unit,'(A,I4)') ' specify method [no partitioning] : ',i
            write(output_unit,'(A)') '--------------------------------------'
            write(output_unit,'(A)') ''
        end if

        write(output_unit,'(A)') formula
        write(output_unit,'(A)') ''
        write(output_unit,'(A,1X,*(F27.16,","))') 'jac =',jac
        write(output_unit,'(A,1X,I5)') 'function evaluations:',func_evals
        write(output_unit,'(A)') ''

    end do

    do j=1,size(methods)  ! try different finite diff methods

        i = methods(j)

        call get_finite_diff_formula(i,formula)
        if (formula=='') then
            write(output_unit,*) 'id ', i , 'not available'
            cycle
        end if

        func_evals = 0
        call my_prob%destroy()
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_mode=sparsity_mode,&
                                jacobian_method=i,&
                                partition_sparsity_pattern=.true.,&
                                cache_size=1000) ! use the cache for these cases
        if (my_prob%failed()) then
            call my_prob%get_error_status(error_msg=error_msg)
            write(error_unit,'(A)') error_msg
            stop
        end if

        call my_prob%compute_jacobian(x,jac)
        if (my_prob%failed()) then
            call my_prob%get_error_status(error_msg=error_msg)
            write(error_unit,'(A)') error_msg
            stop
        end if

        if (i==1) then
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_pattern(output_unit)
            write(output_unit,'(A)') ''
            call my_prob%print_sparsity_matrix(output_unit)
            write(output_unit,'(A)') ''

            write(output_unit,'(A)') ''
            write(output_unit,'(A)') '--------------------------------------'
            write(output_unit,'(A,I4)') ' specify method [with partitioning] : ',i
            write(output_unit,'(A)') '--------------------------------------'
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
        call my_prob%destroy()
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_mode=sparsity_mode,&
                                class=i,&
                                cache_size=cache_size)
        if (my_prob%failed()) then
            call my_prob%get_error_status(error_msg=error_msg)
            write(error_unit,'(A)') error_msg
            stop
        end if

        call my_prob%compute_jacobian(x,jac)
        if (my_prob%failed()) then
            call my_prob%get_error_status(error_msg=error_msg)
            write(error_unit,'(A)') error_msg
            stop
        end if
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
        call my_prob%destroy()
        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_mode=sparsity_mode,&
                                class=i,partition_sparsity_pattern=.true.,&
                                cache_size=cache_size)
        if (my_prob%failed()) then
            call my_prob%get_error_status(error_msg=error_msg)
            write(error_unit,'(A)') error_msg
            stop
        end if

        call my_prob%compute_jacobian(x,jac)
        if (my_prob%failed()) then
            call my_prob%get_error_status(error_msg=error_msg)
            write(error_unit,'(A)') error_msg
            stop
        end if
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
    call my_prob%diff_initialize(n,m,xlow,xhigh,my_func,sparsity_mode=1,&  ! use a dense method for this one
                                    cache_size=cache_size)
    if (my_prob%failed()) then
        call my_prob%get_error_status(error_msg=error_msg)
        write(error_unit,'(A)') error_msg
        stop
    end if

    call my_prob%compute_jacobian(x,jac)
    if (my_prob%failed()) then
        call my_prob%get_error_status(error_msg=error_msg)
        write(error_unit,'(A)') error_msg
        stop
    end if
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
    if (any(funcs_to_compute==5)) f(5) = cos(x(8)) + sqrt(abs(x(9)))
    if (any(funcs_to_compute==6)) f(6) = 1.0_wp / (1.0_wp + exp(x(10)))

    func_evals = func_evals + 1

    end subroutine my_func

!*******************************************************************************
    end program test1
!*******************************************************************************
