!*******************************************************************************
!> author: Jacob Williams
!  date: October 31, 2016
!
!  Test1 for the numerical differentiation module.

    program test1

    use iso_fortran_env, only: wp => real64, output_unit
    use numerical_differentiation_module

    implicit none

    type(numdiff_type) :: my_prob

    integer,parameter :: n = 3
    integer,parameter :: m = 2
    real(wp),dimension(n),parameter :: x = [1.0_wp,2.0_wp,3.0_wp]
    real(wp),dimension(n),parameter :: xlow = [-10.0_wp,-10.0_wp,-10.0_wp]
    real(wp),dimension(n),parameter :: xhigh = [10.0_wp,10.0_wp,10.0_wp]
    real(wp),dimension(n),parameter :: dpert = [1.0e-6_wp,1.0e-6_wp,1.0e-6_wp]
    integer,parameter :: perturb_mode = 1

    integer :: i !! counter
    real(wp),dimension(:),allocatable :: jac
    character(len=:),allocatable :: formula
    type(finite_diff_method) :: fd
    logical :: status_ok

    do i=1,5  ! try different finite diff methods

        write(output_unit,'(A)') ''
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ' specify method'
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ''

        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_func=compute_sparsity_random,&
                                jacobian_method=i)  ! 1 = forward diffs

        call get_finite_diff_formula(i,formula)
        write(output_unit,'(A)') ''
        write(*,*) formula
        write(output_unit,'(A)') ''

        call my_prob%compute_jacobian(x,jac)

        write(output_unit,'(A)') ''
        call my_prob%print_sparsity_pattern(output_unit)
        write(output_unit,'(A,1X,*(F17.12,","))') 'jac =',jac
        write(output_unit,'(A)') ''
        call my_prob%print_sparsity_matrix(output_unit)
        write(output_unit,'(A)') ''

    end do

    write(output_unit,'(A)') ''
    write(output_unit,'(A)') 'select_finite_diff_method'
    write(output_unit,'(A)') ''
    call my_prob%select_finite_diff_method(0.0_wp,0.0_wp,1.0_wp,0.001_wp,3,fd,status_ok)
    call fd%get_formula(formula)
    write(output_unit,'(A)') formula
    write(output_unit,'(A)') ''

    call my_prob%select_finite_diff_method(0.9999_wp,0.0_wp,1.0_wp,0.001_wp,3,fd,status_ok)
    call fd%get_formula(formula)
    write(output_unit,'(A)') formula
    write(output_unit,'(A)') ''

    do i=2,3

        write(output_unit,'(A)') ''
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ' specify class'
        write(output_unit,'(A)') '-------------------'
        write(output_unit,'(A)') ''

        call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                                problem_func=my_func,&
                                sparsity_func=compute_sparsity_random,&
                                class=i)

        call my_prob%compute_jacobian(x,jac)
        write(output_unit,'(A,1X,*(F17.12,","))') 'jac =',jac
        write(output_unit,'(A)') ''

    end do

contains

    subroutine my_func(me,x,f,indices_to_compute)

    !! Problem function

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x
    real(wp),dimension(:),intent(out) :: f
    integer,dimension(:),intent(in) :: indices_to_compute

    if (any(indices_to_compute==1)) f(1) = x(1)*x(2) - x(3)
    if (any(indices_to_compute==2)) f(2) = x(3) - 1.0_wp

    end subroutine my_func

!*******************************************************************************
    end program test1
!*******************************************************************************
