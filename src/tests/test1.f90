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

    real(wp),dimension(:),allocatable :: jac

    call my_prob%initialize(n,m,xlow,xhigh,perturb_mode,dpert,&
                            problem_func=my_func,&
                            sparsity_func=compute_sparsity_dense,&
                            jacobian_func=forward_diff)

    call my_prob%compute_jacobian(x,jac)

    write(*,*) ''
    call my_prob%print_sparsity_pattern(output_unit)
    write(output_unit,'(A,1X,*(F10.3,","))') 'jac =',jac

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
