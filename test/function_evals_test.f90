!*******************************************************************************
!>
!  Test that function evaluations are not wasted:
!
!  * sparsity detection (`sparsity_mode=2` and `4`) skips the remaining
!    points for a column once every row is known to be nonzero.
!  * a function value passed in as `fx` is used instead of computing \( f(x) \).

    program function_evals_test

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_kinds_module, only: wp

    implicit none

    integer,parameter :: n = 6 !! number of variables (and functions)
    real(wp),dimension(n),parameter :: xlow  = -10.0_wp
    real(wp),dimension(n),parameter :: xhigh = 10.0_wp
    real(wp),dimension(n),parameter :: dpert = 1.0e-5_wp
    real(wp),dimension(n),parameter :: x = [1.0_wp,2.0_wp,3.0_wp,4.0_wp,5.0_wp,6.0_wp]

    integer :: n_failed = 0  !! number of failed checks
    integer :: func_evals = 0  !! function evaluation counter
    logical :: dense = .false. !! which test function to use (see [[func]])
    integer :: i, k
    type(numdiff_type) :: prob
    integer,dimension(:),allocatable :: irow, icol
    real(wp),dimension(:),allocatable :: jac1, jac2
    real(wp),dimension(n) :: f
    character(len=100) :: label

    ! sparsity detection: dense problem, where every column can stop early:
    dense = .true.
    call check_sparsity(2, n+1)       ! was 2n+1
    call check_sparsity(4, n+1)       ! was 3n+3 (with 3 points)

    ! sparsity detection: tridiagonal problem, where no column can stop early
    ! (a zero element has to be confirmed at every point):
    dense = .false.
    call check_sparsity(2, 2*n+1)
    call check_sparsity(4, 3*n+3)

    ! passing in the function value at x (forward differences, partitioned tridiagonal):
    irow = [1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6]
    icol = [1,2,1,2,3,2,3,4,3,4,5,4,5,6,5,6]
    do k = 1, 2
        write(label,'(A,L1)') 'fx, partitioned=', k==2
        call prob%destroy()
        call prob%initialize(n,n,xlow,xhigh,perturb_mode=1,dpert=dpert,&
                             problem_func=func,sparsity_mode=3,&
                             jacobian_method=1,partition_sparsity_pattern=(k==2))
        call prob%set_sparsity_pattern(irow,icol)

        func_evals = 0
        call prob%compute_jacobian(x,jac1)
        if (func_evals /= merge(4,n+1,k==2)) call fail(trim(label)//': wrong number of evaluations without fx')

        call func(prob,x,f,[(i,i=1,n)])
        func_evals = 0
        call prob%compute_jacobian(x,jac2,fx=f)
        if (func_evals /= merge(3,n,k==2)) call fail(trim(label)//': wrong number of evaluations with fx')

        if (prob%failed()) then
            call fail(trim(label)//': exception raised')
        else if (any(jac1/=jac2)) then
            call fail(trim(label)//': Jacobian differs when fx is given')
        end if
    end do

    if (n_failed==0) then
        write(output_unit,'(A)') 'function_evals_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'function_evals_test: ', n_failed, ' check(s) failed'
        error stop 1
    end if

contains

    subroutine check_sparsity(sparsity_mode, expected_evals)
        !! compute the sparsity pattern and check it and the number of evaluations
        integer,intent(in) :: sparsity_mode
        integer,intent(in) :: expected_evals
        integer :: r, c, nnz
        write(label,'(A,I0,A,L1)') 'sparsity_mode=',sparsity_mode,', dense=',dense
        call prob%destroy()
        call prob%initialize(n,n,xlow,xhigh,perturb_mode=1,dpert=dpert,&
                             problem_func=func,sparsity_mode=sparsity_mode,&
                             jacobian_method=1)
        func_evals = 0
        call prob%compute_sparsity_pattern(x,irow,icol)
        if (prob%failed()) then
            call fail(trim(label)//': exception raised')
            return
        end if
        if (func_evals /= expected_evals) then
            write(error_unit,'(A,I0,A,I0)') 'function evaluations: ', func_evals, ' expected: ', expected_evals
            call fail(trim(label)//': wrong number of evaluations')
        end if
        ! check the pattern:
        nnz = 0
        do c = 1, n
            do r = 1, n
                if (dense .or. abs(r-c)<=1) then
                    nnz = nnz + 1
                    if (.not. any(irow==r .and. icol==c)) call fail(trim(label)//': missing element')
                end if
            end do
        end do
        if (size(irow)/=nnz) call fail(trim(label)//': wrong number of elements')
    end subroutine check_sparsity

    subroutine fail(msg)
        !! record a failed check
        character(len=*),intent(in) :: msg
        n_failed = n_failed + 1
        write(error_unit,'(A)') 'FAILED: '//msg
    end subroutine fail

    subroutine func(me,x,f,funcs_to_compute)
        !! test function: dense (every `f` depends on every `x`) or tridiagonal.
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        integer :: j, r
        func_evals = func_evals + 1
        do j = 1, size(funcs_to_compute)
            r = funcs_to_compute(j)
            if (dense) then
                f(r) = exp(0.01_wp*r*sum(x))
            else
                f(r) = x(r)**2
                if (r>1) f(r) = f(r) + sin(x(r-1))
                if (r<n) f(r) = f(r) + x(r)*x(r+1)
            end if
        end do
    end subroutine func

    end program function_evals_test
!*******************************************************************************
