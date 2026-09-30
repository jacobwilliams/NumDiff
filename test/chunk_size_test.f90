!*******************************************************************************
!>
!  Test that chunk sizes < 1 are handled (they are treated as 1) in
!  [[expand_vector]], [[unique]], the function cache, and [[numdiff_type]].

    program chunk_size_test

    use iso_fortran_env, only: output_unit, error_unit
    use numerical_differentiation_module
    use numdiff_utilities_module, only: expand_vector, unique
    use numdiff_cache_module, only: function_cache
    use numdiff_kinds_module, only: wp

    implicit none

    integer,dimension(*),parameter :: chunk_sizes = [0, -3, 1, 5, 100]
    integer :: n_failed = 0  !! number of failed checks
    integer :: i
    character(len=30) :: label

    do i = 1, size(chunk_sizes)
        write(label,'(A,I0)') 'chunk_size=', chunk_sizes(i)
        call test_expand_vector(chunk_sizes(i), trim(label))
        call test_unique(chunk_sizes(i), trim(label))
        call test_cache(chunk_sizes(i), trim(label))
        call test_numdiff(chunk_sizes(i), trim(label))
    end do

    if (n_failed==0) then
        write(output_unit,'(A)') 'chunk_size_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'chunk_size_test: ', n_failed, ' check(s) failed'
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

    subroutine test_expand_vector(chunk_size, label)
        !! add 10 elements to integer and real vectors
        integer,intent(in) :: chunk_size
        character(len=*),intent(in) :: label
        integer,dimension(:),allocatable :: ivec
        real(wp),dimension(:),allocatable :: rvec
        integer :: k, ni, nr
        ni = 0
        nr = 0
        do k = 1, 10
            call expand_vector(ivec,ni,chunk_size,val=k)
            call expand_vector(rvec,nr,chunk_size,val=real(k,wp))
        end do
        call expand_vector(ivec,ni,chunk_size,finished=.true.)
        call expand_vector(rvec,nr,chunk_size,finished=.true.)
        call check(size(ivec)==10 .and. all(ivec==[(k,k=1,10)]), label//': expand_vector (integer)')
        call check(size(rvec)==10 .and. all(rvec==[(real(k,wp),k=1,10)]), label//': expand_vector (real)')
    end subroutine test_expand_vector

    subroutine test_unique(chunk_size, label)
        !! unique elements of integer and real vectors
        integer,intent(in) :: chunk_size
        character(len=*),intent(in) :: label
        integer,dimension(:),allocatable :: iu
        real(wp),dimension(:),allocatable :: ru
        iu = unique([3,1,3,2,1], chunk_size)
        ru = unique([3.0_wp,1.0_wp,3.0_wp,2.0_wp], chunk_size)
        call check(size(iu)==3, label//': unique (integer) size')
        if (size(iu)==3) call check(all(iu==[1,2,3]), label//': unique (integer) values')
        call check(size(ru)==3, label//': unique (real) size')
        if (size(ru)==3) call check(all(ru==[1.0_wp,2.0_wp,3.0_wp]), label//': unique (real) values')
    end subroutine test_unique

    subroutine test_cache(chunk_size, label)
        !! add two sets of functions for the same x (merged with unique)
        integer,intent(in) :: chunk_size
        character(len=*),intent(in) :: label
        type(function_cache) :: cache
        integer :: i
        real(wp),dimension(3) :: f
        logical :: xfound
        logical,dimension(3) :: ffound
        call cache%initialize(isize=10, n=1, m=3, chunk_size=chunk_size)
        call cache%get([1.0_wp],[1,2,3],i,f,xfound,ffound)
        call cache%put(i,[1.0_wp],[10.0_wp,20.0_wp,30.0_wp],[1,3])
        call cache%put(i,[1.0_wp],[10.0_wp,20.0_wp,30.0_wp],[2])
        call cache%get([1.0_wp],[1,2,3],i,f,xfound,ffound)
        call check(xfound .and. all(ffound), label//': cache merges function indices')
        if (all(ffound)) call check(all(f==[10.0_wp,20.0_wp,30.0_wp]), label//': cache values')
    end subroutine test_cache

    subroutine test_numdiff(chunk_size, label)
        !! compute a sparsity pattern (which is accumulated with expand_vector)
        integer,intent(in) :: chunk_size
        character(len=*),intent(in) :: label
        type(numdiff_type) :: prob
        integer,dimension(:),allocatable :: irow, icol
        call prob%initialize(3,2,[-10.0_wp,-10.0_wp,-10.0_wp],[10.0_wp,10.0_wp,10.0_wp],&
                             perturb_mode=1,dpert=[1.0e-6_wp,1.0e-6_wp,1.0e-6_wp],&
                             problem_func=func,sparsity_mode=2,jacobian_method=3,&
                             chunk_size=chunk_size)
        call prob%compute_sparsity_pattern([1.0_wp,2.0_wp,3.0_wp],irow,icol)
        call check(.not. prob%failed(), label//': numdiff sparsity computed')
        if (allocated(irow) .and. allocated(icol)) then
            call check(size(irow)==3 .and. size(icol)==3, label//': numdiff sparsity size')
            if (size(irow)==3 .and. size(icol)==3) &
                call check(all(irow==[1,2,2]) .and. all(icol==[1,2,3]), label//': numdiff sparsity pattern')
        else
            call check(.false., label//': numdiff sparsity allocated')
        end if
    end subroutine test_numdiff

    subroutine func(me,x,f,funcs_to_compute)
        !! test function: f1 depends on x1, f2 depends on x2 and x3
        class(numdiff_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer,dimension(:),intent(in)   :: funcs_to_compute
        f(1) = x(1)**2
        f(2) = x(2)*x(3)
    end subroutine func

    end program chunk_size_test
!*******************************************************************************
