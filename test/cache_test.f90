!*******************************************************************************
!>
!  Tests for [[function_cache]], including the output of `print`.

    program cache_test

    use iso_fortran_env, only: output_unit, error_unit
    use numdiff_cache_module, only: function_cache
    use numdiff_kinds_module, only: wp

    implicit none

    integer :: n_failed = 0  !! number of failed checks

    call test_uninitialized()
    call test_single_slot()
    call test_all_slots_printed()

    if (n_failed==0) then
        write(output_unit,'(A)') 'cache_test: all tests passed'
    else
        write(error_unit,'(A,I0,A)') 'cache_test: ', n_failed, ' check(s) failed'
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

    subroutine print_entries(cache, entries, n_entries, text)
        !! print the cache to a scratch file and return the slot
        !! numbers of the entries that were printed.
        type(function_cache),intent(inout) :: cache
        integer,dimension(:),intent(out) :: entries  !! printed slot numbers
        integer,intent(out) :: n_entries  !! number of entries printed
        character(len=:),allocatable,intent(out) :: text !! all the printed lines (concatenated)
        integer :: iunit, istat
        character(len=1000) :: line
        open(newunit=iunit, status='scratch', form='formatted')
        call cache%print(iunit)
        rewind(iunit)
        n_entries = 0
        text = ''
        do
            read(iunit,'(A)',iostat=istat) line
            if (istat/=0) exit
            text = text//trim(line)//' | '
            if (line(1:5)=='Entry') then
                n_entries = n_entries + 1
                read(line(6:),*) entries(n_entries)
            end if
        end do
        close(iunit)
    end subroutine print_entries

    subroutine test_uninitialized()
        !! printing a cache that was never initialized
        type(function_cache) :: cache
        integer,dimension(10) :: entries
        integer :: n_entries
        character(len=:),allocatable :: text
        call print_entries(cache, entries, n_entries, text)
        call check(n_entries==0, 'uninitialized: no entries printed')
        call check(index(text,'Cache is not initialized')>0, 'uninitialized: message printed')
    end subroutine test_uninitialized

    subroutine test_single_slot()
        !! a cache with one slot: every x hashes to slot 0,
        !! which is the slot the off-by-one loop skipped.
        type(function_cache) :: cache
        integer :: i
        real(wp),dimension(2) :: f
        logical :: xfound
        logical,dimension(2) :: ffound
        integer,dimension(10) :: entries
        integer :: n_entries
        character(len=:),allocatable :: text

        call cache%initialize(isize=1, n=1, m=2)
        call cache%get([1.5_wp], [1,2], i, f, xfound, ffound)
        call check(i==0, 'single slot: index is 0')
        call check(.not. xfound, 'single slot: empty cache')
        call cache%put(i, [1.5_wp], [5.0_wp, 6.0_wp], [1,2])

        call print_entries(cache, entries, n_entries, text)
        call check(n_entries==1, 'single slot: one entry printed')
        if (n_entries==1) call check(entries(1)==0, 'single slot: entry 0 printed')
        call check(index(text,'1.5000000000000000')>0, 'single slot: x printed')
        call check(index(text,'5.0000000000000000')>0 .and. &
                   index(text,'6.0000000000000000')>0, 'single slot: f printed')
    end subroutine test_single_slot

    subroutine test_all_slots_printed()
        !! fill some slots of a larger table, and check that
        !! exactly the occupied slots are printed.
        integer,parameter :: isize = 4
        type(function_cache) :: cache
        integer :: i, k
        real(wp),dimension(1) :: f
        logical :: xfound
        logical,dimension(1) :: ffound
        logical,dimension(0:isize-1) :: occupied
        integer,dimension(isize+1) :: entries
        integer :: n_entries
        character(len=:),allocatable :: text

        call cache%initialize(isize=isize, n=1, m=1)
        occupied = .false.
        do k = 1, 20
            call cache%get([real(k,wp)], [1], i, f, xfound, ffound)
            call check(i>=0 .and. i<isize, 'all slots: index in range')
            if (i<0 .or. i>=isize) return
            call cache%put(i, [real(k,wp)], [real(10*k,wp)], [1])
            occupied(i) = .true.
        end do

        call print_entries(cache, entries, n_entries, text)
        call check(n_entries==count(occupied), 'all slots: every occupied slot printed')
        do k = 1, n_entries
            call check(entries(k)>=0 .and. entries(k)<isize, 'all slots: printed slot in range')
            if (entries(k)>=0 .and. entries(k)<isize) &
                call check(occupied(entries(k)), 'all slots: printed slot is occupied')
        end do
    end subroutine test_all_slots_printed

    end program cache_test
!*******************************************************************************
