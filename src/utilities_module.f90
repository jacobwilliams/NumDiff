!*******************************************************************************
!> author: Jacob Williams
!
!  Utility routines.

    module numdiff_utilities_module

    use numdiff_kinds_module

    integer,parameter :: max_size_for_insertion_sort = 20 !! max size for using insertion sort.

    private

    interface expand_vector
        module procedure :: expand_vector_int, expand_vector_real
    end interface expand_vector
    public :: expand_vector

    interface unique
        module procedure :: unique_int, unique_real
    end interface unique
    public :: unique

    interface sort_ascending
        module procedure :: sort_ascending_int, sort_ascending_real
    end interface sort_ascending
    public :: sort_ascending

    interface swap
        module procedure :: swap_int, swap_real
    end interface swap

    public :: equal_within_tol
    public :: divide_interval

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Add elements to the integer vector in chunks.

    pure subroutine expand_vector_int(vec,n,chunk_size,val,finished)

    implicit none

    integer,dimension(:),allocatable,intent(inout) :: vec
    integer,intent(inout)       :: n           !! counter for last element added to `vec`.
                                               !! must be initialized to `size(vec)`
                                               !! (or 0 if not allocated) before first call
    integer,intent(in)          :: chunk_size  !! allocate `vec` in blocks of this size (>0)
    integer,intent(in),optional :: val         !! the value to add to `vec`
    logical,intent(in),optional :: finished    !! set to true to return `vec`
                                               !! as its correct size (`n`)

    integer,dimension(:),allocatable :: tmp  !! temporary array

    if (present(val)) then
        if (allocated(vec)) then
            if (n==size(vec)) then
                ! have to add another chunk:
                allocate(tmp(size(vec)+chunk_size))
                tmp(1:size(vec)) = vec
                call move_alloc(tmp,vec)
            end if
            n = n + 1
        else
            ! the first element:
            allocate(vec(chunk_size))
            n = 1
        end if
        vec(n) = val
    end if

    if (present(finished)) then
        if (finished) then
            ! set vec to actual size (n):
            if (allocated(tmp)) deallocate(tmp)
            allocate(tmp(n))
            tmp = vec(1:n)
            call move_alloc(tmp,vec)
        end if
    end if

    end subroutine expand_vector_int
!*******************************************************************************

!*******************************************************************************
!>
!  Add elements to the real vector in chunks.

    pure subroutine expand_vector_real(vec,n,chunk_size,val,finished)

    implicit none

    real(wp),dimension(:),allocatable,intent(inout) :: vec
    integer,intent(inout)       :: n           !! counter for last element added to `vec`.
                                               !! must be initialized to `size(vec)`
                                               !! (or 0 if not allocated) before first call
    integer,intent(in)          :: chunk_size  !! allocate `vec` in blocks of this size (>0)
    real(wp),intent(in),optional :: val        !! the value to add to `vec`
    logical,intent(in),optional :: finished    !! set to true to return `vec`
                                               !! as its correct size (`n`)

    real(wp),dimension(:),allocatable :: tmp  !! temporary array

    if (present(val)) then
        if (allocated(vec)) then
            if (n==size(vec)) then
                ! have to add another chunk:
                allocate(tmp(size(vec)+chunk_size))
                tmp(1:size(vec)) = vec
                call move_alloc(tmp,vec)
            end if
            n = n + 1
        else
            ! the first element:
            allocate(vec(chunk_size))
            n = 1
        end if
        vec(n) = val
    end if

    if (present(finished)) then
        if (finished) then
            ! set vec to actual size (n):
            if (allocated(tmp)) deallocate(tmp)
            allocate(tmp(n))
            tmp = vec(1:n)
            call move_alloc(tmp,vec)
        end if
    end if

    end subroutine expand_vector_real
!*******************************************************************************

!*******************************************************************************
!>
!  Returns only the unique elements of the vector (sorted in ascending order).

    function unique_int(vec,chunk_size) result(ivec_unique)

    implicit none

    integer,dimension(:),intent(in)    :: vec         !! a vector of integers
    integer,intent(in)                 :: chunk_size  !! chunk size for adding to arrays
    integer,dimension(:),allocatable   :: ivec_unique !! unique elements of `ivec`

    integer,dimension(size(vec)) :: ivec !! temp copy of vec
    integer :: i !! counter
    integer :: n !! number of unique elements

    ! first we sort it:
    ivec = vec ! make a copy
    call sort_ascending(ivec)

    ! add the first element:
    n = 1
    ivec_unique = [ivec(1)]

    ! walk through array and get the unique ones:
    if (size(ivec)>1) then
        do i = 2, size(ivec)
            if (ivec(i)/=ivec(i-1)) then
                call expand_vector(ivec_unique,n,chunk_size,val=ivec(i))
            end if
        end do
        call expand_vector(ivec_unique,n,chunk_size,finished=.true.)
    end if

    end function unique_int
!*******************************************************************************

!*******************************************************************************
!>
!  Returns only the unique elements of the vector (sorted in ascending order).

    function unique_real(vec,chunk_size) result(ivec_unique)

    implicit none

    real(wp),dimension(:),intent(in)   :: vec         !! a vector of integers
    integer,intent(in)                 :: chunk_size  !! chunk size for adding to arrays
    real(wp),dimension(:),allocatable  :: ivec_unique !! unique elements of `ivec`

    real(wp),dimension(size(vec)) :: ivec !! temp copy of vec
    integer :: i !! counter
    integer :: n !! number of unique elements

    ! first we sort it:
    ivec = vec ! make a copy
    call sort_ascending(ivec)

    ! add the first element:
    n = 1
    ivec_unique = [ivec(1)]

    ! walk through array and get the unique ones:
    if (size(ivec)>1) then
        do i = 2, size(ivec)
            if (ivec(i)/=ivec(i-1)) then
                call expand_vector(ivec_unique,n,chunk_size,val=ivec(i))
            end if
        end do
        call expand_vector(ivec_unique,n,chunk_size,finished=.true.)
    end if

    end function unique_real
!*******************************************************************************

!*******************************************************************************
!>
!  Sorts an integer array `ivec` in increasing order.
!  Uses a basic recursive quicksort
!  (with insertion sort for partitions with \(\le\) 20 elements).

    subroutine sort_ascending_int(ivec)

    implicit none

    integer,dimension(:),intent(inout) :: ivec

    call quicksort(1,size(ivec))

    contains

        recursive subroutine quicksort(ilow,ihigh)

        !! Sort the array

        implicit none

        integer,intent(in) :: ilow
        integer,intent(in) :: ihigh

        integer :: ipivot !! pivot element
        integer :: i      !! counter
        integer :: j      !! counter

        if ( ihigh-ilow<=max_size_for_insertion_sort .and. ihigh>ilow ) then

            ! do insertion sort:
            do i = ilow + 1,ihigh
                do j = i,ilow + 1,-1
                    if ( ivec(j) < ivec(j-1) ) then
                        call swap(ivec(j),ivec(j-1))
                    else
                        exit
                    end if
                end do
            end do

        elseif ( ihigh-ilow>max_size_for_insertion_sort ) then

            ! do the normal quicksort:
            call partition(ilow,ihigh,ipivot)
            call quicksort(ilow,ipivot - 1)
            call quicksort(ipivot + 1,ihigh)

        end if

        end subroutine quicksort

        subroutine partition(ilow,ihigh,ipivot)

        !! Partition the array, based on the
        !! lexical ivecing comparison.

        implicit none

        integer,intent(in)  :: ilow
        integer,intent(in)  :: ihigh
        integer,intent(out) :: ipivot

        integer :: i,ip

        call swap(ivec(ilow),ivec((ilow+ihigh)/2))
        ip = ilow
        do i = ilow + 1, ihigh
            if ( ivec(i) < ivec(ilow) ) then
                ip = ip + 1
                call swap(ivec(ip),ivec(i))
            end if
        end do
        call swap(ivec(ilow),ivec(ip))
        ipivot = ip

        end subroutine partition

    end subroutine sort_ascending_int
!*******************************************************************************

!*******************************************************************************
!>
!  Sorts a real array `ivec` in increasing order.
!  Uses a basic recursive quicksort
!  (with insertion sort for partitions with \(\le\) 20 elements).

    subroutine sort_ascending_real(ivec)

    implicit none

    real(wp),dimension(:),intent(inout) :: ivec

    call quicksort(1,size(ivec))

    contains

        recursive subroutine quicksort(ilow,ihigh)

        !! Sort the array

        implicit none

        integer,intent(in) :: ilow
        integer,intent(in) :: ihigh

        integer :: ipivot !! pivot element
        integer :: i      !! counter
        integer :: j      !! counter

        if ( ihigh-ilow<=max_size_for_insertion_sort .and. ihigh>ilow ) then

            ! do insertion sort:
            do i = ilow + 1,ihigh
                do j = i,ilow + 1,-1
                    if ( ivec(j) < ivec(j-1) ) then
                        call swap(ivec(j),ivec(j-1))
                    else
                        exit
                    end if
                end do
            end do

        elseif ( ihigh-ilow>max_size_for_insertion_sort ) then

            ! do the normal quicksort:
            call partition(ilow,ihigh,ipivot)
            call quicksort(ilow,ipivot - 1)
            call quicksort(ipivot + 1,ihigh)

        end if

        end subroutine quicksort

        subroutine partition(ilow,ihigh,ipivot)

        !! Partition the array, based on the
        !! lexical ivecing comparison.

        implicit none

        integer,intent(in)  :: ilow
        integer,intent(in)  :: ihigh
        integer,intent(out) :: ipivot

        integer :: i,ip

        call swap(ivec(ilow),ivec((ilow+ihigh)/2))
        ip = ilow
        do i = ilow + 1, ihigh
            if ( ivec(i) < ivec(ilow) ) then
                ip = ip + 1
                call swap(ivec(ip),ivec(i))
            end if
        end do
        call swap(ivec(ilow),ivec(ip))
        ipivot = ip

        end subroutine partition

    end subroutine sort_ascending_real
!*******************************************************************************

!*******************************************************************************
!>
!  Swap two integer values.

    pure elemental subroutine swap_int(i1,i2)

    implicit none

    integer,intent(inout) :: i1
    integer,intent(inout) :: i2

    integer :: tmp

    tmp = i1
    i1  = i2
    i2  = tmp

    end subroutine swap_int
!*******************************************************************************

!*******************************************************************************
!>
!  Swap two integer values.

    pure elemental subroutine swap_real(i1,i2)

    implicit none

    real(wp),intent(inout) :: i1
    real(wp),intent(inout) :: i2

    real(wp) :: tmp

    tmp = i1
    i1  = i2
    i2  = tmp

    end subroutine swap_real
!*******************************************************************************

!*******************************************************************************
!>
!  Returns true if the values in the array are the same
!  (to within the specified absolute tolerance).

    pure function equal_within_tol(vals,tol) result (equal)

    implicit none

    real(wp),dimension(:),intent(in) :: vals  !! a set of values
    real(wp),intent(in)              :: tol   !! a positive tolerance value
    logical                          :: equal !! true if they are equal
                                              !! within the tolerance

    equal = all ( abs(vals - vals(1)) <= abs(tol) )

    end function equal_within_tol
!*******************************************************************************

!*******************************************************************************
!>
!  Returns a set of slightly randomized equally-spaced
!  points that divide an interval.
!
!### Example:
!
!  for `num_points` = 3:
!```
!     o---|---|---|---o
!         1   2   3
!```
!  returns: `[0.25308641972530865, 0.5061728394506173, 0.759259259175926]`.

    function divide_interval(num_points) result(points)

    implicit none

    integer,intent(in)                :: num_points  !! the number of points in the interval
    real(wp),dimension(:),allocatable :: points      !! the resultant vector

    real(wp),parameter :: noise = 1.012345678901234567_wp !! a noise value. Not a round number
                                                          !! so as to avoid freak zeros in the
                                                          !! jacobian
    real(wp),parameter :: min_val = 10.0_wp * epsilon(1.0_wp) !! the minimize distance from the lower bound
    real(wp),parameter :: max_val = 1.0_wp - min_val          !! the minimize distance from the upper bound

    integer :: i !! counter
    real(wp) :: delta !! step size
    real(wp),dimension(:),allocatable :: tmp !! a temp array to hold the values

    delta = 1.0_wp / (num_points + 1)

    allocate(tmp(num_points))
    do i = 1, num_points
        tmp(i) = min(max(min_val,delta*i*noise),max_val)
    end do
    ! this is to protect for the min/max case if there
    ! are enough points so that some are duplicated near
    ! the bounds:
    points = unique(tmp,chunk_size=10)

    end function divide_interval
!*******************************************************************************

!*******************************************************************************
    end module numdiff_utilities_module
!*******************************************************************************
