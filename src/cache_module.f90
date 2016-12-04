!*******************************************************************************
!> author: Jacob Williams
!
!  For caching function evaluations.

    module cache_module

    use kinds_module
    use iso_fortran_env, only: ip => int64  ! same number of bits are real(wp)

    implicit none

    private

    type :: fx
        !! an [x,f(x)] cached pair.
        private
        real(wp),dimension(:),allocatable :: x    !! vector of input values
        real(wp),dimension(:),allocatable :: f    !! vector of output functions `f(x)`
                                                  !! note: only the elements indicated by `ifs`
                                                  !! will have valid values. The others will
                                                  !! be dummy values.
        integer,dimension(:),allocatable  :: ifs  !! elements of `f` present in the cache
                                                  !! (this is just an array of the indices
                                                  !! present in `f`)
    end type fx

    type,public :: function_cache
        !! a vector function cache.
        private
        integer :: n = 0 !! size of `x`
        integer :: m = 0 !! size of `f`
        type(fx),dimension(:),allocatable :: c  !! the cache of `f(x)`
        integer :: chunk_size = 100 !! for resizing vectors
                                    !! in the [[unique]] function
    contains
        private
        procedure,public :: initialize => initialize_cache
        procedure,public :: get        => get_from_cache
        procedure,public :: put        => put_in_cache
        procedure,public :: destroy    => destroy_cache
        procedure,public :: print      => print_cache
    end type function_cache

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Initialize the cache. Must be called first before use.

    subroutine initialize_cache(me,isize,n,m,chunk_size)

    implicit none

    class(function_cache),intent(inout) :: me
    integer,intent(in) :: isize !! the size of the hash table
    integer,intent(in) :: n     !! number of independant variables (x)
    integer,intent(in) :: m     !! number of functions (f)
    integer,intent(in),optional :: chunk_size  !! chunk size to speed up reallocation
                                               !! of arrays. A good value is a guess for
                                               !! the actual number of elements of `f` that
                                               !! will be saved per value of `x` [default is 100]

    call me%destroy()

    allocate(me%c(0:isize-1))
    me%n = n
    me%m = m

    if (present(chunk_size)) then
        me%chunk_size = chunk_size
    else
        me%chunk_size = 100
    end if

    end subroutine initialize_cache
!*******************************************************************************

!*******************************************************************************
!>
!  Print the contents of the cache. Used for debugging.

    subroutine print_cache(me,iunit)

    implicit none

    class(function_cache),intent(inout) :: me
    integer,intent(in) :: iunit !! file unit for writing
                                !! (assumed to be opened)

    integer :: i  !! counter

    write(iunit,'(A)') ''
    write(iunit,'(A)') '------------------------'
    if (allocated(me%c)) then
        do i = 1, size(me%c)
            if (allocated(me%c(i)%x)) then
                write(iunit,'(A)') ''
                write(iunit,'(A,1X,I10)') 'Entry ',i
                write(iunit,'(A,1X,*(F27.16,1X))') 'x  :', me%c(i)%x
                if (allocated(me%c(i)%f)) then
                    write(iunit,'(A,1X,*(I27,1X))')    'ifs:', me%c(i)%ifs
                    write(iunit,'(A,1X,*(F27.16,1X))') 'f  :', me%c(i)%f(me%c(i)%ifs)
                end if
                write(iunit,'(A)') '------------------------'
            end if
        end do
    else
        write(iunit,'(A)') 'Cache is not initialized'
        write(iunit,'(A)') '------------------------'
    end if
    write(iunit,'(A)') ''

    end subroutine print_cache
!*******************************************************************************

!*******************************************************************************
!>
!  Check if the `x` vector is in the cache, if so return `f`.
!  Note that only some of the elements may be present, so it will return
!  the ones there are there, and indicate which ones were found.

    subroutine get_from_cache(me,x,ifs,i,f,xfound,ffound)

    implicit none

    class(function_cache),intent(inout)      :: me
    real(wp),dimension(:),intent(in)         :: x      !! independant variable vector
    integer,dimension(:),intent(in)          :: ifs    !! elements of `f` needed
    integer,intent(out)                      :: i      !! index in the hash table
    real(wp),dimension(:),intent(out)        :: f      !! `f(x)` from the cache (if it was found)
    logical,intent(out)                      :: xfound !! if `x` was found in the cache
    logical,dimension(size(ifs)),intent(out) :: ffound !! which `ifs` were found in the cache

    integer :: j !! counter

    ! initialize:
    xfound = .false.
    ffound = .false.

    if (allocated(me%c)) then

        ! get index in the hash table:
        i = mod( abs(vector_djb_hash(x)), int(size(me%c),ip) )

        ! check the table:
        if (allocated(me%c(i)%x) .and. allocated(me%c(i)%f)) then
            if (size(me%c(i)%x)==size(x) .and. &
                size(me%c(i)%f)==size(f)) then
                if (all(me%c(i)%x==x)) then
                    xfound = .true.
                    ! return the elements that were found in the cache:
                    f(me%c(i)%ifs) = me%c(i)%f(me%c(i)%ifs)
                    ! what indices are in the cache?
                    do j = 1, size(ifs)
                        ffound(j) = any(ifs(j)==me%c(i)%ifs)
                    end do
                end if
            end if
        end if

    else
        error stop 'Error: the cache has not been initialized.'
    end if

    end subroutine get_from_cache
!*******************************************************************************

!*******************************************************************************
!>
!  Put a value into the cache.

    subroutine put_in_cache(me,i,x,f,ifs)

    use utilities_module,   only: unique

    implicit none

    class(function_cache),intent(inout) :: me
    integer,intent(in)                  :: i    !! index in the hash table
    real(wp),dimension(:),intent(in)    :: x    !! independant variable vector (dimension `n`)
    real(wp),dimension(:),intent(in)    :: f    !! function vector `f(x)` (dimension `m`)
    integer,dimension(:),intent(in)     :: ifs  !! elements of `f` to add (should all be `>0, <=m`)

    real(wp),parameter :: null = huge(1.0_wp) !! an unusual value to initialize arrays

    if (allocated(me%c)) then
        if (i<=size(me%c)) then

            if (allocated(me%c(i)%x)) then
                ! we need to check if there is an x already there,
                ! which may already have some function indices present.
                ! if same x, then add the new indices to them.
                ! if a different x, then replace indices.
                if (all(me%c(i)%x==x)) then
                    ! this x is already present in this location.
                    ! so merge the new f,ifs into what is already there.
                    if (allocated(me%c(i)%f)) then
                        me%c(i)%ifs = unique([me%c(i)%ifs,ifs],&
                                                chunk_size=me%chunk_size)
                    else
                        allocate(me%c(i)%f(me%m))
                        me%c(i)%f   = null ! initialize to an unusual value
                        me%c(i)%ifs = ifs
                    end if
                    me%c(i)%f(ifs) = f(ifs)
                else
                    ! replace existing x and f.
                    me%c(i)%x   = x
                    me%c(i)%ifs = ifs
                    if (allocated(me%c(i)%f)) deallocate(me%c(i)%f)
                    allocate(me%c(i)%f(me%m))
                    me%c(i)%f      = null ! initialize to an unusual value
                    me%c(i)%f(ifs) = f(ifs)
                end if
            else
                ! new entry in the cache:
                me%c(i)%x = x
                allocate(me%c(i)%f(me%m))
                me%c(i)%f      = null ! initialize to an unusual value
                me%c(i)%ifs    = ifs
                me%c(i)%f(ifs) = f(ifs)
            end if

        else
            error stop 'Error: invalid index in hash table.'
        end if
    else
        error stop 'Error: the cache has not been initialized.'
    end if

    end subroutine put_in_cache
!*******************************************************************************

!*******************************************************************************
!>
!  Destroy a cache.

    subroutine destroy_cache(me)

    implicit none

    class(function_cache),intent(out) :: me

    end subroutine destroy_cache
!*******************************************************************************

!*******************************************************************************
!>
!  DJB hash algorithm for a `real(wp)` vector.
!
!### See also
!  * J. Shahbazian, Fortran hashing algorithm, July 6, 2013
!   [Fortran Dev](https://fortrandev.wordpress.com/2013/07/06/fortran-hashing-algorithm/)

    pure function vector_djb_hash(r) result(hash)

    real(wp),dimension(:),intent(in) :: r     !! the vector
    integer(ip)                      :: hash  !! the hash value

    integer :: i !! counter

    hash = 5381_ip

    do i=1,size(r)
        hash = ishft(hash,5_ip) + hash + transfer(r(i),1_ip)
    end do

    end function vector_djb_hash
!*******************************************************************************

!*******************************************************************************
    end module cache_module
!*******************************************************************************
