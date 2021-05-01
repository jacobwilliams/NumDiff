!*******************************************************************************
!>
!  Jacobian partitioning using the DSM algorithm.
!
!### Reference
!  * Argonne National Laboratory. MINPACK Project. July 1983.
!    Thomas F. Coleman, Burton S. Garbow, Jorge J. More
!  * Thomas F. Coleman, Burton S. Garbow, Jorge J. More, "Algorithm 618:
!    FORTRAN subroutines for estimating sparse Jacobian Matrices",
!    ACM Transactions on Mathematical Software (TOMS),
!    Volume 10 Issue 3, Sept. 1984, Pages 346-347
!
!### History
!  * Jacob Williams, Nov. 2016, extensive refactoring into modern Fortran.

    module dsm_module

    use numdiff_kinds_module

    implicit none

    private

    public :: dsm
    public :: fdjs

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  The purpose of `dsm` is to determine an optimal or near-
!  optimal consistent partition of the columns of a sparse
!  `m` by `n` matrix `a`.
!
!  the sparsity pattern of the matrix `a` is specified by
!  the arrays `indrow` and `indcol`. on input the indices
!  for the non-zero elements of `a` are
!
!  `indrow(k),indcol(k), k = 1,2,...,npairs`.
!
!  the (`indrow`,`indcol`) pairs may be specified in any order.
!  duplicate input pairs are permitted, but the subroutine
!  eliminates them.
!
!  the subroutine partitions the columns of `a` into groups
!  such that columns in the same group do not have a
!  non-zero in the same row position. a partition of the
!  columns of `a` with this property is consistent with the
!  direct determination of `a`.

    subroutine dsm(m,n,Npairs,Indrow,Indcol,Ngrp,Maxgrp,Mingrp,Info,Ipntr,Jpntr)

    implicit none

      integer,intent(in)  :: m       !! number of rows of `a` (>0)
      integer,intent(in)  :: n       !! number of columns of `a` (>0)
      integer,intent(in)  :: npairs  !! number of (`indrow`,`indcol`) pairs used
                                     !! to describe the sparsity pattern of `a` (>0)
      integer,intent(out) :: maxgrp  !! the number of groups in the partition
                                     !! of the columns of `a`.
      integer,intent(out) :: mingrp  !! a lower bound for the number of groups
                                     !! in any consistent partition of the
                                     !! columns of `a`.
      integer,intent(out) :: info    !! for normal termination `info = 1`.
                                     !! if `m`, `n`, or `npairs` is not positive,
                                     !! then `info = 0`. if the k-th element of
                                     !! `indrow` is not an integer between
                                     !! 1 and m or the k-th element of `indcol`
                                     !! is not an integer between 1 and n,
                                     !! then `info = -k`.
      integer,dimension(npairs),intent(inout) :: indrow  !! an integer array of length `npairs`. on input `indrow`
                                                         !! must contain the row indices of the non-zero elements of `a`.
                                                         !! on output `indrow` is permuted so that the corresponding
                                                         !! column indices are in non-decreasing order. the column
                                                         !! indices can be recovered from the array `jpntr`.
      integer,dimension(npairs),intent(inout) :: indcol  !! an integer array of length `npairs`. on input `indcol`
                                                         !! must contain the column indices of the non-zero elements of
                                                         !! `a`. on output `indcol` is permuted so that the corresponding
                                                         !! row indices are in non-decreasing order. the row indices
                                                         !! can be recovered from the array `ipntr`.
      integer,dimension(n),intent(out) :: ngrp   !! specifies the partition of the columns of `a`.
                                                 !! column `jcol` belongs to group `ngrp(jcol)`.
      integer,dimension(m+1),intent(out) :: ipntr   !! an integer output array of length `m + 1` which
                                                    !! specifies the locations of the column indices in `indcol`.
                                                    !! the column indices for row `i` are
                                                    !! `indcol(k), k = ipntr(i),...,ipntr(i+1)-1`.
                                                    !! note that `ipntr(m+1)-1` is then the number of non-zero
                                                    !! elements of the matrix `a`.
      integer,dimension(n+1),intent(out) :: jpntr   !! jpntr is an integer output array of length n + 1 which
                                                    !! specifies the locations of the row indices in indrow.
                                                    !! the row indices for column j are
                                                    !! `indrow(k), k = jpntr(j),...,jpntr(j+1)-1`.
                                                    !! note that `jpntr(n+1)-1` is then the number of non-zero
                                                    !! elements of the matrix `a`.

    integer,dimension(max(m,6*n)) :: iwa     !! an integer work array
    integer :: i , ir , j , jp , k , maxclq , nnz , numgrp

    !  check the input data.

    Info = 0
    if ( m<1 .or. n<1 .or. Npairs<1 ) return
    do k = 1 , Npairs
        Info = -k
        if ( Indrow(k)<1 .or. Indrow(k)>m .or. Indcol(k)<1 .or. Indcol(k)>n ) return
    enddo
    Info = 1

    !  sort the data structure by columns.

    call srtdat(n,Npairs,Indrow,Indcol,Jpntr,Iwa)

    !  compress the data and determine the number of
    !  non-zero elements of a.

    do i = 1 , m
        Iwa(i) = 0
    enddo
    nnz = 1
    do j = 1 , n
        k = nnz
        do jp = Jpntr(j) , Jpntr(j+1) - 1
            ir = Indrow(jp)
            if ( Iwa(ir)/=j ) then
                Indrow(nnz) = ir
                nnz = nnz + 1
                Iwa(ir) = j
            endif
        enddo
        Jpntr(j) = k
    enddo
    Jpntr(n+1) = nnz

    !  extend the data structure to rows.

    call setr(m,n,Indrow,Jpntr,Indcol,Ipntr,Iwa)

    !  determine a lower bound for the number of groups.

    Mingrp = 0
    do i = 1 , m
        Mingrp = max(Mingrp,Ipntr(i+1)-Ipntr(i))
    enddo

    !  determine the degree sequence for the intersection
    !  graph of the columns of a.

    call degr(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(5*n+1),Iwa(n+1))

    !  color the intersection graph of the columns of a
    !  with the smallest-last (sl) ordering.

    call slo(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(5*n+1),Iwa(4*n+1),maxclq,&
                Iwa(1),Iwa(n+1),Iwa(2*n+1),Iwa(3*n+1))
    call seq(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(4*n+1),Ngrp,Maxgrp,      &
                Iwa(n+1))
    Mingrp = max(Mingrp,maxclq)

    !  exit if the smallest-last ordering is optimal.

    if ( Maxgrp==Mingrp ) return

    !  color the intersection graph of the columns of a
    !  with the incidence-degree (id) ordering.

    call ido(m,n,Indrow,Jpntr,Indcol,Ipntr,Iwa(5*n+1),Iwa(4*n+1),     &
                maxclq,Iwa(1),Iwa(n+1),Iwa(2*n+1),Iwa(3*n+1))
    call seq(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(4*n+1),Iwa(1),numgrp,    &
                Iwa(n+1))
    Mingrp = max(Mingrp,maxclq)

    !  retain the better of the two orderings so far.

    if ( numgrp<Maxgrp ) then
        Maxgrp = numgrp
        do j = 1 , n
            Ngrp(j) = Iwa(j)
        enddo

        !  exit if the incidence-degree ordering is optimal.

        if ( Maxgrp==Mingrp ) return
    endif

    !  color the intersection graph of the columns of a
    !  with the largest-first (lf) ordering.

    call numsrt(n,n-1,Iwa(5*n+1),-1,Iwa(4*n+1),Iwa(2*n+1),Iwa(n+1))
    call seq(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(4*n+1),Iwa(1),numgrp,    &
                Iwa(n+1))

    !  retain the best of the three orderings and exit.

    if ( numgrp<Maxgrp ) then
        Maxgrp = numgrp
        do j = 1 , n
            Ngrp(j) = Iwa(j)
        enddo
    endif

    end subroutine dsm
!*******************************************************************************

!*******************************************************************************
!>
!  Given the sparsity pattern of an `m` by `n` matrix `a`,
!  this subroutine determines the degree sequence for
!  the intersection graph of the columns of `a`.
!
!  In graph-theory terminology, the intersection graph of
!  the columns of `a` is the loopless graph `g` with vertices
!  `a(j), j = 1,2,...,n` where `a(j)` is the `j`-th column of `a`
!  and with edge `(a(i),a(j))` if and only if columns `i` and `j`
!  have a non-zero in the same row position.
!
!@note The value of `m` is not needed by `degr` and is
!      therefore not present in the subroutine statement.

    subroutine degr(n,Indrow,Jpntr,Indcol,Ipntr,Ndeg,Iwa)

    implicit none

    integer,intent(in)                :: n       !! a positive integer input variable set to the number
                                                 !! of columns of `a`.
    integer,dimension(*),intent(in)   :: indrow  !! an integer input array which contains the row
                                                 !! indices for the non-zeroes in the matrix `a`.
    integer,dimension(n+1),intent(in) :: jpntr   !! an integer input array of length `n + 1` which
                                                 !! specifies the locations of the row indices in `indrow`.
                                                 !! the row indices for column `j` are
                                                 !! `indrow(k), k = jpntr(j),...,jpntr(j+1)-1`.
                                                 !! **note** that `jpntr(n+1)-1` is then the number of non-zero
                                                 !! elements of the matrix `a`.
    integer,dimension(*),intent(in)   :: indcol  !! an integer input array which contains the
                                                 !! column indices for the non-zeroes in the matrix `a`.
    integer,dimension(*),intent(in)   :: ipntr   !! an integer input array of length `m + 1` which
                                                 !! specifies the locations of the column indices in `indcol`.
                                                 !! the column indices for row `i` are
                                                 !! `indcol(k), k = ipntr(i),...,ipntr(i+1)-1`.
                                                 !! **note** that `ipntr(m+1)-1` is then the number of non-zero
                                                 !! elements of the matrix `a`.
    integer,dimension(n),intent(out)   :: ndeg   !! an integer output array of length `n` which
                                                 !! specifies the degree sequence. the degree of the
                                                 !! `j`-th column of `a` is `ndeg(j)`.
    integer,dimension(n)               :: iwa    !! an integer work array of length `n`

    integer :: ic , ip , ir , jcol , jp

    ! initialization block.

    do jp = 1 , n
        Ndeg(jp) = 0
        Iwa(jp) = 0
    enddo

    ! compute the degree sequence by determining the contributions
    ! to the degrees from the current(jcol) column and further
    ! columns which have not yet been considered.

    do jcol = 2 , n
        Iwa(jcol) = n

        ! determine all positions (ir,jcol) which correspond
        ! to non-zeroes in the matrix.

        do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
            ir = Indrow(jp)

            ! for each row ir, determine all positions (ir,ic)
            ! which correspond to non-zeroes in the matrix.

            do ip = Ipntr(ir) , Ipntr(ir+1) - 1
                ic = Indcol(ip)

                ! array iwa marks columns which have contributed to
                ! the degree count of column jcol. update the degree
                ! counts of these columns as well as column jcol.

                if ( Iwa(ic)<jcol ) then
                    Iwa(ic) = jcol
                    Ndeg(ic) = Ndeg(ic) + 1
                    Ndeg(jcol) = Ndeg(jcol) + 1
                endif

            enddo
        enddo
    enddo

    end subroutine degr
!*******************************************************************************

!*******************************************************************************
!>
!  given the sparsity pattern of an `m` by `n` matrix `a`, this
!  subroutine determines an incidence-degree ordering of the
!  columns of `a`.
!
!  the incidence-degree ordering is defined for the loopless
!  graph `g` with vertices `a(j), j = 1,2,...,n` where `a(j)` is the
!  `j`-th column of `a` and with edge `(a(i),a(j))` if and only if
!  columns `i` and `j` have a non-zero in the same row position.
!
!  the incidence-degree ordering is determined recursively by
!  letting `list(k), k = 1,...,n` be a column with maximal
!  incidence to the subgraph spanned by the ordered columns.
!  among all the columns of maximal incidence, `ido` chooses a
!  column of maximal degree.

    subroutine ido(m,n,Indrow,Jpntr,Indcol,Ipntr,Ndeg,List,Maxclq, &
                   Iwa1,Iwa2,Iwa3,Iwa4)

    implicit none

    integer,intent(in)       :: m         !! a positive integer input variable set to the number
                                          !! of rows of `a`.
    integer,intent(in)       :: n         !! a positive integer input variable set to the number
                                          !! of columns of `a`.
    integer,intent(out)      :: Maxclq    !! an integer output variable set to the size
                                          !! of the largest clique found during the ordering.
    integer,dimension(*),intent(in) :: Indrow    !! an integer input array which contains the row
                                                 !! indices for the non-zeroes in the matrix `a`.
    integer,dimension(n+1),intent(in) :: Jpntr   !! an integer input array of length `n + 1` which
                                                 !! specifies the locations of the row indices in `indrow`.
                                                 !! the row indices for column `j` are
                                                 !! `indrow(k), k = jpntr(j),...,jpntr(j+1)-1`.
                                                 !! **note** that `jpntr(n+1)-1` is then the number of non-zero
                                                 !! elements of the matrix `a`.
    integer,dimension(*),intent(in)     :: Indcol   !! an integer input array which contains the
                                                    !! column indices for the non-zeroes in the matrix `a`.
    integer,dimension(m+1),intent(in)   :: Ipntr    !! an integer input array of length `m + 1` which
                                                    !! specifies the locations of the column indices in `indcol`.
                                                    !! the column indices for row `i` are
                                                    !! `indcol(k), k = ipntr(i),...,ipntr(i+1)-1`.
                                                    !! **note** that `ipntr(m+1)-1` is then the number of non-zero
                                                    !! elements of the matrix `a`.
    integer,dimension(n),intent(in)     :: Ndeg     !! an integer input array of length `n` which specifies
                                                    !! the degree sequence. the degree of the `j`-th column
                                                    !! of `a` is `ndeg(j)`.
    integer,dimension(n),intent(out)    :: List     !! an integer output array of length `n` which specifies
                                                    !! the incidence-degree ordering of the columns of `a`. the `j`-th
                                                    !! column in this order is `list(j)`.
    integer,dimension(0:n-1) :: Iwa1      !! integer work array of length `n`.
    integer,dimension(n)     :: Iwa2      !! integer work array of length `n`.
    integer,dimension(n)     :: Iwa3      !! integer work array of length `n`.
    integer,dimension(n)     :: Iwa4      !! integer work array of length `n`.

    integer :: ic , ip , ir , jcol , jp , maxinc , maxlst , ncomp , &
               numinc , numlst , numord , numwgt

    ! sort the degree sequence.

    call numsrt(n,n-1,Ndeg,-1,Iwa4,Iwa2,Iwa3)

    ! initialization block.
    !
    ! create a doubly-linked list to access the incidences of the
    ! columns. the pointers for the linked list are as follows.
    !
    ! each un-ordered column ic is in a list (the incidence list)
    ! of columns with the same incidence.
    !
    ! iwa1(numinc) is the first column in the numinc list
    ! unless iwa1(numinc) = 0. in this case there are
    ! no columns in the numinc list.
    !
    ! iwa2(ic) is the column before ic in the incidence list
    ! unless iwa2(ic) = 0. in this case ic is the first
    ! column in this incidence list.
    !
    ! iwa3(ic) is the column after ic in the incidence list
    ! unless iwa3(ic) = 0. in this case ic is the last
    ! column in this incidence list.
    !
    ! if ic is an un-ordered column, then list(ic) is the
    ! incidence of ic to the graph induced by the ordered
    ! columns. if jcol is an ordered column, then list(jcol)
    ! is the incidence-degree order of column jcol.

    maxinc = 0
    do jp = n , 1 , -1
        ic = Iwa4(jp)
        Iwa1(n-jp) = 0
        Iwa2(ic) = 0
        Iwa3(ic) = Iwa1(0)
        if ( Iwa1(0)>0 ) Iwa2(Iwa1(0)) = ic
        Iwa1(0) = ic
        Iwa4(jp) = 0
        List(jp) = 0
    enddo

    ! DETERMINE THE MAXIMAL SEARCH LENGTH FOR THE LIST
    ! OF COLUMNS OF MAXIMAL INCIDENCE.

    maxlst = 0
    do ir = 1 , m
        maxlst = maxlst + (Ipntr(ir+1)-Ipntr(ir))**2
    enddo
    maxlst = maxlst/n
    Maxclq = 0
    numord = 1

    ! BEGINNING OF ITERATION LOOP.

    ! UPDATE THE SIZE OF THE LARGEST CLIQUE
    ! FOUND DURING THE ORDERING.

100 if ( maxinc==0 ) ncomp = 0
    ncomp = ncomp + 1
    if ( maxinc+1==ncomp ) Maxclq = max(Maxclq,ncomp)

    ! CHOOSE A COLUMN JCOL OF MAXIMAL DEGREE AMONG THE
    ! COLUMNS OF MAXIMAL INCIDENCE MAXINC.

200 jp = Iwa1(maxinc)
    if ( jp>0 ) then
        numwgt = -1
        do numlst = 1 , maxlst
            if ( Ndeg(jp)>numwgt ) then
                numwgt = Ndeg(jp)
                jcol = jp
            endif
            jp = Iwa3(jp)
            if ( jp<=0 ) exit
        enddo
        List(jcol) = numord
        numord = numord + 1

        ! TERMINATION TEST.

        if ( numord>n ) then

            ! INVERT THE ARRAY LIST.

            do jcol = 1 , n
                Iwa2(List(jcol)) = jcol
            enddo
            do jp = 1 , n
                List(jp) = Iwa2(jp)
            enddo

        else

            ! DELETE COLUMN JCOL FROM THE MAXINC LIST.

            if ( Iwa2(jcol)==0 ) then
                Iwa1(maxinc) = Iwa3(jcol)
            else
                Iwa3(Iwa2(jcol)) = Iwa3(jcol)
            endif
            if ( Iwa3(jcol)>0 ) Iwa2(Iwa3(jcol)) = Iwa2(jcol)

            ! FIND ALL COLUMNS ADJACENT TO COLUMN JCOL.

            Iwa4(jcol) = n

            ! DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
            ! TO NON-ZEROES IN THE MATRIX.

            do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
                ir = Indrow(jp)

                ! FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
                ! WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.

                do ip = Ipntr(ir) , Ipntr(ir+1) - 1
                    ic = Indcol(ip)

                    ! ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO
                    ! COLUMN JCOL.

                    if ( Iwa4(ic)<numord ) then
                        Iwa4(ic) = numord

                        ! UPDATE THE POINTERS TO THE CURRENT INCIDENCE LISTS.

                        numinc = List(ic)
                        List(ic) = List(ic) + 1
                        maxinc = max(maxinc,List(ic))

                        ! DELETE COLUMN IC FROM THE NUMINC LIST.

                        if ( Iwa2(ic)==0 ) then
                            Iwa1(numinc) = Iwa3(ic)
                        else
                            Iwa3(Iwa2(ic)) = Iwa3(ic)
                        endif
                        if ( Iwa3(ic)>0 ) Iwa2(Iwa3(ic)) = Iwa2(ic)

                        ! ADD COLUMN IC TO THE NUMINC+1 LIST.

                        Iwa2(ic) = 0
                        Iwa3(ic) = Iwa1(numinc+1)
                        if ( Iwa1(numinc+1)>0 ) Iwa2(Iwa1(numinc+1)) = ic
                        Iwa1(numinc+1) = ic
                    endif
                enddo
            enddo

            ! END OF ITERATION LOOP.

            goto 100
        endif
    else
        maxinc = maxinc - 1
        goto 200
    endif

    end subroutine ido
!*******************************************************************************

!*******************************************************************************
!>
!  Given a sequence of integers, this subroutine groups
!  together those indices with the same sequence value
!  and, optionally, sorts the sequence into either
!  ascending or descending order.
!
!  The sequence of integers is defined by the array `num`,
!  and it is assumed that the integers are each from the set
!  `0,1,...,nmax`. on output the indices `k` such that `num(k) = l`
!  for any `l = 0,1,...,nmax` can be obtained from the arrays
!  last and next as follows.
!```fortran
!  k = last(l)
!  while (k /= 0) k = next(k)
!```
!  Optionally, the subroutine produces an array index so that
!  the sequence `num(index(i)), i = 1,2,...,n` is sorted.

    subroutine numsrt(n,Nmax,Num,Mode,Index,Last,Next)

    implicit none

    integer,intent(in)        :: n      !! a positive integer input variable.
    integer                   :: Nmax   !! a positive integer input variable.
    integer                   :: Mode   !! an integer input variable. the sequence `num` is
                                        !! sorted in ascending order if `mode` is positive and in
                                        !! descending order if `mode` is negative. if `mode` is 0,
                                        !! no sorting is done.
    integer,dimension(n)      :: Num    !! an input array of length `n` which contains the
                                        !! sequence of integers to be grouped and sorted. it
                                        !! is assumed that the integers are each from the set
                                        !! `0,1,...,nmax`.
    integer,dimension(n)      :: Index  !! an integer output array of length `n` set so
                                        !! that the sequence
                                        !! `num(index(i)), i = 1,2,...,n`
                                        !! is sorted according to the setting of mode.
                                        !! if `mode` is 0, `index` is not referenced.
    integer,dimension(0:Nmax) :: Last   !! an integer output array of length `nmax + 1`. the
                                        !! index of `num` for the last occurrence of `l` is `last(l)`
                                        !! for any `l = 0,1,...,nmax` unless `last(l) = 0`. in
                                        !! this case `l` does not appear in `num`.
    integer,dimension(n)      :: Next   !! an integer output array of length `n`. if
                                        !! `num(k) = l`, then the index of `num` for the previous
                                        !! occurrence of `l` is `next(k)` for any `l = 0,1,...,nmax`
                                        !! unless `next(k) = 0`. in this case there is no previous
                                        !! occurrence of `l` in `num`.

    integer :: i , j , jinc , jl , ju , k , l

    ! determine the arrays next and last.

    do i = 0 , Nmax
        Last(i) = 0
    enddo

    do k = 1 , n
        l = Num(k)
        Next(k) = Last(l)
        Last(l) = k
    enddo

    if ( Mode/=0 ) then

        ! store the pointers to the sorted array in index.
        i = 1
        if ( Mode>0 ) then
            jl = 0
            ju = Nmax
            jinc = 1
        else
            jl = Nmax
            ju = 0
            jinc = -1
        endif
        do j = jl , ju , jinc
            k = Last(j)
            do
                if (k==0) exit
                Index(i) = k
                i = i + 1
                k = Next(k)
            enddo
        enddo

    end if

    end subroutine numsrt
!*******************************************************************************

!*******************************************************************************
!>
!  given the sparsity pattern of an `m` by `n` matrix `a`, this
!  subroutine determines a consistent partition of the
!  columns of `a` by a sequential algorithm.
!
!  a consistent partition is defined in terms of the loopless
!  graph `g` with vertices `a(j), j = 1,2,...,n` where `a(j)` is the
!  `j`-th column of `a` and with edge `(a(i),a(j))` if and only if
!  columns `i` and `j` have a non-zero in the same row position.
!
!  a partition of the columns of a into groups is consistent
!  if the columns in any group are not adjacent in the graph `g`.
!  in graph-theory terminology, a consistent partition of the
!  columns of a corresponds to a coloring of the graph `g`.
!
!  the subroutine examines the columns in the order specified
!  by the array list, and assigns the current column to the
!  group with the smallest possible number.
!
!  note that the value of `m` is not needed by `seq` and is
!  therefore not present in the subroutine statement.

    subroutine seq(n,Indrow,Jpntr,Indcol,Ipntr,List,Ngrp,Maxgrp,Iwa)

    implicit none

    integer                :: n         !! a positive integer input variable set to the number
                                        !! of columns of `a`.
    integer                :: Maxgrp    !! an integer output variable which specifies the
                                        !! number of groups in the partition of the columns of `a`.
    integer,dimension(*)   :: Indrow    !! an integer input array which contains the row
                                        !! indices for the non-zeroes in the matrix `a`.
    integer,dimension(n+1) :: Jpntr     !! an integer input array of length `n + 1` which
                                        !! specifies the locations of the row indices in `indrow`.
                                        !! the row indices for column `j` are
                                        !! `indrow(k), k = jpntr(j),...,jpntr(j+1)-1`.
                                        !! **note** that `jpntr(n+1)-1` is then the number of non-zero
                                        !! elements of the matrix `a`.
    integer,dimension(*)   :: Indcol    !! an integer input array which contains the
                                        !! column indices for the non-zeroes in the matrix `a`.
    integer,dimension(*)   :: Ipntr     !! an integer input array of length `m + 1` which
                                        !! specifies the locations of the column indices in `indcol`.
                                        !! the column indices for row `i` are
                                        !! `indcol(k), k = ipntr(i),...,ipntr(i+1)-1`.
                                        !! **note** that `ipntr(m+1)-1` is then the number of non-zero
                                        !! elements of the matrix `a`.
    integer,dimension(n)   :: List      !! an integer input array of length `n` which specifies
                                        !! the order to be used by the sequential algorithm.
                                        !! the `j`-th column in this order is `list(j)`.
    integer,dimension(n)   :: Ngrp      !! an integer output array of length `n` which specifies
                                        !! the partition of the columns of `a`. column `jcol` belongs
                                        !! to group `ngrp(jcol)`.
    integer,dimension(n)   :: Iwa       !! an integer work array of length `n`

    integer :: ic , ip , ir , j , jcol , jp

    ! initialization block.

    Maxgrp = 0
    do jp = 1 , n
        Ngrp(jp) = n
        Iwa(jp) = 0
    enddo

    ! beginning of iteration loop.
    do j = 1 , n

        jcol = List(j)

        ! find all columns adjacent to column jcol.
        !
        ! determine all positions (ir,jcol) which correspond
        ! to non-zeroes in the matrix.
        do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
            ir = Indrow(jp)
            ! for each row ir, determine all positions (ir,ic)
            ! which correspond to non-zeroes in the matrix.
            do ip = Ipntr(ir) , Ipntr(ir+1) - 1
                ic = Indcol(ip)
                ! array iwa marks the group numbers of the
                ! columns which are adjacent to column jcol.
                Iwa(Ngrp(ic)) = j
            enddo
        enddo

        ! assign the smallest un-marked group number to jcol.
        do jp = 1 , Maxgrp
            if ( Iwa(jp)/=j ) then
                Maxgrp = Maxgrp - 1
                exit
            end if
        enddo
        Maxgrp = Maxgrp + 1
        Ngrp(jcol) = jp

    enddo
    ! end of iteration loop.

    end subroutine seq
!*******************************************************************************

!*******************************************************************************
!>
!  given a column-oriented definition of the sparsity pattern
!  of an `m` by `n` matrix `a`, this subroutine determines a
!  row-oriented definition of the sparsity pattern of `a`.
!
!  on input the column-oriented definition is specified by
!  the arrays `indrow` and `jpntr`. on output the row-oriented
!  definition is specified by the arrays `indcol` and `ipntr`.

    subroutine setr(m,n,Indrow,Jpntr,Indcol,Ipntr,Iwa)

    implicit none

    integer,intent(in)     :: m       !! a positive integer input variable set to the number
                                      !! of rows of `a`.
    integer,intent(in)     :: n       !! a positive integer input variable set to the number
                                      !! of columns of `a`.
    integer,dimension(*)   :: Indrow  !! an integer input array which contains the row
                                      !! indices for the non-zeroes in the matrix `a`.
    integer,dimension(n+1) :: Jpntr   !! an integer input array of length `n + 1` which
                                      !! specifies the locations of the row indices in `indrow`.
                                      !! the row indices for column `j` are
                                      !! `indrow(k), k = jpntr(j),...,jpntr(j+1)-1`.
                                      !! **note** that `jpntr(n+1)-1` is then the number of non-zero
                                      !! elements of the matrix `a`.
    integer,dimension(*)   :: Indcol  !! an integer output array which contains the
                                      !! column indices for the non-zeroes in the matrix `a`.
    integer,dimension(m+1) :: Ipntr   !! an integer output array of length `m + 1` which
                                      !! specifies the locations of the column indices in `indcol`.
                                      !! the column indices for row `i` are
                                      !! `indcol(k), k = ipntr(i),...,ipntr(i+1)-1`.
                                      !! **note** that `ipntr(1)` is set to 1 and that `ipntr(m+1)-1` is
                                      !! then the number of non-zero elements of the matrix `a`.
    integer,dimension(m)   :: Iwa     !! an integer work array of length `m`.

    integer :: ir , jcol , jp

    ! store in array iwa the counts of non-zeroes in the rows.

    do ir = 1 , m
        Iwa(ir) = 0
    enddo
    do jp = 1 , Jpntr(n+1) - 1
        Iwa(Indrow(jp)) = Iwa(Indrow(jp)) + 1
    enddo

    ! set pointers to the start of the rows in indcol.

    Ipntr(1) = 1
    do ir = 1 , m
        Ipntr(ir+1) = Ipntr(ir) + Iwa(ir)
        Iwa(ir) = Ipntr(ir)
    enddo

    ! fill indcol.

    do jcol = 1 , n
        do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
            ir = Indrow(jp)
            Indcol(Iwa(ir)) = jcol
            Iwa(ir) = Iwa(ir) + 1
        enddo
    enddo

    end subroutine setr
!*******************************************************************************

!*******************************************************************************
!>
!  given the sparsity pattern of an `m` by `n` matrix `a`, this
!  subroutine determines the smallest-last ordering of the
!  columns of `a`.
!
!  the smallest-last ordering is defined for the loopless
!  graph `g` with vertices `a(j), j = 1,2,...,n` where `a(j)` is the
!  `j`-th column of `a` and with edge `(a(i),a(j))` if and only if
!  columns `i` and `j` have a non-zero in the same row position.
!
!  the smallest-last ordering is determined recursively by
!  letting `list(k), k = n,...,1` be a column with least degree
!  in the subgraph spanned by the un-ordered columns.
!
!  note that the value of `m` is not needed by `slo` and is
!  therefore not present in the subroutine statement.

    subroutine slo(n,Indrow,Jpntr,Indcol,Ipntr,Ndeg,List,Maxclq,Iwa1, &
                   Iwa2,Iwa3,Iwa4)

    implicit none

    integer                  :: n         !! a positive integer input variable set to the number
                                          !! of columns of `a`.
    integer                  :: Maxclq    !! an integer output variable set to the size
                                          !! of the largest clique found during the ordering.
    integer,dimension(*)     :: Indrow    !! an integer input array which contains the row
                                          !! indices for the non-zeroes in the matrix `a`.
    integer,dimension(n+1)   :: Jpntr     !! an integer input array of length `n + 1` which
                                          !! specifies the locations of the row indices in `indrow`.
                                          !! the row indices for column `j` are
                                          !! `indrow(k), k = jpntr(j),...,jpntr(j+1)-1`.
                                          !! **note** that `jpntr(n+1)-1` is then the number of non-zero
                                          !! elements of the matrix `a`.
    integer,dimension(*)     :: Indcol    !! an integer input array which contains the
                                          !! column indices for the non-zeroes in the matrix `a`.
    integer,dimension(*)     :: Ipntr     !! an integer input array of length `m + 1` which
                                          !! specifies the locations of the column indices in `indcol`.
                                          !! the column indices for row `i` are
                                          !! `indcol(k), k = ipntr(i),...,ipntr(i+1)-1`.
                                          !! **note** that `ipntr(m+1)-1` is then the number of non-zero
                                          !! elements of the matrix `a`.
    integer,dimension(n)     :: Ndeg      !! an integer input array of length `n` which specifies
                                          !! the degree sequence. the degree of the `j`-th column
                                          !! of `a` is `ndeg(j)`.
    integer,dimension(n)     :: List      !! an integer output array of length `n` which specifies
                                          !! the smallest-last ordering of the columns of `a`. the `j`-th
                                          !! column in this order is `list(j)`.
    integer,dimension(0:n-1) :: Iwa1      !! integer work array of length `n`
    integer,dimension(n)     :: Iwa2      !! integer work array of length `n`
    integer,dimension(n)     :: Iwa3      !! integer work array of length `n`
    integer,dimension(n)     :: Iwa4      !! integer work array of length `n`

    integer :: ic , ip , ir , jcol , jp , mindeg , numdeg , numord

    ! INITIALIZATION BLOCK.

    mindeg = n
    do jp = 1 , n
        Iwa1(jp-1) = 0
        Iwa4(jp) = n
        List(jp) = Ndeg(jp)
        mindeg = min(mindeg,Ndeg(jp))
    enddo

    ! CREATE A DOUBLY-LINKED LIST TO ACCESS THE DEGREES OF THE
    ! COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS.
    !
    ! EACH UN-ORDERED COLUMN IC IS IN A LIST (THE DEGREE LIST)
    ! OF COLUMNS WITH THE SAME DEGREE.
    !
    ! IWA1(NUMDEG) IS THE FIRST COLUMN IN THE NUMDEG LIST
    ! UNLESS IWA1(NUMDEG) = 0. IN THIS CASE THERE ARE
    ! NO COLUMNS IN THE NUMDEG LIST.
    !
    ! IWA2(IC) IS THE COLUMN BEFORE IC IN THE DEGREE LIST
    ! UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST
    ! COLUMN IN THIS DEGREE LIST.
    !
    ! IWA3(IC) IS THE COLUMN AFTER IC IN THE DEGREE LIST
    ! UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST
    ! COLUMN IN THIS DEGREE LIST.
    !
    ! IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE
    ! DEGREE OF IC IN THE GRAPH INDUCED BY THE UN-ORDERED
    ! COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL)
    ! IS THE SMALLEST-LAST ORDER OF COLUMN JCOL.

    do jp = 1 , n
        numdeg = Ndeg(jp)
        Iwa2(jp) = 0
        Iwa3(jp) = Iwa1(numdeg)
        if ( Iwa1(numdeg)>0 ) Iwa2(Iwa1(numdeg)) = jp
        Iwa1(numdeg) = jp
    enddo
    Maxclq = 0
    numord = n

    !  BEGINNING OF ITERATION LOOP.
    !
    !
    ! MARK THE SIZE OF THE LARGEST CLIQUE
    ! FOUND DURING THE ORDERING.

100 if ( mindeg+1==numord .and. Maxclq==0 ) Maxclq = numord

    ! CHOOSE A COLUMN JCOL OF MINIMAL DEGREE MINDEG.

200 jcol = Iwa1(mindeg)
    if ( jcol>0 ) then
        List(jcol) = numord
        numord = numord - 1

        ! TERMINATION TEST.

        if ( numord==0 ) then

            ! INVERT THE ARRAY LIST.

            do jcol = 1 , n
                Iwa2(List(jcol)) = jcol
            enddo
            do jp = 1 , n
                List(jp) = Iwa2(jp)
            enddo

        else

            ! DELETE COLUMN JCOL FROM THE MINDEG LIST.

            Iwa1(mindeg) = Iwa3(jcol)
            if ( Iwa3(jcol)>0 ) Iwa2(Iwa3(jcol)) = 0

            ! FIND ALL COLUMNS ADJACENT TO COLUMN JCOL.

            Iwa4(jcol) = 0

            ! DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
            ! TO NON-ZEROES IN THE MATRIX.

            do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
                ir = Indrow(jp)

                ! FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
                ! WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.

                do ip = Ipntr(ir) , Ipntr(ir+1) - 1
                    ic = Indcol(ip)

                    ! ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO
                    ! COLUMN JCOL.

                    if ( Iwa4(ic)>numord ) then
                        Iwa4(ic) = numord

                        ! UPDATE THE POINTERS TO THE CURRENT DEGREE LISTS.

                        numdeg = List(ic)
                        List(ic) = List(ic) - 1
                        mindeg = min(mindeg,List(ic))

                        ! DELETE COLUMN IC FROM THE NUMDEG LIST.

                        if ( Iwa2(ic)==0 ) then
                            Iwa1(numdeg) = Iwa3(ic)
                        else
                            Iwa3(Iwa2(ic)) = Iwa3(ic)
                        endif
                        if ( Iwa3(ic)>0 ) Iwa2(Iwa3(ic)) = Iwa2(ic)

                        ! ADD COLUMN IC TO THE NUMDEG-1 LIST.

                        Iwa2(ic) = 0
                        Iwa3(ic) = Iwa1(numdeg-1)
                        if ( Iwa1(numdeg-1)>0 ) Iwa2(Iwa1(numdeg-1)) = ic
                        Iwa1(numdeg-1) = ic
                    endif
                enddo
            enddo

            ! END OF ITERATION LOOP.

            goto 100
        endif
    else
        mindeg = mindeg + 1
        goto 200
    endif

    end subroutine slo
!*******************************************************************************

!*******************************************************************************
!>
!  given the non-zero elements of an `m` by `n` matrix `a` in
!  arbitrary order as specified by their row and column
!  indices, this subroutine permutes these elements so
!  that their column indices are in non-decreasing order.
!
!  on input it is assumed that the elements are specified in
!
!  `indrow(k),indcol(k), k = 1,...,nnz`.
!
!  on output the elements are permuted so that `indcol` is
!  in non-decreasing order. in addition, the array `jpntr`
!  is set so that the row indices for column `j` are
!
!  `indrow(k), k = jpntr(j),...,jpntr(j+1)-1`.
!
!  note that the value of `m` is not needed by srtdat and is
!  therefore not present in the subroutine statement.

    subroutine srtdat(n,Nnz,Indrow,Indcol,Jpntr,Iwa)

    implicit none

    integer                :: n       !! a positive integer input variable set to the number
                                      !! of columns of `a`.
    integer                :: Nnz     !! a positive integer input variable set to the number
                                      !! of non-zero elements of `a`.
    integer,dimension(Nnz) :: Indrow  !! an integer array of length `nnz`. on input `indrow`
                                      !! must contain the row indices of the non-zero elements of `a`.
                                      !! on output `indrow` is permuted so that the corresponding
                                      !! column indices of `indcol` are in non-decreasing order.
    integer,dimension(Nnz) :: Indcol  !! an integer array of length `nnz`. on input `indcol`
                                      !! must contain the column indices of the non-zero elements
                                      !! of `a`. on output `indcol` is permuted so that these indices
                                      !! are in non-decreasing order.
    integer,dimension(n+1) :: Jpntr   !! an integer output array of length `n + 1` which
                                      !! specifies the locations of the row indices in the output
                                      !! `indrow`. the row indices for column `j` are
                                      !! `indrow(k), k = jpntr(j),...,jpntr(j+1)-1`.
                                      !! **note** that `jpntr(1)` is set to 1 and that `jpntr(n+1)-1`
                                      !! is then `nnz`.
    integer,dimension(n)   :: Iwa     !! an integer work array of length `n`.

    integer :: i , j , k , l

    ! store in array iwa the counts of non-zeroes in the columns.

    do j = 1 , n
        Iwa(j) = 0
    enddo
    do k = 1 , Nnz
        Iwa(Indcol(k)) = Iwa(Indcol(k)) + 1
    enddo

    ! set pointers to the start of the columns in indrow.

    Jpntr(1) = 1
    do j = 1 , n
        Jpntr(j+1) = Jpntr(j) + Iwa(j)
        Iwa(j) = Jpntr(j)
    enddo
    k = 1

    ! begin in-place sort.

    do
        j = Indcol(k)
        if ( k>=Jpntr(j) ) then
            ! current element is in position. now examine the
            ! next element or the first un-sorted element in
            ! the j-th group.
            k = max(k+1,Iwa(j))
        else
            ! current element is not in position. place element
            ! in position and make the displaced element the
            ! current element.
            l = Iwa(j)
            Iwa(j) = Iwa(j) + 1
            i = Indrow(k)
            Indrow(k) = Indrow(l)
            Indcol(k) = Indcol(l)
            Indrow(l) = i
            Indcol(l) = j
        endif
        if ( k>Nnz ) exit
    end do

    end subroutine srtdat
!*******************************************************************************

!*******************************************************************************
!>
!  Given a consistent partition of the columns of an `m` by `n`
!  jacobian matrix into groups, this subroutine computes
!  approximations to those columns in a given group.  the
!  approximations are stored into either a column-oriented
!  or a row-oriented pattern.
!
!  a partition is consistent if the columns in any group
!  do not have a non-zero in the same row position.
!
!  approximations to the columns of the jacobian matrix in a
!  given group can be obtained by specifying a difference
!  parameter array `d` with `d(jcol)` non-zero if and only if
!  `jcol` is a column in the group, and an approximation to
!  `jac*d` where `jac` denotes the jacobian matrix of a mapping f.
!
!  `d` can be defined with the following segment of code.
!```fortran
!  do jcol = 1, n
!  d(jcol) = 0.0
!  if (ngrp(jcol) == numgrp) d(jcol) = eta(jcol)
!  end do
!```
!  in the above code `numgrp` is the given group number,
!  `ngrp(jcol)` is the group number of column `jcol`, and
!  `eta(jcol)` is the difference parameter used to
!  approximate column `jcol` of the jacobian matrix.
!  suitable values for the array `eta` must be provided.
!
!  as mentioned above, an approximation to `jac*d` must
!  also be provided. for example, the approximation
!```fortran
!  f(x+d) - f(x)
!```
!  corresponds to the forward difference formula at `x`.

    subroutine fdjs(m,n,Col,Ind,Npntr,Ngrp,Numgrp,d,Fjacd,Fjac)

    implicit none

    integer,intent(in)    :: m        !! a positive integer input variable set to the number
                                      !! of rows of the jacobian matrix.
    integer,intent(in)    :: n        !! a positive integer input variable set to the number
                                      !! of columns of the jacobian matrix.
    integer               :: Numgrp   !! a positive integer input variable set to a group
                                      !! number in the partition. the columns of the jacobian
                                      !! matrix in this group are to be estimated on this call.
    integer,dimension(*)  :: Ind      !! an integer input array which contains the row
                                      !! indices for the non-zeroes in the jacobian matrix
                                      !! if `col` is true, and contains the column indices for
                                      !! the non-zeroes in the jacobian matrix if `col` is false.
    integer,dimension(*)  :: Npntr    !! an integer input array which specifies the
                                      !! locations of the row indices in `ind` if `col` is true, and
                                      !! specifies the locations of the column indices in `ind` if
                                      !! `col` is false. if `col` is true, the indices for column `j` are
                                      !!       `ind(k), k = npntr(j),...,npntr(j+1)-1`.
                                      !! if `col` is false, the indices for row `i` are
                                      !!       `ind(k), k = npntr(i),...,npntr(i+1)-1`.
                                      !! ***Note*** that `npntr(n+1)-1` if `col` is true, or `npntr(m+1)-1`
                                      !! if `col` is false, is then the number of non-zero elements
                                      !! of the jacobian matrix.
    integer,dimension(n)  :: Ngrp     !! an integer input array of length `n` which specifies
                                      !! the partition of the columns of the jacobian matrix.
                                      !! column `jcol` belongs to group `ngrp(jcol)`.
    real(wp),dimension(n) :: d        !! an input array of length `n` which contains the
                                      !! difference parameter vector for the estimate of
                                      !! the jacobian matrix columns in group `numgrp`.
    real(wp),dimension(m) :: Fjacd    !! an input array of length `m` which contains
                                      !! an approximation to the difference vector `jac*d`,
                                      !! where `jac` denotes the jacobian matrix.
    real(wp),dimension(*) :: Fjac     !! an output array of length `nnz`, where `nnz` is the
                                      !! number of its non-zero elements. at each call of `fdjs`,
                                      !! `fjac` is updated to include the non-zero elements of the
                                      !! jacobian matrix for those columns in group `numgrp`. `fjac`
                                      !! should not be altered between successive calls to `fdjs`.
    logical,intent(in)    :: Col      !! a logical input variable. if `col` is set true, then the
                                      !! jacobian approximations are stored into a column-oriented
                                      !! pattern. if `col` is set false, then the jacobian
                                      !! approximations are stored into a row-oriented pattern.

    integer :: ip , irow , jcol , jp

    ! compute estimates of jacobian matrix columns in group
    ! numgrp. the array fjacd must contain an approximation
    ! to jac*d, where jac denotes the jacobian matrix and d
    ! is a difference parameter vector with d(jcol) non-zero
    ! if and only if jcol is a column in group numgrp.

    if ( Col ) then  ! column orientation.

        do jcol = 1 , n
            if ( Ngrp(jcol)==Numgrp ) then
                do jp = Npntr(jcol) , Npntr(jcol+1) - 1
                    irow = Ind(jp)
                    Fjac(jp) = Fjacd(irow)/d(jcol)
                enddo
            endif
        enddo

    else  ! row orientation.

        do irow = 1 , m
            do ip = Npntr(irow) , Npntr(irow+1) - 1
                jcol = Ind(ip)
                if ( Ngrp(jcol)==Numgrp ) then
                    Fjac(ip) = Fjacd(irow)/d(jcol)
                    exit
                endif
            enddo
        enddo

    endif

    end subroutine fdjs
!*******************************************************************************

!*******************************************************************************
    end module dsm_module
!*******************************************************************************
