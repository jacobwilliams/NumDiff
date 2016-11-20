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

    module dsm_module

    use kinds_module

    implicit none

    private

    type,public :: dsm_partition
    end type dsm_partition

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
!        `indrow(k),indcol(k), k = 1,2,...,npairs`.
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

!     check the input data.

      Info = 0
      if ( m<1 .or. n<1 .or. Npairs<1 ) return
      do k = 1 , Npairs
         Info = -k
         if ( Indrow(k)<1 .or. Indrow(k)>m .or. Indcol(k)<1 .or. Indcol(k)>n ) return
      enddo
      Info = 1

!     sort the data structure by columns.

      call srtdat(n,Npairs,Indrow,Indcol,Jpntr,Iwa)

!     compress the data and determine the number of
!     non-zero elements of a.

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

!     extend the data structure to rows.

      call setr(m,n,Indrow,Jpntr,Indcol,Ipntr,Iwa)

!     determine a lower bound for the number of groups.

      Mingrp = 0
      do i = 1 , m
         Mingrp = max(Mingrp,Ipntr(i+1)-Ipntr(i))
      enddo

!     determine the degree sequence for the intersection
!     graph of the columns of a.

      call degr(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(5*n+1),Iwa(n+1))

!     color the intersection graph of the columns of a
!     with the smallest-last (sl) ordering.

      call slo(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(5*n+1),Iwa(4*n+1),maxclq,&
               Iwa(1),Iwa(n+1),Iwa(2*n+1),Iwa(3*n+1))
      call seq(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(4*n+1),Ngrp,Maxgrp,      &
               Iwa(n+1))
      Mingrp = max(Mingrp,maxclq)

!     exit if the smallest-last ordering is optimal.

      if ( Maxgrp==Mingrp ) return

!     color the intersection graph of the columns of a
!     with the incidence-degree (id) ordering.

      call ido(m,n,Indrow,Jpntr,Indcol,Ipntr,Iwa(5*n+1),Iwa(4*n+1),     &
               maxclq,Iwa(1),Iwa(n+1),Iwa(2*n+1),Iwa(3*n+1))
      call seq(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(4*n+1),Iwa(1),numgrp,    &
               Iwa(n+1))
      Mingrp = max(Mingrp,maxclq)

!     retain the better of the two orderings so far.

      if ( numgrp<Maxgrp ) then
         Maxgrp = numgrp
         do j = 1 , n
            Ngrp(j) = Iwa(j)
         enddo

!        exit if the incidence-degree ordering is optimal.

         if ( Maxgrp==Mingrp ) return
      endif

!     color the intersection graph of the columns of a
!     with the largest-first (lf) ordering.

      call numsrt(n,n-1,Iwa(5*n+1),-1,Iwa(4*n+1),Iwa(2*n+1),Iwa(n+1))
      call seq(n,Indrow,Jpntr,Indcol,Ipntr,Iwa(4*n+1),Iwa(1),numgrp,    &
               Iwa(n+1))

!     retain the best of the three orderings and exit.

      if ( numgrp<Maxgrp ) then
         Maxgrp = numgrp
         do j = 1 , n
            Ngrp(j) = Iwa(j)
         enddo
      endif

    end subroutine dsm
!*******************************************************************************

      subroutine degr(n,Indrow,Jpntr,Indcol,Ipntr,Ndeg,Iwa)

      implicit none

      integer n
      integer Indrow(*) , Jpntr(n+1) , Indcol(*) , Ipntr(*) , Ndeg(n) , &
              Iwa(n)

!     SUBROUTINE DEGR
!
!     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A,
!     THIS SUBROUTINE DETERMINES THE DEGREE SEQUENCE FOR
!     THE INTERSECTION GRAPH OF THE COLUMNS OF A.
!
!     IN GRAPH-THEORY TERMINOLOGY, THE INTERSECTION GRAPH OF
!     THE COLUMNS OF A IS THE LOOPLESS GRAPH G WITH VERTICES
!     A(J), J = 1,2,...,N WHERE A(J) IS THE J-TH COLUMN OF A
!     AND WITH EDGE (A(I),A(J)) IF AND ONLY IF COLUMNS I AND J
!     HAVE A NON-ZERO IN THE SAME ROW POSITION.
!
!     NOTE THAT THE VALUE OF M IS NOT NEEDED BY DEGR AND IS
!     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE DEGR(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,IWA)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       NDEG IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH
!         SPECIFIES THE DEGREE SEQUENCE. THE DEGREE OF THE
!         J-TH COLUMN OF A IS NDEG(J).
!
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH N.

      integer ic , ip , ir , jcol , jp
!
!     INITIALIZATION BLOCK.
!
      do jp = 1 , n
         Ndeg(jp) = 0
         Iwa(jp) = 0
      enddo
!
!     COMPUTE THE DEGREE SEQUENCE BY DETERMINING THE CONTRIBUTIONS
!     TO THE DEGREES FROM THE CURRENT(JCOL) COLUMN AND FURTHER
!     COLUMNS WHICH HAVE NOT YET BEEN CONSIDERED.
!
      do jcol = 2 , n
         Iwa(jcol) = n
!
!        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
!        TO NON-ZEROES IN THE MATRIX.
!
         do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
            ir = Indrow(jp)
!
!           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
!           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.
!
            do ip = Ipntr(ir) , Ipntr(ir+1) - 1
               ic = Indcol(ip)
!
!              ARRAY IWA MARKS COLUMNS WHICH HAVE CONTRIBUTED TO
!              THE DEGREE COUNT OF COLUMN JCOL. UPDATE THE DEGREE
!              COUNTS OF THESE COLUMNS AS WELL AS COLUMN JCOL.
!
               if ( Iwa(ic)<jcol ) then
                  Iwa(ic) = jcol
                  Ndeg(ic) = Ndeg(ic) + 1
                  Ndeg(jcol) = Ndeg(jcol) + 1
               endif
            enddo
         enddo
      enddo

      end subroutine degr

      subroutine ido(m,n,Indrow,Jpntr,Indcol,Ipntr,Ndeg,List,Maxclq,    &
                     Iwa1,Iwa2,Iwa3,Iwa4)

      implicit none

      integer m , n , Maxclq
      integer Indrow(*) , Jpntr(n+1) , Indcol(*) , Ipntr(m+1) , Ndeg(n) &
              , List(n) , Iwa1(0:n-1) , Iwa2(n) , Iwa3(n) , Iwa4(n)

!     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS
!     SUBROUTINE DETERMINES AN INCIDENCE-DEGREE ORDERING OF THE
!     COLUMNS OF A.
!
!     THE INCIDENCE-DEGREE ORDERING IS DEFINED FOR THE LOOPLESS
!     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE
!     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF
!     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION.
!
!     THE INCIDENCE-DEGREE ORDERING IS DETERMINED RECURSIVELY BY
!     LETTING LIST(K), K = 1,...,N BE A COLUMN WITH MAXIMAL
!     INCIDENCE TO THE SUBGRAPH SPANNED BY THE ORDERED COLUMNS.
!     AMONG ALL THE COLUMNS OF MAXIMAL INCIDENCE, IDO CHOOSES A
!     COLUMN OF MAXIMAL DEGREE.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE IDO(M,N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST,
!                      MAXCLQ,IWA1,IWA2,IWA3,IWA4)
!
!     WHERE
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF ROWS OF A.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       NDEG IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE DEGREE SEQUENCE. THE DEGREE OF THE J-TH COLUMN
!         OF A IS NDEG(J).
!
!       LIST IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE INCIDENCE-DEGREE ORDERING OF THE COLUMNS OF A. THE J-TH
!         COLUMN IN THIS ORDER IS LIST(J).
!
!       MAXCLQ IS AN INTEGER OUTPUT VARIABLE SET TO THE SIZE
!         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING.
!
!       IWA1,IWA2,IWA3, AND IWA4 ARE INTEGER WORK ARRAYS OF LENGTH N.

      integer ic , ip , ir , jcol , jp , maxinc , maxlst , ncomp ,      &
              numinc , numlst , numord , numwgt
!
!     SORT THE DEGREE SEQUENCE.
!
      call numsrt(n,n-1,Ndeg,-1,Iwa4,Iwa2,Iwa3)
!
!     INITIALIZATION BLOCK.
!
!     CREATE A DOUBLY-LINKED LIST TO ACCESS THE INCIDENCES OF THE
!     COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS.
!
!     EACH UN-ORDERED COLUMN IC IS IN A LIST (THE INCIDENCE LIST)
!     OF COLUMNS WITH THE SAME INCIDENCE.
!
!     IWA1(NUMINC) IS THE FIRST COLUMN IN THE NUMINC LIST
!     UNLESS IWA1(NUMINC) = 0. IN THIS CASE THERE ARE
!     NO COLUMNS IN THE NUMINC LIST.
!
!     IWA2(IC) IS THE COLUMN BEFORE IC IN THE INCIDENCE LIST
!     UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST
!     COLUMN IN THIS INCIDENCE LIST.
!
!     IWA3(IC) IS THE COLUMN AFTER IC IN THE INCIDENCE LIST
!     UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST
!     COLUMN IN THIS INCIDENCE LIST.
!
!     IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE
!     INCIDENCE OF IC TO THE GRAPH INDUCED BY THE ORDERED
!     COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL)
!     IS THE INCIDENCE-DEGREE ORDER OF COLUMN JCOL.
!
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
!
!     DETERMINE THE MAXIMAL SEARCH LENGTH FOR THE LIST
!     OF COLUMNS OF MAXIMAL INCIDENCE.
!
      maxlst = 0
      do ir = 1 , m
         maxlst = maxlst + (Ipntr(ir+1)-Ipntr(ir))**2
      enddo
      maxlst = maxlst/n
      Maxclq = 0
      numord = 1
!
!     BEGINNING OF ITERATION LOOP.
!
!
!        UPDATE THE SIZE OF THE LARGEST CLIQUE
!        FOUND DURING THE ORDERING.
!
 100  if ( maxinc==0 ) ncomp = 0
      ncomp = ncomp + 1
      if ( maxinc+1==ncomp ) Maxclq = max(Maxclq,ncomp)
!
!        CHOOSE A COLUMN JCOL OF MAXIMAL DEGREE AMONG THE
!        COLUMNS OF MAXIMAL INCIDENCE MAXINC.
!
 200  jp = Iwa1(maxinc)
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
!
!        TERMINATION TEST.
!
         if ( numord>n ) then
!
!     INVERT THE ARRAY LIST.
!
            do jcol = 1 , n
               Iwa2(List(jcol)) = jcol
            enddo
            do jp = 1 , n
               List(jp) = Iwa2(jp)
            enddo
         else
!
!        DELETE COLUMN JCOL FROM THE MAXINC LIST.
!
            if ( Iwa2(jcol)==0 ) then
               Iwa1(maxinc) = Iwa3(jcol)
            else
               Iwa3(Iwa2(jcol)) = Iwa3(jcol)
            endif
            if ( Iwa3(jcol)>0 ) Iwa2(Iwa3(jcol)) = Iwa2(jcol)
!
!        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL.
!
            Iwa4(jcol) = n
!
!        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
!        TO NON-ZEROES IN THE MATRIX.
!
            do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
               ir = Indrow(jp)
!
!           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
!           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.
!
               do ip = Ipntr(ir) , Ipntr(ir+1) - 1
                  ic = Indcol(ip)
!
!              ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO
!              COLUMN JCOL.
!
                  if ( Iwa4(ic)<numord ) then
                     Iwa4(ic) = numord
!
!                 UPDATE THE POINTERS TO THE CURRENT INCIDENCE LISTS.
!
                     numinc = List(ic)
                     List(ic) = List(ic) + 1
                     maxinc = max(maxinc,List(ic))
!
!                 DELETE COLUMN IC FROM THE NUMINC LIST.
!
                     if ( Iwa2(ic)==0 ) then
                        Iwa1(numinc) = Iwa3(ic)
                     else
                        Iwa3(Iwa2(ic)) = Iwa3(ic)
                     endif
                     if ( Iwa3(ic)>0 ) Iwa2(Iwa3(ic)) = Iwa2(ic)
!
!                 ADD COLUMN IC TO THE NUMINC+1 LIST.
!
                     Iwa2(ic) = 0
                     Iwa3(ic) = Iwa1(numinc+1)
                     if ( Iwa1(numinc+1)>0 ) Iwa2(Iwa1(numinc+1)) = ic
                     Iwa1(numinc+1) = ic
                  endif
               enddo
            enddo
!
!        END OF ITERATION LOOP.
!
            goto 100
         endif
      else
         maxinc = maxinc - 1
         goto 200
      endif

      end subroutine ido

      subroutine numsrt(n,Nmax,Num,Mode,Index,Last,Next)

      implicit none

      integer n , Nmax , Mode
      integer Num(n) , Index(n) , Last(0:Nmax) , Next(n)

!     GIVEN A SEQUENCE OF INTEGERS, THIS SUBROUTINE GROUPS
!     TOGETHER THOSE INDICES WITH THE SAME SEQUENCE VALUE
!     AND, OPTIONALLY, SORTS THE SEQUENCE INTO EITHER
!     ASCENDING OR DESCENDING ORDER.
!
!     THE SEQUENCE OF INTEGERS IS DEFINED BY THE ARRAY NUM,
!     AND IT IS ASSUMED THAT THE INTEGERS ARE EACH FROM THE SET
!     0,1,...,NMAX. ON OUTPUT THE INDICES K SUCH THAT NUM(K) = L
!     FOR ANY L = 0,1,...,NMAX CAN BE OBTAINED FROM THE ARRAYS
!     LAST AND NEXT AS FOLLOWS.
!
!           K = LAST(L)
!           WHILE (K .NE. 0) K = NEXT(K)
!
!     OPTIONALLY, THE SUBROUTINE PRODUCES AN ARRAY INDEX SO THAT
!     THE SEQUENCE NUM(INDEX(I)), I = 1,2,...,N IS SORTED.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE NUMSRT(N,NMAX,NUM,MODE,INDEX,LAST,NEXT)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE.
!
!       NMAX IS A POSITIVE INTEGER INPUT VARIABLE.
!
!       NUM IS AN INPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         SEQUENCE OF INTEGERS TO BE GROUPED AND SORTED. IT
!         IS ASSUMED THAT THE INTEGERS ARE EACH FROM THE SET
!         0,1,...,NMAX.
!
!       MODE IS AN INTEGER INPUT VARIABLE. THE SEQUENCE NUM IS
!         SORTED IN ASCENDING ORDER IF MODE IS POSITIVE AND IN
!         DESCENDING ORDER IF MODE IS NEGATIVE. IF MODE IS 0,
!         NO SORTING IS DONE.
!
!       INDEX IS AN INTEGER OUTPUT ARRAY OF LENGTH N SET SO
!         THAT THE SEQUENCE
!
!               NUM(INDEX(I)), I = 1,2,...,N
!
!         IS SORTED ACCORDING TO THE SETTING OF MODE. IF MODE
!         IS 0, INDEX IS NOT REFERENCED.
!
!       LAST IS AN INTEGER OUTPUT ARRAY OF LENGTH NMAX + 1. THE
!         INDEX OF NUM FOR THE LAST OCCURRENCE OF L IS LAST(L)
!         FOR ANY L = 0,1,...,NMAX UNLESS LAST(L) = 0. IN
!         THIS CASE L DOES NOT APPEAR IN NUM.
!
!       NEXT IS AN INTEGER OUTPUT ARRAY OF LENGTH N. IF
!         NUM(K) = L, THEN THE INDEX OF NUM FOR THE PREVIOUS
!         OCCURRENCE OF L IS NEXT(K) FOR ANY L = 0,1,...,NMAX
!         UNLESS NEXT(K) = 0. IN THIS CASE THERE IS NO PREVIOUS
!         OCCURRENCE OF L IN NUM.

      integer i , j , jinc , jl , ju , k , l
!
!     DETERMINE THE ARRAYS NEXT AND LAST.
!
      do i = 0 , Nmax
         Last(i) = 0
      enddo
      do k = 1 , n
         l = Num(k)
         Next(k) = Last(l)
         Last(l) = k
      enddo
      if ( Mode==0 ) return
!
!     STORE THE POINTERS TO THE SORTED ARRAY IN INDEX.
!
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
 50      if ( k/=0 ) then
            Index(i) = k
            i = i + 1
            k = Next(k)
            goto 50
         endif
      enddo

      end subroutine numsrt

      subroutine seq(n,Indrow,Jpntr,Indcol,Ipntr,List,Ngrp,Maxgrp,Iwa)

      implicit none

      integer n , Maxgrp
      integer Indrow(*) , Jpntr(n+1) , Indcol(*) , Ipntr(*) , List(n) , &
              Ngrp(n) , Iwa(n)

!     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS
!     SUBROUTINE DETERMINES A CONSISTENT PARTITION OF THE
!     COLUMNS OF A BY A SEQUENTIAL ALGORITHM.
!
!     A CONSISTENT PARTITION IS DEFINED IN TERMS OF THE LOOPLESS
!     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE
!     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF
!     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION.
!
!     A PARTITION OF THE COLUMNS OF A INTO GROUPS IS CONSISTENT
!     IF THE COLUMNS IN ANY GROUP ARE NOT ADJACENT IN THE GRAPH G.
!     IN GRAPH-THEORY TERMINOLOGY, A CONSISTENT PARTITION OF THE
!     COLUMNS OF A CORRESPONDS TO A COLORING OF THE GRAPH G.
!
!     THE SUBROUTINE EXAMINES THE COLUMNS IN THE ORDER SPECIFIED
!     BY THE ARRAY LIST, AND ASSIGNS THE CURRENT COLUMN TO THE
!     GROUP WITH THE SMALLEST POSSIBLE NUMBER.
!
!     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SEQ AND IS
!     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE SEQ(N,INDROW,JPNTR,INDCOL,IPNTR,LIST,NGRP,MAXGRP,
!                      IWA)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       LIST IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE ORDER TO BE USED BY THE SEQUENTIAL ALGORITHM.
!         THE J-TH COLUMN IN THIS ORDER IS LIST(J).
!
!       NGRP IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE PARTITION OF THE COLUMNS OF A. COLUMN JCOL BELONGS
!         TO GROUP NGRP(JCOL).
!
!       MAXGRP IS AN INTEGER OUTPUT VARIABLE WHICH SPECIFIES THE
!         NUMBER OF GROUPS IN THE PARTITION OF THE COLUMNS OF A.
!
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH N.

      integer ic , ip , ir , j , jcol , jp
!
!     INITIALIZATION BLOCK.
!

      Maxgrp = 0
      do jp = 1 , n
         Ngrp(jp) = n
         Iwa(jp) = 0
      enddo
!
!     BEGINNING OF ITERATION LOOP.
!
      do j = 1 , n
         jcol = List(j)
!
!        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL.
!
!        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
!        TO NON-ZEROES IN THE MATRIX.
!
         do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
            ir = Indrow(jp)
!
!           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
!           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.
!
            do ip = Ipntr(ir) , Ipntr(ir+1) - 1
               ic = Indcol(ip)
!
!              ARRAY IWA MARKS THE GROUP NUMBERS OF THE
!              COLUMNS WHICH ARE ADJACENT TO COLUMN JCOL.
!
               Iwa(Ngrp(ic)) = j
            enddo
         enddo
!
!        ASSIGN THE SMALLEST UN-MARKED GROUP NUMBER TO JCOL.
!
         do jp = 1 , Maxgrp
            if ( Iwa(jp)/=j ) goto 50
         enddo
         Maxgrp = Maxgrp + 1
 50      Ngrp(jcol) = jp
      enddo
!
!        END OF ITERATION LOOP.
!

      end subroutine seq

      subroutine setr(m,n,Indrow,Jpntr,Indcol,Ipntr,Iwa)

      implicit none

      integer m , n
      integer Indrow(*) , Jpntr(n+1) , Indcol(*) , Ipntr(m+1) , Iwa(m)

!     GIVEN A COLUMN-ORIENTED DEFINITION OF THE SPARSITY PATTERN
!     OF AN M BY N MATRIX A, THIS SUBROUTINE DETERMINES A
!     ROW-ORIENTED DEFINITION OF THE SPARSITY PATTERN OF A.
!
!     ON INPUT THE COLUMN-ORIENTED DEFINITION IS SPECIFIED BY
!     THE ARRAYS INDROW AND JPNTR. ON OUTPUT THE ROW-ORIENTED
!     DEFINITION IS SPECIFIED BY THE ARRAYS INDCOL AND IPNTR.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE SETR(M,N,INDROW,JPNTR,INDCOL,IPNTR,IWA)
!
!     WHERE
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF ROWS OF A.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       INDCOL IS AN INTEGER OUTPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       IPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!
!         NOTE THAT IPNTR(1) IS SET TO 1 AND THAT IPNTR(M+1)-1 IS
!         THEN THE NUMBER OF NON-ZERO ELEMENTS OF THE MATRIX A.
!
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH M.

      integer ir , jcol , jp
!
!     STORE IN ARRAY IWA THE COUNTS OF NON-ZEROES IN THE ROWS.
!
      do ir = 1 , m
         Iwa(ir) = 0
      enddo
      do jp = 1 , Jpntr(n+1) - 1
         Iwa(Indrow(jp)) = Iwa(Indrow(jp)) + 1
      enddo
!
!     SET POINTERS TO THE START OF THE ROWS IN INDCOL.
!
      Ipntr(1) = 1
      do ir = 1 , m
         Ipntr(ir+1) = Ipntr(ir) + Iwa(ir)
         Iwa(ir) = Ipntr(ir)
      enddo
!
!     FILL INDCOL.
!
      do jcol = 1 , n
         do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
            ir = Indrow(jp)
            Indcol(Iwa(ir)) = jcol
            Iwa(ir) = Iwa(ir) + 1
         enddo
      enddo

      end subroutine setr

      subroutine slo(n,Indrow,Jpntr,Indcol,Ipntr,Ndeg,List,Maxclq,Iwa1, &
                     Iwa2,Iwa3,Iwa4)

      implicit none

      integer n , Maxclq
      integer Indrow(*) , Jpntr(n+1) , Indcol(*) , Ipntr(*) , Ndeg(n) , &
              List(n) , Iwa1(0:n-1) , Iwa2(n) , Iwa3(n) , Iwa4(n)

!     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS
!     SUBROUTINE DETERMINES THE SMALLEST-LAST ORDERING OF THE
!     COLUMNS OF A.
!
!     THE SMALLEST-LAST ORDERING IS DEFINED FOR THE LOOPLESS
!     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE
!     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF
!     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION.
!
!     THE SMALLEST-LAST ORDERING IS DETERMINED RECURSIVELY BY
!     LETTING LIST(K), K = N,...,1 BE A COLUMN WITH LEAST DEGREE
!     IN THE SUBGRAPH SPANNED BY THE UN-ORDERED COLUMNS.
!
!     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SLO AND IS
!     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE SLO(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST,
!                      MAXCLQ,IWA1,IWA2,IWA3,IWA4)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!
!       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!
!       NDEG IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE DEGREE SEQUENCE. THE DEGREE OF THE J-TH COLUMN
!         OF A IS NDEG(J).
!
!       LIST IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE SMALLEST-LAST ORDERING OF THE COLUMNS OF A. THE J-TH
!         COLUMN IN THIS ORDER IS LIST(J).
!
!       MAXCLQ IS AN INTEGER OUTPUT VARIABLE SET TO THE SIZE
!         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING.
!
!       IWA1,IWA2,IWA3, AND IWA4 ARE INTEGER WORK ARRAYS OF LENGTH N.

      integer ic , ip , ir , jcol , jp , mindeg , numdeg , numord
!
!     INITIALIZATION BLOCK.
!
      mindeg = n
      do jp = 1 , n
         Iwa1(jp-1) = 0
         Iwa4(jp) = n
         List(jp) = Ndeg(jp)
         mindeg = min(mindeg,Ndeg(jp))
      enddo
!
!     CREATE A DOUBLY-LINKED LIST TO ACCESS THE DEGREES OF THE
!     COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS.
!
!     EACH UN-ORDERED COLUMN IC IS IN A LIST (THE DEGREE LIST)
!     OF COLUMNS WITH THE SAME DEGREE.
!
!     IWA1(NUMDEG) IS THE FIRST COLUMN IN THE NUMDEG LIST
!     UNLESS IWA1(NUMDEG) = 0. IN THIS CASE THERE ARE
!     NO COLUMNS IN THE NUMDEG LIST.
!
!     IWA2(IC) IS THE COLUMN BEFORE IC IN THE DEGREE LIST
!     UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST
!     COLUMN IN THIS DEGREE LIST.
!
!     IWA3(IC) IS THE COLUMN AFTER IC IN THE DEGREE LIST
!     UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST
!     COLUMN IN THIS DEGREE LIST.
!
!     IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE
!     DEGREE OF IC IN THE GRAPH INDUCED BY THE UN-ORDERED
!     COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL)
!     IS THE SMALLEST-LAST ORDER OF COLUMN JCOL.
!
      do jp = 1 , n
         numdeg = Ndeg(jp)
         Iwa2(jp) = 0
         Iwa3(jp) = Iwa1(numdeg)
         if ( Iwa1(numdeg)>0 ) Iwa2(Iwa1(numdeg)) = jp
         Iwa1(numdeg) = jp
      enddo
      Maxclq = 0
      numord = n
!
!     BEGINNING OF ITERATION LOOP.
!
!
!        MARK THE SIZE OF THE LARGEST CLIQUE
!        FOUND DURING THE ORDERING.
!
 100  if ( mindeg+1==numord .and. Maxclq==0 ) Maxclq = numord
!
!        CHOOSE A COLUMN JCOL OF MINIMAL DEGREE MINDEG.
!
 200  jcol = Iwa1(mindeg)
      if ( jcol>0 ) then
         List(jcol) = numord
         numord = numord - 1
!
!        TERMINATION TEST.
!
         if ( numord==0 ) then
!
!     INVERT THE ARRAY LIST.
!
            do jcol = 1 , n
               Iwa2(List(jcol)) = jcol
            enddo
            do jp = 1 , n
               List(jp) = Iwa2(jp)
            enddo
         else
!
!        DELETE COLUMN JCOL FROM THE MINDEG LIST.
!
            Iwa1(mindeg) = Iwa3(jcol)
            if ( Iwa3(jcol)>0 ) Iwa2(Iwa3(jcol)) = 0
!
!        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL.
!
            Iwa4(jcol) = 0
!
!        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
!        TO NON-ZEROES IN THE MATRIX.
!
            do jp = Jpntr(jcol) , Jpntr(jcol+1) - 1
               ir = Indrow(jp)
!
!           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
!           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.
!
               do ip = Ipntr(ir) , Ipntr(ir+1) - 1
                  ic = Indcol(ip)
!
!              ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO
!              COLUMN JCOL.
!
                  if ( Iwa4(ic)>numord ) then
                     Iwa4(ic) = numord
!
!                 UPDATE THE POINTERS TO THE CURRENT DEGREE LISTS.
!
                     numdeg = List(ic)
                     List(ic) = List(ic) - 1
                     mindeg = min(mindeg,List(ic))
!
!                 DELETE COLUMN IC FROM THE NUMDEG LIST.
!
                     if ( Iwa2(ic)==0 ) then
                        Iwa1(numdeg) = Iwa3(ic)
                     else
                        Iwa3(Iwa2(ic)) = Iwa3(ic)
                     endif
                     if ( Iwa3(ic)>0 ) Iwa2(Iwa3(ic)) = Iwa2(ic)
!
!                 ADD COLUMN IC TO THE NUMDEG-1 LIST.
!
                     Iwa2(ic) = 0
                     Iwa3(ic) = Iwa1(numdeg-1)
                     if ( Iwa1(numdeg-1)>0 ) Iwa2(Iwa1(numdeg-1)) = ic
                     Iwa1(numdeg-1) = ic
                  endif
               enddo
            enddo
!
!        END OF ITERATION LOOP.
!
            goto 100
         endif
      else
         mindeg = mindeg + 1
         goto 200
      endif

      end subroutine slo

      subroutine srtdat(n,Nnz,Indrow,Indcol,Jpntr,Iwa)

      implicit none

      integer n , Nnz
      integer Indrow(Nnz) , Indcol(Nnz) , Jpntr(n+1) , Iwa(n)

!     GIVEN THE NON-ZERO ELEMENTS OF AN M BY N MATRIX A IN
!     ARBITRARY ORDER AS SPECIFIED BY THEIR ROW AND COLUMN
!     INDICES, THIS SUBROUTINE PERMUTES THESE ELEMENTS SO
!     THAT THEIR COLUMN INDICES ARE IN NON-DECREASING ORDER.
!
!     ON INPUT IT IS ASSUMED THAT THE ELEMENTS ARE SPECIFIED IN
!
!           INDROW(K),INDCOL(K), K = 1,...,NNZ.
!
!     ON OUTPUT THE ELEMENTS ARE PERMUTED SO THAT INDCOL IS
!     IN NON-DECREASING ORDER. IN ADDITION, THE ARRAY JPNTR
!     IS SET SO THAT THE ROW INDICES FOR COLUMN J ARE
!
!           INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!
!     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SRTDAT AND IS
!     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE SRTDAT(N,NNZ,INDROW,INDCOL,JPNTR,IWA)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!
!       NNZ IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF NON-ZERO ELEMENTS OF A.
!
!       INDROW IS AN INTEGER ARRAY OF LENGTH NNZ. ON INPUT INDROW
!         MUST CONTAIN THE ROW INDICES OF THE NON-ZERO ELEMENTS OF A.
!         ON OUTPUT INDROW IS PERMUTED SO THAT THE CORRESPONDING
!         COLUMN INDICES OF INDCOL ARE IN NON-DECREASING ORDER.
!
!       INDCOL IS AN INTEGER ARRAY OF LENGTH NNZ. ON INPUT INDCOL
!         MUST CONTAIN THE COLUMN INDICES OF THE NON-ZERO ELEMENTS
!         OF A. ON OUTPUT INDCOL IS PERMUTED SO THAT THESE INDICES
!         ARE IN NON-DECREASING ORDER.
!
!       JPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN THE OUTPUT
!         INDROW. THE ROW INDICES FOR COLUMN J ARE
!
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!
!         NOTE THAT JPNTR(1) IS SET TO 1 AND THAT JPNTR(N+1)-1
!         IS THEN NNZ.
!
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH N.

      integer i , j , k , l
!
!     STORE IN ARRAY IWA THE COUNTS OF NON-ZEROES IN THE COLUMNS.
!
      do j = 1 , n
         Iwa(j) = 0
      enddo
      do k = 1 , Nnz
         Iwa(Indcol(k)) = Iwa(Indcol(k)) + 1
      enddo
!
!     SET POINTERS TO THE START OF THE COLUMNS IN INDROW.
!
      Jpntr(1) = 1
      do j = 1 , n
         Jpntr(j+1) = Jpntr(j) + Iwa(j)
         Iwa(j) = Jpntr(j)
      enddo
      k = 1
!
!     BEGIN IN-PLACE SORT.
!
 100  j = Indcol(k)
      if ( k>=Jpntr(j) ) then
!
!           CURRENT ELEMENT IS IN POSITION. NOW EXAMINE THE
!           NEXT ELEMENT OR THE FIRST UN-SORTED ELEMENT IN
!           THE J-TH GROUP.
!
         k = max(k+1,Iwa(j))
      else
!
!           CURRENT ELEMENT IS NOT IN POSITION. PLACE ELEMENT
!           IN POSITION AND MAKE THE DISPLACED ELEMENT THE
!           CURRENT ELEMENT.
!
         l = Iwa(j)
         Iwa(j) = Iwa(j) + 1
         i = Indrow(k)
         Indrow(k) = Indrow(l)
         Indcol(k) = Indcol(l)
         Indrow(l) = i
         Indcol(l) = j
      endif
      if ( k<=Nnz ) goto 100

      end subroutine srtdat

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
!     d(jcol) = 0.0
!     if (ngrp(jcol) == numgrp) d(jcol) = eta(jcol)
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
