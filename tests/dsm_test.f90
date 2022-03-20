    program dsm_test

    !! This is a test program for subroutines [[dsm]] and [[fdjs]].
    !! the test data represents a neutron kinetics problem.
    !!
    !!### Reference
    !!  * Argonne National Laboratory. MINPACK Project. July 1983.
    !!    Thomas F. Coleman, Burton S. Garbow, Jorge J. More
    !!  * Thomas F. Coleman, Burton S. Garbow, Jorge J. More, "Algorithm 618:
    !!    FORTRAN subroutines for estimating sparse Jacobian Matrices",
    !!    ACM Transactions on Mathematical Software (TOMS),
    !!    Volume 10 Issue 3, Sept. 1984, Pages 346-347

    use dsm_module
    use iso_fortran_env, only: output_unit
    use numdiff_kinds_module, only: wp

    implicit none

      integer,parameter :: nwrite = output_unit  !! unit for printing

      integer i , info , ip , j , jp , l , m , maxgrp , maxrow , &
              mingrp , minrow , n , nnz , numgrp
      integer indcol(6000) , indrow(6000) , ipntr(1201) , jpntr(1201) , &
              ngrp(1200)
      logical col
      real(wp) :: dnsm , errij , errmax , fjact , fjactr , sum
      real(wp) :: d(1200) , fjac(6000) , fjacd(1200) , fvec(1200) , x(1200) ,  &
                  xd(1200)

      col = .true.
!
!     TEST FOR DSM AND FDJS.
!
      write (nwrite,99001)
!
!     FORMAT STATEMENTS.
!
99001 format (//' TESTS FOR DSM AND FDJS - NEUTRON KINETICS PROBLEM'//  &
             &' STATISTICS GENERATED '//'       N - NUMBER OF COLUMNS '/&
             &'     NNZ - NUMBER OF NON-ZERO ELEMENTS'/                 &
             &'    DNSM - MATRIX DENSITY (PERCENTAGE)'/                 &
             &'  MINROW - MINIMUM NUMBER OF NON-ZEROS IN ANY ROW'/      &
             &'  MAXROW - MAXIMUM NUMBER OF NON-ZEROS IN ANY ROW'/      &
             &'  MINGRP - LOWER BOUND ON NUMBER OF GROUPS'/             &
             &'  MAXGRP - NUMBER OF GROUPS DETERMINED BY DSM'//)
      do n = 300 , 1200 , 300
         write (nwrite,99002)
99002    format (//3x,'N',6x,'NNZ',5x,'DNSM',5x,'MINROW',4x,'MAXROW',4x,&
                &'MINGRP',4x,'MAXGRP'//)
!
!        DEFINITION OF SPARSITY PATTERN.
!
         m = n
         l = n/3
         nnz = 0
         do j = 1 , n
            nnz = nnz + 1
            indrow(nnz) = j
            indcol(nnz) = j
            if ( mod(j,l)/=0 ) then
               nnz = nnz + 1
               indrow(nnz) = j + 1
               indcol(nnz) = j
            endif
            if ( j<=2*l ) then
               nnz = nnz + 1
               indrow(nnz) = j + l
               indcol(nnz) = j
               if ( mod(j,l)/=1 ) then
                  nnz = nnz + 1
                  indrow(nnz) = j - 1
                  indcol(nnz) = j
               endif
            endif
            nnz = nnz + 1
            if ( j>l ) then
               indrow(nnz) = j - l
            else
               indrow(nnz) = j + 2*l
            endif
            indcol(nnz) = j
         enddo
!
!        CALL DSM.
!
         call dsm(m,n,nnz,indrow,indcol,ngrp,maxgrp,mingrp,info,ipntr,jpntr)
         if ( info<=0 ) write (nwrite,99003) info
99003    format (//' *** MISTAKE IN INPUT DATA, INFO IS ***',i6)
!
!        STATISTICS FOR THE MATRIX.
!
         maxrow = 0
         minrow = n
         do i = 1 , m
            maxrow = max(maxrow,ipntr(i+1)-ipntr(i))
            minrow = min(minrow,ipntr(i+1)-ipntr(i))
         enddo
         dnsm = real(100*nnz,wp)/real(m*n,wp)
         write (nwrite,99004) n , nnz , dnsm , minrow , maxrow , &
                            & mingrp , maxgrp
99004    format (2(i5,3x),f6.2,4x,4(i5,5x))
!
!        TEST FOR FDJS.
!
         do j = 1 , n
            x(j) = real(j,wp)/real(n,wp)
         enddo
         call fcn(n,x,indcol,ipntr,fvec)
!
!        APPROXIMATE THE JACOBIAN MATRIX.
!
         do numgrp = 1 , maxgrp
            do j = 1 , n
               d(j) = 0.0_wp
               if ( ngrp(j)==numgrp ) d(j) = 1.0e-6_wp !d(j) = 0.001_wp
               xd(j) = x(j) + d(j)
            enddo
            call fcn(n,xd,indcol,ipntr,fjacd)
            do i = 1 , m
               fjacd(i) = fjacd(i) - fvec(i)
            enddo
            if ( col ) then
               call fdjs(m,n,col,indrow,jpntr,ngrp,numgrp,d,fjacd,fjac)
            else
               call fdjs(m,n,col,indcol,ipntr,ngrp,numgrp,d,fjacd,fjac)
            endif
         enddo
!
!        TEST THE APPROXIMATION TO THE JACOBIAN.
!
         errmax = 0.0_wp
         if ( col ) then
!
!           TEST FOR THE COLUMN-ORIENTED DEFINITION OF
!           THE SPARSITY PATTERN.
!
            do j = 1 , n
               do jp = jpntr(j) , jpntr(j+1) - 1
                  i = indrow(jp)
                  sum = 0.0_wp
                  do ip = ipntr(i) , ipntr(i+1) - 1
                     sum = sum + x(indcol(ip))
                  enddo
                  sum = sum + x(i)
                  fjact = 1.0_wp + 2.0_wp*sum
                  if ( i==j ) fjact = 2.0_wp*fjact
                  errij = fjac(jp) - fjact
                  if ( fjact/=0.0_wp ) errij = errij/fjact
                  errmax = max(errmax,abs(errij))
               enddo
            enddo
         else
!
!           TEST FOR THE ROW-ORIENTED DEFINITION OF
!           THE SPARSITY PATTERN.
!
            do i = 1 , m
               sum = 0.0_wp
               do ip = ipntr(i) , ipntr(i+1) - 1
                  sum = sum + x(indcol(ip))
               enddo
               sum = sum + x(i)
               fjactr = 1.0_wp + 2.0_wp*sum
               do ip = ipntr(i) , ipntr(i+1) - 1
                  j = indcol(ip)
                  fjact = fjactr
                  if ( i==j ) fjact = 2.0_wp*fjact
                  errij = fjac(ip) - fjact
                  if ( fjact/=0.0_wp ) errij = errij/fjact
                  errmax = max(errmax,abs(errij))
               enddo
            enddo
         endif
         write (nwrite,99005) errmax
99005    format (//' LARGEST RELATIVE ERROR OF APPROXIMATION IS',e10.2)
         col = .not.col
      enddo
      stop


  contains

      subroutine fcn(n,x,Indcol,Ipntr,Fvec)

      !! Function subroutine for testing [[fdjs]].

      implicit none

      integer n
      integer Indcol(*) , Ipntr(n+1)
      real(wp) ::x(n) , Fvec(n)

      integer i , ip
      real(wp) :: sum

      do i = 1 , n
         sum = 0.0_wp
         do ip = Ipntr(i) , Ipntr(i+1) - 1
            sum = sum + x(Indcol(ip))
         enddo
         sum = sum + x(i)
         Fvec(i) = sum*(1.0_wp+sum) + 1.0_wp
      enddo

      end subroutine fcn

    end program dsm_test
