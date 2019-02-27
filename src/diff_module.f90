!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!  Numerical differentiation of a 1D function `f(x)` using Neville's process.
!
!## Authors
!   * J. Oliver, "An algorithm for numerical differentiation of a function
!     of one real variable", Journal of Computational and Applied Mathematics
!     6 (2) (1980) 145–160. [Algol 60 source in original paper]
!   * David Kahaner, Fortran 77 code from
!     [NIST](ftp://math.nist.gov/pub/repository/diff/src/DIFF)
!   * Jacob Williams : 2/17/2013 : Converted to modern Fortran.
!     Some refactoring, addition of test cases.

    module diff_module

    use numdiff_kinds_module

    implicit none

    private

    type,public :: diff_func
        !! class to define the function for [[diff]]
        private
        logical :: stop = .false.
        procedure(func),pointer :: f => null()
        contains
        private
        procedure :: faccur
        procedure,public :: set_function
        procedure,public :: compute_derivative => diff
        procedure,public :: terminate
    end type diff_func

    abstract interface
        function func(me,x) result(fx)
            !! interface to function for [[diff]]
            import :: diff_func,wp
            implicit none
            class(diff_func),intent(inout) :: me
            real(wp),intent(in) :: x
            real(wp) :: fx
        end function func
    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date:  12/27/2015
!
!  Set the function in a [[diff_func]].
!  Must be called before [[diff]].

    subroutine set_function(me,f)

    implicit none

    class(diff_func),intent(inout) :: me
    procedure(func) :: f

    me%f => f

    end subroutine set_function
!*****************************************************************************************

!*****************************************************************************************
!>
!  Can be called by the user in the function to terminate the computation.
!  This will set `ifail=-1`.

    subroutine terminate(me)

    implicit none

    class(diff_func),intent(inout) :: me

    me%stop = .true.

    end subroutine terminate
!*****************************************************************************************

!*****************************************************************************************
!>
!  the procedure `diff` calculates the first, second or
!  third order derivative of a function by using neville's process to
!  extrapolate from a sequence of simple polynomial approximations based on
!  interpolating points distributed symmetrically about `x0` (or lying only on
!  one side of `x0` should this be necessary).  if the specified tolerance is
!  non-zero then the procedure attempts to satisfy this absolute or relative
!  accuracy requirement, while if it is unsuccessful or if the tolerance is
!  set to zero then the result having the minimum achievable estimated error
!  is returned instead.

    subroutine diff(me,iord,x0,xmin,xmax,eps,accr,deriv,error,ifail)

    implicit none

    class(diff_func),intent(inout) :: me
    integer,intent(in)    :: iord   !! 1, 2 or 3 specifies that the first, second or third order
                                    !! derivative,respectively, is required.
    real(wp), intent(in)  :: x0     !! the point at which the derivative of the function is to be calculated.
    real(wp), intent(in)  :: xmin   !! `xmin`, `xmax` restrict the interpolating points to lie in [`xmin`, `xmax`], which
                                    !! should be the largest interval including `x0` in which the function is
                                    !! calculable and continuous.
    real(wp), intent(in)  :: xmax   !! `xmin`, `xmax` restrict the interpolating points to lie in [`xmin`, `xmax`], which
                                    !! should be the largest interval including `x0` in which the function is
                                    !! calculable and continuous.
    real(wp), intent(in)  :: eps    !! denotes the tolerance, either absolute or relative. `eps=0` specifies that
                                    !! the error is to be minimised, while `eps>0` or `eps<0` specifies that the
                                    !! absolute or relative error, respectively, must not exceed `abs(eps)` if
                                    !! possible.  the accuracy requirement should not be made stricter than
                                    !! necessary, since the amount of computation tends to increase as
                                    !! the magnitude of `eps` decreases, and is particularly high when `eps=0`.
    real(wp), intent(in)  :: accr   !! denotes that the absolute (`accr>0`) or relative (`accr<0`) errors in the
                                    !! computed values of the function are most unlikely to exceed `abs(accr)`, which
                                    !! should be as small as possible.  if the user cannot estimate `accr` with
                                    !! complete confidence, then it should be set to zero.
    real(wp), intent(out) :: deriv  !! the calculated value of the derivative
    real(wp), intent(out) :: error  !! an estimated upper bound on the magnitude of the absolute error in
                                    !! the calculated result.  it should always be examined, since in extreme case
                                    !! may indicate that there are no correct significant digits in the value
                                    !! returned for derivative.
    integer, intent(out)  :: ifail  !! will have one of the following values on exit:
                                    !!  *0* the procedure was successful.
                                    !!  *1* the estimated error in the result exceeds the (non-zero) requested
                                    !!      error, but the most accurate result possible has been returned.
                                    !!  *2* input data incorrect (derivative and error will be undefined).
                                    !!  *3* the interval [`xmin`, `xmax`] is too small (derivative and error will be
                                    !!      undefined).
                                    !!  *-1* stopped by the user.

    real(wp) :: acc,beta,beta4,h,h0,h1,h2, &
        newh1,newh2,heval,hprev,baseh,hacc1,hacc2,nhacc1, &
        nhacc2,minh,maxh,maxh1,maxh2,tderiv,f0,twof0,f1,f2,f3,f4,fmax, &
        maxfun,pmaxf,df1,deltaf,pdelta,z,zpower,c0f0,c1,c2,c3,dnew,dprev, &
        re,te,newerr,temerr,newacc,pacc1,pacc2,facc1,facc2,acc0, &
        acc1,acc2,relacc,twoinf,twosup,s, &
        d(10),denom(10),e(10),minerr(10),maxf(0:10),save(0:13),storef(-45:45),factor
    integer :: i,j,k,n,nmax,method,signh,fcount,init
    logical :: ignore(10),contin,saved
    real(wp) :: dummy1,dummy2

    integer,parameter :: eta = digits(1.0_wp) - 1       !! minimum number of significant binary digits (apart from the
                                                        !! sign digit) used to represent the mantissa of real(wp) numbers. it should
                                                        !! be decreased by one if the computer truncates rather than rounds.
    integer,parameter :: inf = -minexponent(1.0_wp) - 2 !! the largest possible positive integers subject to
                                                        !! 2**(-inf) and -2**(-inf) being representable real(wp) numbers.
    integer,parameter :: sup = maxexponent(1.0_wp) - 1  !! the largest possible positive integers subject to
                                                        !! 2**sup and -2**sup being representable real(wp) numbers.

    real(wp),parameter :: sqrt2 = sqrt(2.0_wp)  !! \( \sqrt(2) \)
    real(wp),parameter :: sqrt3 = sqrt(3.0_wp)  !! \( \sqrt(3) \)

    me%stop = .false.

    ! execution commences with examination of input parameters
    if (iord<1 .or. iord>3 .or. xmax<=xmin .or. &
        x0>xmax .or. x0<xmin .or. .not. associated(me%f)) then

      ifail = 2

    else

        acc = accr
        twoinf = 2.0_wp**(-inf)
        twosup = 2.0_wp**sup
        factor = 2.0_wp**(real((inf+sup),wp)/30.0_wp)
        if (factor < 256.0_wp)factor=256.0_wp
        maxh1 = xmax - x0
        signh = 1
        if (x0-xmin <= maxh1) then
          maxh2 = x0 - xmin
        else
          maxh2 = maxh1
          maxh1 = x0 - xmin
          signh = -1
        end if
        relacc = 2.0_wp**(1.0_wp-eta)
        maxh1 = (1.0_wp-relacc)*maxh1
        maxh2 = (1.0_wp-relacc)*maxh2
        s=128.0_wp*twoinf
        if (abs(x0) > 128.0_wp*twoinf*2.0_wp**eta) s = abs(x0)*2.0_wp**(-eta)
        if (maxh1 < s) then
          ! interval too small
          ifail =3
          return
        end if
        if (acc < 0.0_wp) then
          if (-acc > relacc) relacc = -acc
          acc = 0.0_wp
        end if

        ! determine the smallest spacing at which the calculated
        ! function values are unequal near x0.

        f0 = me%f(x0)
        if (me%stop) then
          ifail = -1
          return
        end if
        twof0 = f0 + f0
        if (abs(x0) > twoinf*2.0_wp**eta) then
          h = abs(x0)*2.0_wp**(-eta)
          z = 2.0_wp
        else
          h = twoinf
          z = 64.0_wp
        end if
        df1 = me%f(x0+signh*h) - f0
        if (me%stop) then
          ifail = -1
          return
        end if
        do
            if (df1 /= 0.0_wp .or. z*h > maxh1) exit
            h = z*h
            df1 = me%f(x0+signh*h) - f0
            if (me%stop) then
              ifail = -1
              return
            end if
            if (z /= 2.0_wp) then
              if (df1 /= 0.0_wp) then
                h = h/z
                z = 2.0_wp
                df1 = 0.0_wp
              else
                if (z*h > maxh1) z = 2.0_wp
              end if
            end if
        end do

        if (df1 == 0.0_wp) then
        ! constant function
          deriv = 0.0_wp
          error = 0.0_wp
          ifail = 0
          return
        end if
        if (h > maxh1/128.0_wp) then
        ! minimum h too large
          ifail = 3
          return
        end if

        h = 8.0_wp*h
        h1 = signh*h
        h0 = h1
        h2 = -h1
        minh = 2.0_wp**(-min(inf,sup)/iord)
        if (minh < h) minh = h
        select case (iord)
        case(1)
            s = 8.0_wp
        case(2)
            s = 9.0_wp*sqrt3
        case(3)
            s = 27.0_wp
        end select
        if (minh > maxh1/s) then
          ifail = 3
          return
        end if
        if (minh > maxh2/s .or. maxh2 < 128.0_wp*twoinf) then
          method = 1
        else
          method = 2
        end if

        ! method 1 uses 1-sided formulae, and method 2 symmetric.
        ! now estimate accuracy of calculated function values.

        if (method /= 2 .or. iord == 2) then
          if (x0 /= 0.0_wp) then
            dummy1 = 0.0_wp
            dummy2 = -h1
            call me%faccur(dummy1,dummy2,acc0,x0,twoinf,f0,f1,ifail)
            if (ifail==-1) return
          else
            acc0 = 0.0_wp
          end if
        end if

        if (abs(h1) > twosup/128.0_wp) then
          hacc1 = twosup
        else
          hacc1 = 128.0_wp*h1
        end if

        if (abs(hacc1)/4.0_wp < minh) then
          hacc1 = 4.0_wp*signh*minh
        else if (abs(hacc1) > maxh1) then
          hacc1 = signh*maxh1
        end if
        f1 = me%f(x0+hacc1)
        if (me%stop) then
          ifail = -1
          return
        end if
        call me%faccur(hacc1,h1,acc1,x0,twoinf,f0,f1,ifail)
        if (ifail==-1) return
        if (method == 2) then
          hacc2 = -hacc1
          if (abs(hacc2) > maxh2) hacc2 = -signh * maxh2
          f1 = me%f(x0 + hacc2)
          if (me%stop) then
            ifail = -1
            return
          end if
          call me%faccur(hacc2,h2,acc2,x0,twoinf,f0,f1,ifail)
          if (ifail==-1) return
        end if
        nmax = 8
        if (eta > 36) nmax = 10
        n = -1
        fcount = 0
        deriv = 0.0_wp
        error = twosup
        init = 3
        contin = .true.

        do

            n = n+1
            if (.not. contin) exit

            if (init == 3) then
            ! calculate coefficients for differentiation
            ! formulae and neville extrapolation algorithm
              if (iord == 1) then
                beta=2.0_wp
              else if (method == 2) then
                beta = sqrt2
              else
                beta = sqrt3
              end if
              beta4 = beta**4
              z = beta
              if (method == 2) z = z**2
              zpower = 1.0_wp
              do k = 1,nmax
                zpower = z*zpower
                denom(k) = zpower-1
              end do
              if (method == 2 .and. iord == 1) then
                e(1) = 5.0_wp
                e(2) = 6.3_wp
                do i = 3,nmax
                  e(i) = 6.81_wp
                end do
              else if ((method /= 2 .and. iord == 1) .or. &
                       (method == 2 .and. iord == 2)) then
                e(1) = 10.0_wp
                e(2) = 16.0_wp
                e(3) = 20.36_wp
                e(4) = 23.0_wp
                e(5) = 24.46_wp
                do i = 6,nmax
                  e(i) = 26.0_wp
                end do
                if (method == 2 .and. iord == 2) then
                  do i = 1,nmax
                    e(i)=2.0_wp*e(i)
                  end do
                end if
              else if (method /= 2 .and. iord == 2) then
                e(1) = 17.78_wp
                e(2) = 30.06_wp
                e(3) = 39.66_wp
                e(4) = 46.16_wp
                e(5) = 50.26_wp
                do i = 6,nmax
                  e(i) = 55.0_wp
                end do
              else if (method == 2 .and. iord == 3) then
                e(1) = 25.97_wp
                e(2) = 41.22_wp
                e(3) = 50.95_wp
                e(4) = 56.4_wp
                e(5) = 59.3_wp
                do i = 6,nmax
                  e(i) = 62.0_wp
                end do
              else
                e(1) = 24.5_wp
                e(2) = 40.4_wp
                e(3) = 52.78_wp
                e(4) = 61.2_wp
                e(5) = 66.55_wp
                do i = 6,nmax
                  e(i) = 73.0_wp
                end do
                c0f0 = -twof0/(3.0_wp*beta)
                c1 = 3.0_wp/(3.0_wp*beta-1.0_wp)
                c2 = -1.0_wp/(3.0_wp*(beta-1.0_wp))
                c3 = 1.0_wp/(3.0_wp*beta*(5.0_wp-2.0_wp*beta))
              end if
            end if

            if (init >= 2) then
            ! initialization of steplengths, accuracy and other parameters

              heval = signh*minh
              h = heval
              baseh = heval
              maxh = maxh2
              if (method == 1)maxh = maxh1
              do k = 1,nmax
                minerr(k) = twosup
                ignore(k) = .false.
              end do
              if (method == 1) newacc = acc1
              if (method == -1) newacc = acc2
              if (method == 2) newacc = (acc1+acc2)/2.0_wp
              if (newacc < acc) newacc = acc
              if ((method /= 2 .or. iord == 2) .and. newacc < acc0) newacc = acc0
              if (method /= -1) then
                facc1 = acc1
                nhacc1 = hacc1
                newh1 = h1
              end if
              if (method /= 1) then
                facc2 = acc2
                nhacc2 = hacc2
                newh2 = h2
              else
                facc2 = 0.0_wp
                nhacc2 = 0.0_wp
              end if
              init = 1
              j = 0
              saved = .false.
            end if

            ! calculate new or initial function values

            if (init == 1 .and. (n == 0 .or. iord == 1) .and. &
                  .not.(method == 2 .and. fcount >= 45)) then
              if (method == 2) then
                fcount = fcount + 1
                f1 = me%f(x0+heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
                storef(fcount) = f1
                f2 = me%f(x0-heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
                storef(-fcount) = f2
              else
                j = j+1
                if (j <= fcount) then
                  f1 = storef(j*method)
                else
                  f1 = me%f(x0+heval)
                  if (me%stop) then
                    ifail = -1
                    return
                  end if
                end if
              end if
            else
              f1 = me%f(x0+heval)
              if (me%stop) then
                ifail = -1
                return
              end if
              if (method == 2) then
                f2 = me%f(x0-heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
              end if
            end if
            if (n == 0) then
              if (method == 2 .and. iord == 3) then
                pdelta = f1-f2
                pmaxf = (abs(f1)+abs(f2))/2.0_wp
                heval = beta*heval
                f1 = me%f(x0+heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
                f2 = me%f(x0-heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
                deltaf = f1-f2
                maxfun = (abs(f1)+abs(f2))/2.0_wp
                heval = beta*heval
                f1 = me%f(x0+heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
                f2 = me%f(x0-heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
              else if (method /= 2 .and. iord >= 2) then
                if (iord == 2) then
                  f3 = f1
                else
                  f4 = f1
                  heval = beta*heval
                  f3 = me%f(x0+heval)
                  if (me%stop) then
                    ifail = -1
                    return
                  end if
                end if
                heval = beta*heval
                f2 = me%f(x0+heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
                heval = beta*heval
                f1 = me%f(x0+heval)
                if (me%stop) then
                  ifail = -1
                  return
                end if
              end if
            end if

            ! evaluate a new approximation dnew to the derivative

            if (n > nmax) then
              n = nmax
              do i = 1,n
                maxf(i-1) = maxf(i)
              end do
            end if
            if (method == 2) then
              maxf(n) = (abs(f1)+abs(f2))/2.0_wp
              if (iord == 1) then
                dnew = (f1-f2)/2.0_wp
              else if (iord == 2) then
                dnew = f1+f2-twof0
              else
                dnew = -pdelta
                pdelta = deltaf
                deltaf = f1-f2
                dnew = dnew + 0.5_wp*deltaf
                if (maxf(n) < pmaxf) maxf(n) = pmaxf
                pmaxf = maxfun
                maxfun = (abs(f1)+abs(f2))/2.0_wp
              end if
            else
              maxf(n) = abs(f1)
              if (iord == 1) then
                dnew = f1-f0
              else if (iord == 2) then
                dnew = (twof0-3.0_wp*f3+f1)/3.0_wp
                if (maxf(n) < abs(f3)) maxf(n) = abs(f3)
                f3 = f2
                f2 = f1
              else
                dnew = c3*f1+c2*f2+c1*f4+c0f0
                if (maxf(n) < abs(f2)) maxf(n) = abs(f2)
                if (maxf(n) < abs(f4)) maxf(n) = abs(f4)
                f4 = f3
                f3 = f2
                f2 = f1
              end if
            end if
            if (abs(h) > 1) then
              dnew = dnew/h**iord
            else
              if (128.0_wp*abs(dnew) > twosup*abs(h)**iord) then
                dnew = twosup/128.0_wp
              else
                dnew = dnew/h**iord
              end if
            end if

            if (init == 0) then
            ! update estimated accuracy of function values
              newacc = acc
              if ((method /= 2 .or. iord == 2) .and. newacc < acc0) newacc = acc0
              if (method /= -1 .and. abs(nhacc1) <= 1.125_wp*abs(heval)/beta4) then
                nhacc1 = heval
                pacc1 = facc1
                call me%faccur(nhacc1,newh1,facc1,x0,twoinf,f0,f1,ifail)
                if (ifail==-1) return
                if (facc1 < pacc1) facc1=(3.0_wp*facc1+pacc1)/4.0_wp
              end if
              if (method /= 1 .and. abs(nhacc2) <= 1.125_wp*abs(heval)/beta4) then
                if (method == 2) then
                  f1 = f2
                  nhacc2 = -heval
                else
                  nhacc2 = heval
                end if
                pacc2 = facc2
                call me%faccur(nhacc2,newh2,facc2,x0,twoinf,f0,f1,ifail)
                if (ifail==-1) return
                if (facc2 < pacc2) facc2 = (3.0_wp*facc2+pacc2)/4.0_wp
              end if
              if (method == 1 .and. newacc < facc1) newacc = facc1
              if (method == -1 .and. newacc < facc2) newacc = facc2
              if (method == 2 .and. newacc < (facc1+facc2)/2.0_wp) &
                      newacc = (facc1+facc2)/2.0_wp
            end if

            ! evaluate successive elements of the current row in the neville
            ! array, estimating and examining the truncation and rounding
            ! errors in each

            contin = n < nmax
            hprev = abs(h)
            fmax = maxf(n)
            if ((method /= 2 .or. iord == 2) .and. fmax < abs(f0)) fmax = abs(f0)

            do k = 1,n
              dprev = d(k)
              d(k) = dnew
              dnew = dprev+(dprev-dnew)/denom(k)
              te = abs(dnew-d(k))
              if (fmax < maxf(n-k)) fmax = maxf(n-k)
              hprev = hprev/beta
              if (newacc >= relacc*fmax) then
                re = newacc*e(k)
              else
                re = relacc*fmax*e(k)
              end if
              if (re /= 0.0_wp) then
                if (hprev > 1) then
                  re = re/hprev**iord
                else if (2.0_wp*re > twosup*hprev**iord) then
                  re = twosup/2.0_wp
                else
                  re = re/hprev**iord
                end if
              end if
              newerr = te+re
              if (te > re) newerr = 1.25_wp*newerr
              if (.not. ignore(k)) then
                if ((init == 0 .or. (k == 2 .and. .not.ignore(1)))  &
                      .and. newerr < error) then
                  deriv = d(k)
                  error = newerr
                end if
                if (init == 1 .and. n == 1) then
                  tderiv = d(1)
                  temerr = newerr
                end if
                if (minerr(k) < twosup/4.0_wp) then
                  s = 4.0_wp*minerr(k)
                else
                  s = twosup
                end if
                if (te > re .or. newerr > s) then
                  ignore(k) = .true.
                else
                  contin = .true.
                end if
                if (newerr < minerr(k)) minerr(k) = newerr
                if (init == 1 .and. n == 2 .and. k == 1 .and. .not.ignore(1)) then
                  if (newerr < temerr) then
                    tderiv = d(1)
                    temerr = newerr
                  end if
                  if (temerr < error) then
                    deriv = tderiv
                    error = temerr
                  end if
                end if
              end if
            end do

            if (n < nmax) d(n+1) = dnew
            if (eps < 0.0_wp) then
              s = abs(eps*deriv)
            else
              s = eps
            end if
            if (error <= s) then
              contin = .false.
            else if (init == 1 .and. (n == 2 .or. ignore(1))) then
              if ((ignore(1) .or. ignore(2)) .and. saved) then
                saved = .false.
                n = 2
                h = beta * save(0)
                heval = beta*save(1)
                maxf(0) = save(2)
                maxf(1) = save(3)
                maxf(2) = save(4)
                d(1) = save(5)
                d(2) = save(6)
                d(3) = save(7)
                minerr(1) = save(8)
                minerr(2) = save(9)
                if (method == 2 .and. iord == 3) then
                  pdelta = save(10)
                  deltaf = save(11)
                  pmaxf = save(12)
                  maxfun = save(13)
                else if (method /= 2 .and. iord >= 2) then
                  f2 = save(10)
                  f3 = save(11)
                  if (iord == 3) f4 = save(12)
                end if
                init = 0
                ignore(1) = .false.
                ignore(2) = .false.
              else if (.not. (ignore(1) .or. ignore(2)) .and. n == 2  &
                    .and. beta4*factor*abs(heval) <= maxh) then
                ! save all current values in case of return to current point
                saved = .true.
                save(0) = h
                save(1) = heval
                save(2) = maxf(0)
                save(3) = maxf(1)
                save(4) = maxf(2)
                save(5) = d(1)
                save(6) = d(2)
                save(7) = d(3)
                save(8) = minerr(1)
                save(9) = minerr (2)
                if (method == 2 .and. iord == 3) then
                  save(10) = pdelta
                  save(11) = deltaf
                  save(12) = pmaxf
                  save(13) = maxfun
                else if (method /= 2 .and. iord >= 2) then
                  save(10) = f2
                  save(11) = f3
                  if (iord == 3) save(12) = f4
                end if
                h = factor*baseh
                heval = h
                baseh = h
                n = -1
              else
                init = 0
                h = beta*h
                heval = beta*heval
              end if
            else if (contin .and. beta*abs(heval) <= maxh) then
              h = beta*h
              heval = beta*heval
            else if (method /= 1) then
              contin = .true.
              if (method == 2) then
                init = 3
                method = -1
                if (iord /= 2) then
                  if (x0 /= 0.0_wp) then
                    dummy1 = 0.0_wp
                    dummy2 = -h0
                    call me%faccur(dummy1,dummy2,acc0,x0,twoinf,f0,f1,ifail)
                    if (ifail==-1) return
                  else
                    acc0 = 0.0_wp
                  end if
                end if
              else
                init = 2
                method = 1
              end if
              n = -1
              signh = -signh
            else
              contin = .false.
            end if

        end do

        if (eps < 0.0_wp) then
          s = abs(eps*deriv)
        else
          s = eps
        end if
        ifail = 0
        if (eps /= 0.0_wp .and. error > s) ifail = 1

    end if

    end subroutine diff
!*****************************************************************************************

!*****************************************************************************************
!>
! This procedure attempts to estimate the level of rounding errors in
! the calculated function values near the point `x0+h0` by fitting a
! least-squares straight-line approximation to the function at the
! six points `x0+h0-j*h1`, (`j = 0,1,3,5,7,9`), and then setting `facc` to
! twice the largest deviation of the function values from this line.
! `hi` is adjusted if necessary so that it is approximately 8 times the
! smallest spacing at which the function values are unequal near `x0+h0`.

    subroutine faccur(me,h0,h1,facc,x0,twoinf,f0,f1,ifail)

    implicit none

    class(diff_func),intent(inout) :: me
    real(wp), intent(inout)  :: h0
    real(wp), intent(inout)  :: h1
    real(wp), intent(out)    :: facc
    real(wp), intent(in)     :: x0
    real(wp), intent(in)     :: twoinf
    real(wp), intent(in)     :: f0
    real(wp), intent(in)     :: f1
    integer, intent(out)     :: ifail  !! 0 if no error, -1 if user termination.

    real(wp) :: a0,a1,f00,f001,f2,deltaf,t0,t1,df(5)
    integer :: j

    ifail = 0
    t0 = 0.0_wp
    t1 = 0.0_wp
    if (h0 /= 0.0_wp) then
      if (x0+h0 /= 0.0_wp) then
        f00 = f1
      else
        h0 = 0.875_wp*h0
        f00 = me%f(x0+h0)
        if (me%stop) then
          ifail = -1
          return
        end if
      end if
      if (abs(h1) >= 32.0_wp*twoinf) h1 = h1/8.0_wp
      if (16.0_wp*abs(h1) > abs(h0)) h1 = sign(h1,1.0_wp)*abs(h0)/16.0_wp
      f001 = me%f(x0+h0-h1)
      if (me%stop) then
        ifail = -1
        return
      end if
      if (f001 == f00) then
        if (256.0_wp*abs(h1) <= abs(h0)) then
          h1 = 2.0_wp*h1
          do
              f001 = me%f(x0+h0-h1)
              if (me%stop) then
                ifail = -1
                return
              end if
              if (f001 /= f00 .or. 256.0_wp*abs(h1) > abs(h0)) exit
              h1 = 2.0_wp*h1
          end do
          h1 = 8.0_wp*h1
        else
          h1 = sign(h1,1.0_wp)*abs(h0)/16.0_wp
        end if
      else
        if (256.0_wp*twoinf <= abs(h0)) then
          do
              f001 = me%f(x0+h0-h1/2.0_wp)
              if (me%stop) then
                ifail = -1
                return
              end if
              if (f001 == f00 .or. abs(h1) < 4.0_wp*twoinf) exit
              h1 = h1/2.0_wp
          end do
          h1 = 8.0_wp*h1
          if (16.0_wp*abs(h1) > abs(h0)) h1 = sign(h1,1.0_wp)*abs(h0)/16.0_wp
        else
          h1 = sign(h1,1.0_wp)*abs(h0)/16.0_wp
        end if
      end if
    else
      f00 = f0
    end if

    do j = 1,5
      f2 = me%f(x0+h0-real(2*j-1,wp)*h1)
      if (me%stop) then
        ifail = -1
        return
      end if
      df(j) = f2 - f00
      t0 = t0+df(j)
      t1 = t1+real(2*j-1,wp)*df(j)
    end do
    a0 = (33.0_wp*t0-5.0_wp*t1)/73.0_wp
    a1 = (-5.0_wp*t0+1.2_wp*t1)/73.0_wp
    facc = abs(a0)
    do j = 1,5
      deltaf = abs(df(j)-(a0+real(2*j-1,wp)*a1))
      if (facc < deltaf) facc = deltaf
    end do
    facc = 2.0_wp*facc

    end subroutine faccur
!*****************************************************************************************

!*****************************************************************************************
    end module diff_module
!*****************************************************************************************
