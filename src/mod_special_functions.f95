! Module containing some special functions like erf(x), erfc(x), phi(x) etc...
!
! References:
!
!   [1] 'Computer Approximations', Hart, F. J., Cheney, E. W and Lawson,
!       Charles, Wiley 1968.
!
!   [2] http://mathworld.wolfram.com/NormalDistribution.html
!
!   [3] `Computer Approximations', Hart, F. J., Cheney, E. W and Lawson,
!       Charles L., Wiley 1968 (Reprint 1978 w/corrections).
!
!   [4] Press, William H., Teukolsky, Saul A., Vetterling, William T.,
!       Flannery, Brian P., `Numerical Recipes in Fortran 77', second
!       edition, Cambridge University Press, 1992, ISBN 0-521-43064-X.
!
!   [5] `A precision approximation of the gamma function', Lanczos, C.,
!       SIAM Journal on Numerical Analysis, Ser. B., Vol. 1, pp 86--96, 1964.
!
!   [6] http://home.online.no/~pjacklam/notes/invnorm/
!
! Notes:
!   * Maybe the following books might also be useful:
!
!       `Elementary Functions, Algorithms and Implementation' by Jean-Michel
!       Muller, 2nd ed., 2006, XXII, 266 p. 36 illus., Hardcover,
!       ISBN: 0-8176-4372-9, Birkhäuser Boston.
!   
!   * `Software Manual for the Elementary Functions' by W. J. Cody, Jr. and
!      W. Waite, Prentice-Hall 1980, ISBN 0-13-822064-6.

module mod_special_functions

  use numeric_kinds
  use mod_constants

  private

  public :: phi
  public :: erfc
  public :: erf
  public :: gamma
  public :: norminv

  contains


    ! Normal distribution probabilities (Gaussian probability integral)
    ! accurate to 1d-15.  This is the cumulative distribution function of the
    ! standard normal distribution!
    !
    !   x = no. of standard deviations from the mean.
    !   y = probability to the left of x.  
    !   expntl = the probability density.
    !
    ! Based upon algorithm 5666 for the error function, from [1], page 292-293.
    ! erf and phi are related as follows:
    !
    !    erf(x) = 2*(phi(sqrt(2)*x)-0.5)
    !    or also
    !    phi(sqrt(2)*x) = 0.5 + 0.5*erf(x)
    !
    ! Adapted from an implementation by Alan Miller.
    !
    ! Latest revision - 30 March 1986
    !
    ! See also [2]
    !
    ! Note: another algorithm may be better and can be found at
    !
    ! http://groups.google.be/group/sci.math.num-analysis/browse_thread/thread/fea68d30def1de89/b53c1942ed97a236?lnk=st&q=erf+function+in+c+group%3Asci.math&rnum=1&hl=nl#b53c1942ed97a236
    !
    function phi(x) result (y)

      real(kind=qp), intent(in) :: x
      real(kind=qp)             :: y
      
      real(kind=qp)            :: xabs
      real(kind=qp)            :: expntl
      real(kind=qp), parameter ::   &
        p0 = 220.2068679123761_qp,  & ! 0.5*P00/sqrt(2)^0
        p1 = 221.2135961699311_qp,  & ! 0.5*P01/sqrt(2)^1
        p2 = 112.0792914978709_qp,  & ! 0.5*P02/sqrt(2)^2
        p3 = 33.91286607838300_qp,  & ! 0.5*P03/sqrt(2)^3
        p4 = 6.373962203531650_qp,  & ! 0.5*P04/sqrt(2)^4
        p5 = 0.7003830644436881_qp, & ! 0.5*P05/sqrt(2)^5
        p6 = 0.03526249659989109_qp   ! 0.5*P06/sqrt(2)^6
      real(kind=qp), parameter ::   &
        q0 = 440.4137358247522_qp,  & ! Q00/sqrt(2)^0
        q1 = 793.8265125199484_qp,  & ! Q01/sqrt(2)^1
        q2 = 637.3336333788311_qp,  & ! Q02/sqrt(2)^2
        q3 = 296.5642487796737_qp,  & ! Q03/sqrt(2)^3
        q4 = 86.78073220294608_qp,  & ! Q04/sqrt(2)^4
        q5 = 16.06417757920695_qp,  & ! Q05/sqrt(2)^5
        q6 = 1.755667163182642_qp,  & ! Q06/sqrt(2)^6
        q7 = 0.08838834764831844_qp   ! Q07/sqrt(2)^7
      real(kind=qp), parameter :: rootpi = 2.506628274631001_qp

      xabs = abs(x)

      ! |x| > 37
      if ( xabs > 37 ) then
        y = 0
      else

        expntl = exp(-xabs**2/2.0_qp)

        ! 0 <= |x| < 7
        if ( xabs < 7 ) then

          y = expntl*((((((p6*xabs + p5)*xabs + p4)*xabs + p3)*xabs  &
             + p2)*xabs + p1)*xabs + p0)/(((((((q7*xabs + q6)*xabs   &
             + q5)*xabs + q4)*xabs + q3)*xabs + q2)*xabs + q1)*xabs  &
             + q0)
        else

          ! 7 <= |x| <= 37
          y = expntl/(xabs + 1/(xabs + 2/(xabs + 3/(xabs + 4/  &
             (xabs + 0.65_qp)))))/rootpi
          
        end if

      end if

      if (x>0) then
        y = 1-y
      end if
        
    end function phi


    ! The complementary error function with precision 23.45, based on
    ! algorithm 5669 from [3], page 140 and 259-260.
    !
    ! Note that normally, the range is divided into two segments:
    ! 0 <= x <= 0.47 and 0.47 <= x < infinity.  For maximum accuracy
    ! other expressions than the ones here should be used within the
    ! range 0 <= x <= 0.47, but since [3] only gives the approximation
    ! for the range 0.47 <= x < infinity for convenience and simplicity,
    ! this is also what we use here.
    !
    elemental function erfc(x) result (y)
      real(kind=qp), intent(in) :: x
      real(kind=qp)             :: y

      real(kind=qp), parameter :: &
        P00 = 0.3762551871671301287283999844e+5_qp, &
        P01 = 0.7119025029224385863368043391e+5_qp, &
        P02 = 0.6792872095742636785330078287e+5_qp, &
        P03 = 0.411844859393225182268385918e+5_qp,  &
        P04 = 0.1722306325664886621380541951e+5_qp, &
        P05 = 0.5109973157155850523263765065e+4_qp, &
        P06 = 0.107157229472910126088636486e+4_qp,  &
        P07 = 0.1529918622334993008225057905e+3_qp, &
        P08 = 0.135107959589163542056048391e+2_qp,  &
        P09 = 0.5641895302730246255773946519e+0_qp, &
        Q00 = 0.3762551871671301287284013134e+5_qp, &
        Q01 = 0.1136461017633451072139295828e+6_qp, &
        Q02 = 0.1585390958920885698921632216e+6_qp, &
        Q03 = 0.1347344981315020759515349637e+6_qp, &
        Q04 = 0.773934045544434795855136954e+5_qp,  &
        Q05 = 0.3146456522263080448309756973e+5_qp, &
        Q06 = 0.919229873452846334367372023e+4_qp,  &
        Q07 = 0.1911284250203667015612859808e+4_qp, &
        Q08 = 0.2716711288118661667406454959e+3_qp, &
        Q09 = 0.2394725766908084566523751295e+2_qp, &
        Q10 = 0.1e+1_qp

      real(kind=qp) :: R_mn
      real(kind=qp) :: a

        ! Fix up for negative values (TODO: CHECK IF THIS IS OK.  I HAVE
        ! NOT CHECKED HOW NEGATIVE ARGUMENTS ARE HANDLED IN [3]).
        if (x<0) then
          a = -x
        else
          a = x
        end if

        R_mn = ((((((((((P09*a+P08)*a+P07)*a+P06)*a+P05)*a+P04)*a+P03)*a+P02)*a+P01)*a)+P00) &
                 /((((((((((((Q10*a+Q09)*a)+Q08)*a+Q07)*a+Q06)*a+Q05)*a+Q04)*a+Q03)*a+Q02)*a+Q01)*a)+Q00)

        y = exp(-a*a)*R_mn

        ! Fix up for negative values (TODO: CHECK IF THIS IS OK.  I HAVE
        ! NOT CHECKED HOW NEGATIVE ARGUMENTS ARE HANDLED IN [3]).
        if (x<0) then
          ! Apply erfc(-x) = 2-erfc(x)
          y = 2-y
        end if

    end function erfc


    ! The error function with precision 23.45, based on
    ! algorithm 5669 from [3], page 140 and 259-260.
    !
    ! Note that normally, the range is divided into two segments:
    ! 0 <= x <= 0.47 and 0.47 <= x < infinity.  For maximum accuracy
    ! other expressions than the ones here should be used within the
    ! range 0 <= x <= 0.47, but since [3] only gives the approximation
    ! for the range 0.47 <= x < infinity for convenience and simplicity,
    ! this is also what we use here.  The user should be made aware
    ! that the relative error of erf(x) for x near zero may be
    ! large.
    !
    function erf(x) result (y)
      real(kind=qp), intent(in) :: x
      real(kind=qp)             :: y

      y = 1-erfc(x)

    end function erf


    ! The gamma function.
    ! Using the approximation as described in [5] and implemented in [4] (but
    ! with slight changes).
    !
    ! TODO: see if this is the best method to compute the gamma function,
    !       check if better methods exist.
    !
    function gamma(xx) result (res)
      real(kind=qp), intent(in) :: xx
      real(kind=qp)             :: res

      integer(kind=i4b)             :: j
      real(kind=qp)                 :: ser, stp, tmp, x, y
      real(kind=qp), dimension(1:6) :: cof

      cof = (/ 76.18009172947146_qp, -86.50532032941677_qp, &
               24.01409824083091_qp, -1.231739572450155_qp, &
               0.001208650973866179_qp, -0.000005395239384953_qp /)
      stp = 2.5066282746310005_qp

      x = xx
      y = x
      tmp = x+5.5_qp
      tmp = (x+0.5_qp)*log(tmp)-tmp
      ser = 1.000000000190015_qp
      do j=1,6
        y = y + 1.0_qp
        ser = ser + cof(j)/y
      end do
      res = tmp + log(stp*ser/x)
      res = exp(res)

    end function gamma


    ! The inverse normal cumulative distribution function, based on
    ! Peter Jackal's algorithm with refinement with Halley's method.
    ! See [6].
    !
    elemental function norminv(x) result (res)
      real(kind=qp), intent(in) :: x
      real(kind=qp)             :: res

      real(kind=qp), dimension(1:6) :: a
      real(kind=qp), dimension(1:5) :: b
      real(kind=qp), dimension(1:6) :: c
      real(kind=qp), dimension(1:4) :: d
      real(kind=qp)                 :: x_low, x_high, q, r, e, u

      ! Coefficients in the rational approximations.

      a = (/ -3.969683028665376e+01_qp, &
              2.209460984245205e+02_qp, &
             -2.759285104469687e+02_qp, &
              1.383577518672690e+02_qp, &
             -3.066479806614716e+01_qp, &
              2.506628277459239e+00_qp /)

      b = (/ -5.447609879822406e+01_qp, &
              1.615858368580409e+02_qp, &
             -1.556989798598866e+02_qp, &
              6.680131188771972e+01_qp, &
             -1.328068155288572e+01_qp /)

      c = (/ -7.784894002430293e-03_qp, &
             -3.223964580411365e-01_qp, &
             -2.400758277161838e+00_qp, &
             -2.549732539343734e+00_qp, &
              4.374664141464968e+00_qp, &
              2.938163982698783e+00_qp /)

      d = (/ 7.784695709041462e-03_qp, &
             3.224671290700398e-01_qp, &
             2.445134137142996e+00_qp, &
             3.754408661907416e+00_qp /)

      ! Define break-points.

      x_low = 0.02425_qp
      x_high = 1 - x_low

      if ((0 < x) .and. (x < x_low)) then

        ! Rational approximation for lower region.
        q = sqrt(-2*log(x))
        res = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
               ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1)

      else if ((x_low <= x) .and. (x <= x_high)) then

        ! Rational approximation for central region
        q = x - 0.5_qp
        r = q*q
        res = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
              (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1)

      else if ((x_high < x) .and. (x < 1)) then

        ! Rational approximation for upper region.
        q = sqrt(-2*log(1-x))
        res = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
                ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1)

      end if

      ! The relative error of the approximation has absolute value less
      ! than 1.15*10^{-9}.  One iteration of Halley's rational method
      ! (third order) gives full machine precision.
      if ((0 < x) .and. (x < 1)) then
        e = 0.5_qp*erfc(-res/sqrt(2.0_qp)) - x
        u = e*sqrt(2*pi)*exp(0.5_qp*res*res)
        res = res - u/(1+res*0.5_qp*u)
      end if

    end function norminv


end module mod_special_functions
