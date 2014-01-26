!  Module that implements various changes of variables, used as periodizing
!  transformations.
!
!  References:
!
!    [1] `Periodizing transformations for numerical integration', Dirk P.
!        Laurie, Journal of Computational and Applied Mathematics 66, 1996,
!        pages 337-344.
!
!    [2] `Lattice Methods for Multiple Integration', Sloan, I. H. and Joe, S.,
!        Oxford Science Publications, 1994, ISBN 0 19 853472 8.
!
!    [3] `A new variable transformation for numerical integration', Avram Sidi,
!        in 'Numerical Integration IV: Proceedings of the Conference at the
!        Mathematical Research Institute, Oberwolfach, November 8-14, 1992',
!        Birkhauser Verlag, Basel, 1993.
!
!    [4] `Transformation of integrands for lattice rules', Beckers Marc and
!        Haegemans, Ann in `Numerical Integration: Recent Developments,
!        Software and Applications', Eds. Espelid, T. and Genz, A. C.,
!        ISBN 0-7923-1583-9, pages 329-340, 1992.
!
!    [5] Mori, M., `An IMT-type double exponential formula for numerical
!        integration', Publ. Res. Inst. Math. Sci. Kyoto Univ., vol. 14,
!        pages 713-729, 1978.
!
!    [6] Sag, T. W. and Szekeres, G., `Numerical Evaluation of High-Dimensional
!        Integrals', Mathematics of Computation, Vol. 18, No. 86, april, 1964,
!        pages 245--253.
!
!    [7] Robinson, Ian and de Doncker, Elise, `Algorithm 45: Automatic
!        Computation of Improper Integrals over a Bounded or Unbounded Planar
!        Region', Computing, Vol. 27, 1981, pages 253--284.
!
!    [8] Haegemans, A., `Algorithm 34: An Algorithm for the Automatic
!        Integration Over a Triangle', Computing, Vol. 19, 1977, pages 179--187.
!
!    [9] Elise De Doncker and  R. Piessens, `Algorithm 32 - Automatic
!        Computation of Integrals with Singular Integrand, over a Finite or an
!        Infinite Range', Computing, Vol. 17, 1976, pages 265-279.
!
!   [10] Masao Iri, Sigeiti Moriguti and Yoshimitsu Takasawa, `On a certain
!        quadrature formula', Journal of Computational and Applied Mathematics,
!        Vol. 17, 1987, pages 3-20.
!
!   [11] Murota, Kazuo and Iri, Masao, `Parameter Tuning and Repeated
!        Application of the IMT-type Transformation in Numerical Quadrature',
!        Numerische Mathematik, vol. 38, pages 347-363, 1982.
!
!   [12] Takahasi, H. and Mori, M, `Double Exponential Formulas for Numerical
!        Integration', Publ. Res. Inst. Math. Sci. Kyoto University, vol. 9, 
!        pages 721-741, 1974.
!
!   [13] Lyness, J. N. and Delves, L. M., `On the implementation of a modified
!        Sag-Szekeres quadrature method', Journal of Computational and
!        Applied Mathematics, Vol. 112, pages 189-200, 1999.
!
!  TODO:
!
!    * Check why Mori's transformation often results in arithmetic exceptions
!      and check if this can be solved somehow.
!
!    * Add the 'method of complete symmetrization' by Zaremba (1972), which is
!      however computationally intensive...
!
!    * Add other transformations from [4].
!
!    * Do *NOT* add the method with the Bernouilli polynomials, because a
!      numerical implementation would require the computation of a lot of
!      partial derivatives by finite difference methods, and one has to be
!      very careful for the accuracy of the calculations (see also [4],
!      page 333 for this).

module mod_periodize

  use numeric_kinds
  use mod_function
  use mod_constants

  private

  public :: periform

  private :: f_none
  private :: f_sidi_m2
  private :: f_korobov_m1
  private :: f_korobov_m2
  private :: f_korobov_m3
  private :: f_korobov_m4
  private :: f_korobov_m5
  private :: f_korobov_m6
  private :: f_korobov_m7
  private :: f_korobov_m8
  private :: f_laurie_m2
  public :: f_mori
  public :: f_tanh
  public :: f_de
  public :: f_imt

  private :: f_imt_1d

  ! The maximum length for the name of a periodizer.
  integer, parameter, public :: MAX_LENGTH_PERIODIZER_NAME = 10

contains

  subroutine periform(f, params, x, method, y_periodized)

    interface
      function f(x, fparams) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function f
    end interface
    type(functionparams), intent(in), optional :: params
    real(kind=qp), dimension(:), intent(in)    :: x
    character(len=*), intent(in)               :: method
    real(kind=qp), intent(out)                 :: y_periodized

    real(kind=qp), dimension(size(x)) :: y
    real(kind=qp), dimension(size(x)) :: y_prime

    select case (method)

      case ("none")

        call f_none(x, y, y_prime)

      case ("Sidi_m2")
 
        call f_sidi_m2(x, y, y_prime)

      case ("Korobov_m1")

        call f_korobov_m1(x, y, y_prime)

      case ("Korobov_m2")

        call f_korobov_m2(x, y, y_prime)

      case ("Korobov_m3")

        call f_korobov_m3(x, y, y_prime)

      case ("Korobov_m4")

        call f_korobov_m4(x, y, y_prime)

      case ("Korobov_m5")

        call f_korobov_m5(x, y, y_prime)

      case ("Korobov_m6")

        call f_korobov_m6(x, y, y_prime)

      case ("Korobov_m7")

        call f_korobov_m7(x, y, y_prime)

      case ("Korobov_m8")

        call f_korobov_m8(x, y, y_prime)

      case ("Laurie_m2")

        call f_laurie_m2(x, y, y_prime)

      case ("Mori")

        call f_mori(x, y, y_prime)

      case ("TANH")

        call f_tanh(x, y, y_prime)

      case ("IMT")

        call f_imt(x, y, y_prime)

      case ("DE")
      
        call f_de(x, y, y_prime)

      case default

        print *, "ERROR: unknown periodizing method!"

    end select
      
    y_periodized = f(y, params)*product(y_prime)

  end subroutine periform


  ! Dummy transformation for the case of no periodization.
  ! 
  !   phi(t) = t
  !
  !   phi'(t) = 1
  !
  pure subroutine f_none(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t

    y_prime = 1.0_qp

  end subroutine f_none


  ! Sidi's trigonometric transformation for m=2 (See [3]).
  ! 
  !   phi(t) = t - sin(2*pi*t)/(2*pi)
  !
  !   phi'(t) = 1 - cos(2*pi*t)
  !
  pure subroutine f_sidi_m2(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t - (sin(2*pi*t))/(2*pi)

    y_prime = 1 - cos(2*pi*t)

  end subroutine f_sidi_m2


  ! Korobov's polynomial transformation for m=1.  This is also the first
  ! polynomial transformation from Sloan and Joe's book (See [2], page 33).
  !
  !                3      2    2
  !   phi(t) = -2 t  + 3 t  = t  (3 - 2 t)
  !
  !                 2
  !   phi'(t) = -6 t  + 6 t = 6 t (1 - t)
  !
  pure subroutine f_korobov_m1(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t*t*(3-2*t)

    y_prime = 6*t*(1-t)

  end subroutine f_korobov_m1


  ! Korobov's polynomial transformation for m=2.  This is also the second
  ! polynomial transformation from Sloan and Joe's book (See [2], page 33).
  !
  !               5       4       3
  !   phi(t) = 6 t  - 15 t  + 10 t
  !
  !                 4       3       2       2        2
  !   phi'(t) = 30 t  - 60 t  + 30 t  = 30 t  (t - 1)
  !
  pure subroutine f_korobov_m2(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = ((6*t-15)*t+10)*t**3

    y_prime = 30*t**2*(1-t)**2

  end subroutine f_korobov_m2


  ! Korobov's polynomial transformation for m=3, computed with symbolic
  ! manipulation software.
  !
  !                   7       6       5       4
  !   phi(t) =  - 20 t  + 70 t  - 84 t  + 35 t
  !
  !   
  !                    6        5        4        3
  !   phi'(t) = - 140 t  + 420 t  - 420 t  + 140 t
  !
  !               3
  !           =  t  (t ((420 - 140 t) t - 420) + 140)
  !
  pure subroutine f_korobov_m3(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t**4*(t*((70-20*t)*t-84)+35)

    y_prime = t**3*(t*((420-140*t)*t-420)+140)

  end subroutine f_korobov_m3


  ! Korobov's polynomial transformation for m=4, computed with symbolic
  ! manipulation software.
  !
  pure subroutine f_korobov_m4(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t**5*(t*(t*(t*(70*t-315)+540)-420)+126)

    y_prime = t**4*(t*(t*(t*(630*t-2520)+3780)-2520)+630)

  end subroutine f_korobov_m4


  ! Korobov's polynomial transformation for m=5, computed with symbolic
  ! manipulation software.
  !
  pure subroutine f_korobov_m5(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t**6*(t*(t*(t*((1386-252*t)*t-3080)+3465)-1980)+462)

    y_prime = t**5*(t*(t*(t*((13860-2772*t)*t-27720)+27720)-13860)+2772)

  end subroutine f_korobov_m5


  ! Korobov's polynomial transformation for m=6, computed with symbolic
  ! manipulation software.
  !
  pure subroutine f_korobov_m6(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t**7*(t*(t*(t*(t*(t*(924*t-6006)+16380)-24024)+20020)-9009)+1716)

    y_prime = t**6*(t*(t*(t*(t*(t*(12012* &
                           t-72072)+180180)-240240)+180180)-72072)+12012)

  end subroutine f_korobov_m6


  ! Korobov's polynomial transformation for m=7, computed with symbolic
  ! manipulation software.
  !
  pure subroutine f_korobov_m7(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t**8*(t*(t*(t*(t*(t*((25740-3432*t)* &
         t-83160)+150150)-163800)+108108)-40040)+6435)

    y_prime = t**7*(t*(t*(t*(t*(t*((360360-51480*t)* &
         t-1081080)+1801800)-1801800)+1081080)-360360)+51480)

  end subroutine f_korobov_m7


  ! Korobov's polynomial transformation for m=8, computed with symbolic
  ! manipulation software.
  !
  pure subroutine f_korobov_m8(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t**9*(t*(t*(t*(t*(t*(t*(t*(12870* &
         t-109395)+408408)-875160)+1178100)-1021020)+556920)- &
              175032)+24310)

    y_prime = t**8*(t*(t*(t*(t*(t*(t*(t*(218790* &
         t-1750320)+6126120)-12252240)+15315300)-12252240)+6126120)- &
              1750320)+218790)

  end subroutine f_korobov_m8


  ! Laurie's periodizing transformation for m=2 (See [1]).
  !
  !               3       5       6      7
  !   phi(t) = 7 t  - 21 t  + 21 t  - 6 t  
  !
  !                 2        4        5       6
  !   phi'(t) = 21 t  - 105 t  + 126 t  - 42 t
  !
  pure subroutine f_laurie_m2(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    y = t**3*(7-t**2*(21-t*(21-t*6)))

    y_prime = t**2*(21-t**2*(105-t*(126-42*t)))

  end subroutine f_laurie_m2


  ! Mori's IMT-type double exponential transformation (See [5] and also [4]).
  !
  !                              /  1         \
  !   phi(t) = 1/2 tanh(a sinh(b |------ - 1/t|)) + 1/2
  !                              \-t + 1      /
  !   phi(0) = 0
  !   phi(1) = 0.5
  !
  !   phi'(t) = ... quite complicated expression ...
  !
  ! Notes:
  !   * Use this transformation only in the interval [0,1).  If you plot
  !     it for a wider range, you will see why.
  !
  ! TODO:
  !   * Check if the robust method can be implemented in a cleaner way.
  !   * Check why near 0 some f_mori values are 1.
  ! 
  subroutine f_mori(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    integer(kind=i4b) :: i
    real(kind=qp)     :: a, b
    real(kind=qp)     :: sinh_arg, cosh_res, tanh_arg

    a = pi/2.0_qp
    b = pi/4.0_qp


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The version that tries to be robust agains arithmetic exceptions... !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,size(t)

      if (t(i) /= 1.0_qp .and. t(i) /= 0.0_qp) then

        sinh_arg = b*(1.0_qp/(1.0_qp-t(i)) - 1.0_qp/t(i))
        if (abs(sinh_arg) >= log(huge(1.0_qp)/2)) then
          tanh_arg = huge(1.0_qp)/(2*a)
        else
          if (abs(sinh(sinh_arg)) >= huge(1.0_qp)/(2*a) ) then
            tanh_arg = huge(1.0_qp)/2
          else
            tanh_arg = a*sinh(sinh_arg)
          end if
        end if
    
        y(i) = 0.5_qp*tanh( tanh_arg ) + 0.5_qp
    
        if (abs(sinh_arg) >= log(huge(1.0_qp)/2)) then
          cosh_res = huge(1.0_qp)/2
        else
          cosh_res = cosh(sinh_arg)
        end if

        y_prime(i)= 0.5_qp*(1 - tanh(tanh_arg)**2)*a*cosh_res*b*(1/(1 - t(i))**2 + 1/t(i)**2)

      else

        if (t(i) == 0.0_qp) then
          y(i) = 0.0_qp
          y_prime(i) = 0.0_qp ! CHECK THIS!
        else
          y(i) = 0.5_qp
          y_prime(i) = 0.0_qp
        end if

      end if

    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The non-robust version... !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    y = 0.5_qp*tanh( a*sinh(b*(1.0_qp/(1.0_qp-t) - 1.0_qp/t)) ) + 0.5_qp
!
    ! This expression can easily be found by hand or using Maple...
!    y_prime = 0.5_qp*(1 - tanh(a*sinh(b * (1/(1-t) - 1/t)))**2) &
!                               * a * cosh(b*(1/(1-t) - 1/t))    &
!                               * b * (1/(1 - t)**2 + 1/t**2)

  end subroutine f_mori


  ! The TANH-transformation proposed by Sag and Szekeres (See [6] and
  ! also [11] and [12]).
  !
  ! phi(t) = 0.5 + 0.5*tanh( -0.5*c*(1/t-1/(1-t))) 
  ! phi(0) = 0
  ! phi(1) = 1
  !
  ! phi'(t) = ... quite complicated expression ...
  ! phi'(0) = 0
  ! phi'(1) = 0
  ! (see also [13], page 192)
  !
  ! Note:
  !   * This implementation should normally properly handle
  !     the cases t=0 and t=1.
  !   * Some floating underflow might occur.  I'm not sure
  !     in how much this might influence numerical results.
  !
  subroutine f_tanh(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    real(kind=qp)                     :: c
    real(kind=qp), dimension(size(t)) :: tanh_temp

    ! Sag and Szekeres take c=1/2 in their paper [6] on page 245,
    ! while Lyness and Delves take c=2 in [13] on page 190.
    !c = 0.5_qp
    c = 2.0_qp

    y = 0
    y_prime = 0


    if (all(t /= 0) .and. all(t /= 1)) then ! There are no singularities at 0 or 1

      tanh_temp = tanh(-0.5*c*(1/t-1/(1-t)))

      y = 0.5 + 0.5*tanh_temp

      y_prime = -0.25*(1-tanh_temp**2)*c*(-1/t**2-1/(1-t)**2)

    else ! there are singularities at either 0 or 1

      ! If there is a singularity at 0, then both y and y_prime are 0,
      ! so we check if we have such a singularity.
      if (all(t /= 0.0_qp)) then ! if there's only a singularity at 1...
        y = 1 ! phi(1) = 1
        where (t /= 1)
          tanh_temp = tanh(-0.5*c*(1/t-1/(1-t)))
          y = 0.5 + 0.5*tanh_temp
        end where
      end if

    end if

  end subroutine f_tanh


  ! The Double Exponential transformation (See [12] and also [11]).
  !
  ! phi(t) = 0.5*tanh(0.5*Pi*sinh(t)) + 0.5
  !
  ! phi'(t) = 0.25*Pi*cosh(t)/(cosh(0.5*Pi*sinh(t)))^2
  !
  ! Note that this is a transformation which maps [0,1] onto [-infty, +infty].
  !
  subroutine f_de(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    real(kind=qp), dimension(size(t)) :: tanh_temp

    tanh_temp = tanh(0.5*pi*sinh(t))

    y = 0.5*tanh_temp + 0.5

    ! This can easily be found by hand or using Maple...
    y_prime = 0.25*pi*cosh(t)*(1-tanh_temp**2)

  end subroutine f_de


  ! The multi-dimensional IMT-transformation (See [7,8,9,10] for papers
  ! on the one-dimensional version of this transformation).
  !
  !                         t[i]
  !                        /
  !                       |                1
  !                       |      exp(- ---------) du
  !               s       |            u (1 - u)
  !            --------' /
  !           '  |  |      0
  !  phi(t) =    |  |    ---------------------------
  !              |  |        1
  !              |  |       /
  !             i = 1      |             1
  !                        |   exp(- ---------) du
  !                        |         u (1 - u)
  !                       /
  !                         0
  !
  subroutine f_imt(t, y, y_prime)
    real(kind=qp), dimension(:), intent(in)  :: t
    real(kind=qp), dimension(:), intent(out) :: y
    real(kind=qp), dimension(:), intent(out) :: y_prime

    integer(kind=i4b) :: i
    real(kind=qp)     :: fi, si

    do i=1,size(t)
      call f_imt_1d(t(i), fi, si)
      y(i) = fi
      y_prime(i) = si
    end do

  end subroutine f_imt


  ! The one-dimensional IMT-transformation, accurate to 16 digits
  ! (See [7,8,9,10] for papers on the one-dimensional version of this
  ! transformation).
  !
  !                            t
  !                           /
  !                          |             1
  !                          |   exp(- ---------) du
  !                          |         u (1 - u)
  !                         /
  !                           0
  !              phi(t) = ---------------------------
  !                            1
  !                           /
  !                          |             1
  !                          |   exp(- ---------) du
  !                          |         u (1 - u)
  !                         /
  !                           0
  !
  ! The implementation is based on the routine in [7].  Three adaptations
  ! have been made, namely
  !
  !   * First of all, we have adapted the routine from [7] so that it
  !     also works for arguments larger than zero.
  !
  !   * Secondly, we made sure that the function is normalized by multiplying
  !     with the appropriate normalizing constant.
  !
  !   * Lastly, since the routine from [7] is based on the interval [-1..1], we
  !     have to perform a substitution.  If PHI_THERE(t) is the routine
  !     from [7] and PHI_HERE(t) is our routine, then we have that
  !
  !             PHI_HERE(t) = PHI_THERE(2t-1)
  !
  subroutine f_imt_1d(w, fi, si)
    real(kind=qp), intent(in)  :: w
    real(kind=qp), intent(out) :: fi
    real(kind=qp), intent(out) :: si
  
    ! The normalizing constant as calculated with Maple with Digits := 30,
    ! namely
    !
    !                            1
    !                           /
    !                          |           4
    !                          |   exp(- ------) dv
    !                          |              2
    !                         /          1 - v
    !                           -1
    ! 
    real(kind=qp), parameter :: normalizing_constant = &
                                  0.14059716813219312478482541060707912e-1_qp

    ! Q contains the coefficients of the Chebyshev expansion of FI.  These
    ! coefficients are taken from [7]. To achieve greater than 16-digit
    ! accuracy, it is necessary to increase the number of digits in these
    ! coefficients and also the number of coefficients used in the Chebyshev
    ! series.  This requires changes to the dimension of Q, Q0300, Q0603,
    ! Q0806 and Q1008, to the values of ARG.
    ! Refer to appendix 3 of [8] for details.  The number of digits in A1(3)
    ! and A1(4) must also be increased.
    real(kind=qp), dimension(68), parameter :: Q =         &

         ! Q1008: Chebyshev coefficients for the approximation of fi(t) in
         ! the interval [-1, -0.8).  There are 19 of them, being Q(1:19).
      (/ 0.1307592530313668e-01_qp,  0.8580903946529252e-02_qp, &
         0.1969565496511163e-02_qp, -0.6785908047534174e-04_qp, &
         0.5157989084218512e-05_qp, -0.3254283578286669e-06_qp, &
         0.3009165432294816e-07_qp, -0.2853960792992456e-08_qp, &
         0.3125316663674964e-09_qp, -0.369071915636499e-10_qp,  &
         0.47226104267733e-11_qp,   -0.6452090127750e-12_qp,    &
         0.935240772542e-13_qp,     -0.142880027692e-13_qp,     &
         0.22888743282e-14_qp,      -0.3827976860e-15_qp,       &
         0.665889416e-16_qp,        -0.120097537e-16_qp,        &
         0.22395647e-17_qp,                                     &
         ! Q0806: Chebyshev coefficients for the approximation of fi(t) in
         ! the interval [-0.8, -0.6].  There are 15 of them, being Q(20:34).
         0.7530091699129213e-01_qp,  0.2215896190278894e-01_qp, &
         0.1517662771681320e-02_qp, -0.1374204580631322e-04_qp, &
         0.2501181491115358e-05_qp, -0.2491206397236787e-07_qp, &
         0.5430948732810670e-08_qp, -0.1406070367942514e-09_qp, &
         0.160397857131033e-10_qp,  -0.7706914621139e-12_qp,    &
         0.656131517151e-13_qp,     -0.43257167630e-14_qp,      &
         0.3459964512e-15_qp,       -0.264283638e-16_qp,        &
         0.21607480e-17_qp,                                     &
         ! Q0603: Chebyshev coefficients for the approximation of fi(t) in
         ! the interval (-0.6, -0.3].  There are 16 of them, being Q(35:50).
         0.2285305746758029E+00_qp,  0.5670149715650661e-01_qp, &
         0.3881611773394649e-02_qp,  0.1483758828946990e-03_qp, &
         0.1986416462810431e-04_qp,  0.1166284710859293e-05_qp, &
         0.1048168134503124e-06_qp,  0.6572636994171403e-08_qp, &
         0.5344588684897204e-09_qp,  0.335115172128537e-10_qp,  &
         0.26225503158527e-11_qp,    0.1606252080762e-12_qp,    &
         0.124685606769e-13_qp,      0.73251000e-15_qp,         &
         0.579829277e-16_qp,         0.31817034e-17_qp,         &
         ! Q0300: Chebyshev coefficients for the approximation of fi(t) in
         ! the interval (-0.3, 0].  There are 18 of them, being Q(51:68)
         0.5404466499651320E+00_qp,  0.1034761954005156E+00_qp, &
         0.9090551216566596e-02_qp,  0.9130257654553327e-03_qp, &
         0.1026685176623610e-03_qp,  0.1035244509501565e-04_qp, &
         0.1034925154934048e-05_qp,  0.1002050044392230e-06_qp, &
         0.9552837592432324e-08_qp,  0.8941334072231475e-09_qp, &
         0.825727197051832e-10_qp,   0.75255600668576e-11_qp,   &
         0.6783483380205e-12_qp,     0.605164239084e-13_qp,     &
         0.53497536298e-14_qp,       0.4689365198e-15_qp,       &
         0.407909118e-16_qp,         0.35230149e-17_qp         /)
    real(kind=qp), dimension(4) :: A1 = (/ 20.0_qp,   20.0_qp,       &
                                           40.0_qp/3, 40.0_qp/3 /)
    real(kind=qp), dimension(4) :: A2 = (/ 18.0_qp, 14.0_qp,         &
                                            6.0_qp,  2.0_qp  /)

    integer(kind=i4b), dimension(5) :: arg = (/ 1, 20, 35, 51, 69 /)

    integer(kind=i4b) :: interval, k, kfirst, klast, krev, ksum
    real(kind=qp)     :: expont, uflow, eoflow
    real(kind=qp)     :: fac, T1, T2, T3
    real(kind=qp)     :: x, y, z


    ! `eoflow' is about equal to, but must not exceed the Naperian logarithm
    !          of the machine's largest real number (it is 87D0 in [7]).
    ! `uflow'  is about equal to, but must not be less than the smallest
    !          positive real number representable by the machine (it
    !          is 1D-38 in [7]).
    eoflow = log(huge(1.0_qp)/2)
    uflow = 2*tiny(1.0_qp)


    ! Set the argument x for evaluation.  This is necessary because the
    ! code in [7] only calculates things for z<0 and over the interval
    ! [-1..1] instead of over [0..1].  Also, since
    !             PHI_HERE(t) = PHI_THERE(2t-1)
    ! we perform this substitution here.
    z = 2*w-1
    if (z <= 0) then
      x = z
    else
      x = -z
    end if

    ! Compute si
    fi = 0.0_qp
    si = 0.0_qp
    if (x > -1.0_qp) then ! Watch out for the singularity at -1...

      expont = 4.0_qp/(1.0_qp-x*x)
      if (expont <= eoflow) then ! Watch out for underflow...

        si = exp(-expont)

        ! Determine coefficients to be used in Chebyshev series
        interval = 4
        if (x < -0.3_qp) then
          interval = 3
        end if
        if (x < -0.6_qp) then
          interval = 2
        end if
        if (x < -0.8_qp) then
          interval = 1
        end if
        y = A1(interval)*x + A2(interval)
        kfirst = arg(interval)
        klast = arg(interval+1)-2

        ! Compute fi
        T2 = 0.0_qp
        T3 = Q(klast+1)
        ksum = kfirst + klast
        do k = kfirst,klast
          krev = ksum - k
          T1 = T2
          T2 = T3
          T3 = y*T2 - T1 + Q(krev)
        end do
        fac = ((T3-T1)/2.0_qp)/normalizing_constant
        if (fac > uflow/si) then
          fi = fac*si
        end if

        ! TODO: check if the multiplication by 2 is correct here.
        si = 2*si/normalizing_constant

      end if

    end if

    if (z > 0.0_qp) then
      ! fi changes to 1-fi for positive z. si stays the same.
      fi = 1.0_qp-fi
    end if

  end subroutine f_imt_1d


end module mod_periodize
