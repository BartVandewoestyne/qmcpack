! Module to compute Bernoulli numbers and Bernoulli polynomials.
!
! References:
!
!  [1] http://en.wikipedia.org/wiki/Bernoulli_number
!
module mod_bernoulli

use numeric_kinds
use mod_utilities
use mod_constants
use mod_number_theory

private

public :: bernoulli
public :: bernoulli1
public :: bernoulli2
public :: bernoulli3
public :: bernoulli_poly
public :: bernoulli_poly1
public :: bernoulli_poly2

interface bernoulli
  module procedure bernoulli1
end interface bernoulli

interface bernoulli_poly
  module procedure bernoulli_poly1
end interface bernoulli_poly

real(kind=qp), dimension(0:20), parameter, public :: bernoulli_number = &
                                            (/      1.0_qp,       & ! B0
                                                   -0.5_qp,       & ! B1
                                                    1.0_qp/6,     & ! B2
                                                    0.0_qp,       & ! B3
                                                   -1.0_qp/30,    & ! B4
                                                    0.0_qp,       & ! B5
                                                    1.0_qp/42,    & ! B6
                                                    0.0_qp,       & ! B7
                                                   -1.0_qp/30,    & ! B8
                                                    0.0_qp,       & ! B9
                                                    5.0_qp/66,    & ! B10
                                                    0.0_qp,       & ! B11
                                                 -691.0_qp/2730,  & ! B12
                                                    0.0_qp,       & ! B13
                                                    7.0_qp/6,     & ! B14
                                                    0.0_qp,       & ! B15
                                                -3617.0_qp/510,   & ! B16
                                                    0.0_qp,       & ! B17
                                                43867.0_qp/798,   & ! B18
                                                    0.0_qp,       & ! B19
                                              -174611.0_qp/330 /)   ! B20


contains


  ! Return the nth Bernoulli number using a table lookup.
  !
  function bernoulli1(n) result (bn)
    integer(kind=i4b), intent(in) :: n
    real(kind=qp)                 :: bn

    bn = bernoulli_number(n)

  end function bernoulli1


  ! Compute the nth Bernoulli number, with n a non-negative integer.
  ! Use a recurrence relation to compute the numbers.
  !
  ! Note:
  !   The method currently implemented is probably not the best, as
  !   n cannot be too high or else we suffer from roundoff...
  !   A similar implementation is at
  !   http://pagesperso-orange.fr/jean-pierre.moreau/Fortran/mbernoa_f90.txt
  !
  subroutine bernoulli2(n, b)
    integer(kind=i4b), intent(in)              :: n
    real(kind=qp), dimension(0:), intent(out) :: b

    integer(kind=i4b), dimension(:), allocatable :: binomial_coef
    integer(kind=i4b)                            :: j

    if (n==0) then

      b = (/ 1.0_qp /)

    else if (n==1) then

      b = (/ 1.0_qp, -0.5_qp /)

    else      

      b(0) = 1.0_qp
      b(1) = -0.5_qp

      do j = 2,n

        if (modulo(j,2)==1) then

          ! Except for B_1, all Bernoulli numbers of odd index are 0.
          b(j) = 0.0_qp

        else

          ! Use the following recursion to compute the numbers
          !
          !                    n - 1
          !                    ====
          !                    \
          !                     >    B(k) binomial(n + 1, k)
          !                    /
          !                    ====
          !                    k = 0
          !           B(n) = - -----------------------------
          !                                n + 1
          !
          !
          allocate(binomial_coef(0:j+1))
          call binomial_coefficients(j+1, binomial_coef)
          b(j) = -sum(binomial_coef(0:j-1)*b(0:j-1))/binomial_coef(j)
          deallocate(binomial_coef)

        end if

      end do

    end if

  end subroutine bernoulli2


  ! Compute the Bernoulli numbers using a series expansion.
  ! See equation (23.2.16) in Abramowitz and Stegun.
  !
  ! This seems to perform worse than bernoulli2?
  !
  ! Note:
  !     The routine makes use of the factorial(n) function, which only works
  !     well for small integers!!!  Consequently, this routine will also only
  !     work for not so large integers!
  !
  function bernoulli3(n) result(bn)
    integer(kind=i4b), intent(in) :: n
    real(kind=qp)                 :: bn

    integer(kind=i4b) :: i
    real(kind=qp)     :: mysum, term

    if (n<0) then

      bn = 0.0_qp

    else if (n==0) then

      bn = 1.0_qp

    else if (n==1) then

      bn = -0.5_qp

    else if (modulo(n,2)==1) then

      bn = 0.0_qp

    else

      mysum = 0.0_qp
      i = 1
      do
        term = (1.0_qp/i)**n
        !if (term < epsilon(mysum)) then
        if (term < 1.0e-15_qp) then
          exit
        end if
        mysum = mysum + term
        i = i+1
      end do

      bn = 2.0_qp*mysum*real_factorial(n)/((2.0_qp*pi)**n)

      if (modulo(n,4)==0) then
        bn = -bn
      end if

    end if

  end function bernoulli3


  ! Evaluate the Bernoulli polynomial of degree n at x using the
  ! explicit formula
  !                     n
  !                    ====
  !                    \                        n - k
  !            B (x) =  >    b  binomial(n, k) x
  !             n      /      k
  !                    ====
  !
  function bernoulli_poly1(n, x) result(res)
    integer(kind=i4b), intent(in) :: n
    real(kind=qp), intent(in)     :: x
    real(kind=qp)                 :: res

    real(kind=qp)     :: fact
    integer(kind=i4b) :: i

    fact = 1.0_qp
    res = 1.0_qp
    do i = 1,n
      fact = fact*(n+1-i)/i
      res = res*x + fact*bernoulli(i)
    end do

  end function bernoulli_poly1


  ! Evaluate the Bernoulli polynomial of degree n at x using the
  ! explicit formula
  !                     n
  !                    ====
  !                    \                        n - k
  !            B (x) =  >    b  binomial(n, k) x
  !             n      /      k
  !                    ====
  !
  subroutine bernoulli_poly2(n, x, res)
    integer(kind=i4b), intent(in) :: n
    real(kind=qp), intent(in)     :: x
    real(kind=qp), intent(out)    :: res

    real(kind=qp)                 :: fact
    real(kind=qp), dimension(0:n) :: b
    integer(kind=i4b)             :: i

    ! Pre-compute all Bernoulli numbers.
    call bernoulli2(n, b)

    ! Now compute the polynomial.
    fact = 1.0_qp
    res = 1.0_qp
    do i = 1,n
      fact = fact*(n+1-i)/i
      res = res*x + fact*b(i)
    end do

  end subroutine bernoulli_poly2


end module mod_bernoulli
