! Module containing some number-theoretical functions.
!
! References:
!
!   [1] `A Course in Computational Algebraic Number Theory', Cohen, Henri,
!       Graduate Texts in Mathematics, Springer, 2000, ISBN 3-540-55640-0.
!
!   [2] `A survey of quadratic and inversive congruential pseudorandom
!        numbers', Juergen Eichenauer-Herrmann, Eva Herrmann and Stefan
!        Wegenkittl, in `Proceedings of the Second International Conference on
!        Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing',
!        Salzburg, July 9--12, 1996 of Lecture Notes in Statistics.
!        Springer-Verlag, New York, 1997, pages 67-97.
!
!   [3] `The Art of Computer Programming, Volume 2, Seminumerical Algorithms',
!       Donald E. Knuth, 1998, IBSN 0201896842.
!
module mod_number_theory

  use numeric_kinds
  use mod_primes

  private

  public :: integer_factorial
  public :: integer_factorial_recursive
  public :: integer_factorial_nonrecursive
  public :: real_factorial 
  public :: double_factorial
  public :: euler_totient
  public :: gcd
  public :: primitive_root
  public :: primitive_root_prime
  public :: inverse_mod
  public :: factor
  public :: nb_distinct_prime_factors
  public :: power_mod
  public :: cfrac_rational
  public :: warnock_cfrac

  interface integer_factorial
    module procedure integer_factorial_recursive
  end interface integer_factorial

  interface power_mod
    module procedure power_mod_simple
  end interface power_mod

  private :: power_mod_simple
  private :: power_mod_fast

contains


  ! Returns the factorial of the integer n.  Computations are done
  ! recursively and with integers.
  !
  ! Note:
  !   * Watch out with this routine, it only works for small values of n due
  !     to the fact that it works with integers!  For 4-byte integers,
  !     the maximum representable integer is 2147483648 and thus we can
  !     only compute until factorial(12) = 479001600.  From factorial(13) on
  !     we have overflow.
  !
  recursive function integer_factorial_recursive(n) result (res)
    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b)             :: res

    if (n < 1) then
      res = 1
    else
      res = n*integer_factorial_recursive(n-1)
    end if

  end function integer_factorial_recursive


  ! Returns the factorial of the integer n.  Computations are not
  ! done recursively and only with integers.
  !
  ! Note:
  !   * Watch out with this routine, it only works for small values of n due
  !     to the fact that it works with integers!  For 4-byte integers,
  !     the maximum representable integer is 2147483648 and thus we can
  !     only compute until factorial(12) = 479001600.  From factorial(13) on
  !     we have overflow.
  !
  function integer_factorial_nonrecursive(n) result (res)
    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b)             :: res

    integer(kind=i4b) :: i

    if (n < 1) then
      res = 1
    else
      res = product( (/ (i, i=2,n) /) )
    end if

  end function integer_factorial_nonrecursive


  ! Returns the factorial of an integer n as a variable of real kind.
  ! Computations are done with variables of real kind.
  !
  function real_factorial(n) result(res)
    integer(kind=i4b), intent(in) :: n
    real(kind=qp)                 :: res

    integer(kind=i4b) :: i

    res = 1.0_qp
    do i = 1,n
      res = res*i
    end do

  end function real_factorial


  ! Returns the double factorial of an integer n.
  !
  ! See http://mathworld.wolfram.com/DoubleFactorial.html
  !
  recursive function double_factorial(n) result (res)

    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b)             :: res

    if (n==-1 .or. n==0) then
      res = 1
    else
      res = n*double_factorial(n-2)
    end if

  end function double_factorial


  ! Return the Euler totient function of n.  This is equal to the
  ! number of positive integers less than or equal to n that are
  ! coprime to n.
  !
  ! Note:
  !   This algorithm is O(n).  Check if it can be done faster.
  !
  elemental function euler_totient(n) result (phi)
    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b)             :: phi

    integer(kind=i4b) :: i

    phi = 1
    do i = 2,n-1
      if (gcd(i, n) == 1) then
        phi = phi + 1
      end if
    end do

  end function euler_totient


  ! Compute the Greatest Common Divisor of a and b using the Euclidian
  ! algorithm (iterative form) as described in Algorithm 1.3.1 on
  ! page 12 of [1].  If either a or b is less than a given number N, the 
  ! number of Euclidian steps in this algorithm is bounded by a constant
  ! times ln(N), in both the worst case and on average.
  !
  elemental function gcd(a, b) result (res)
    integer(kind=i4b), intent(in) :: a, b
    integer(kind=i4b)             :: res

    integer(kind=i4b) :: x, y, r

    x = a
    y = b

    do

      if ( y == 0) then
        res = x
        exit
      end if

      r = modulo(x, y)
      x = y
      y = r

    end do

  end function gcd


  ! Given an odd prime p, this subroutine finds a primitive root modulo p.
  ! The routine is based on Algorithm 1.4.4 on page 25 of [1].
  !
  subroutine primitive_root_prime(p, r)
    integer(kind=i4b), intent(in)  :: p
    integer(kind=i4b), intent(out) :: r

    integer(kind=i4b)                            :: a, e
    integer(kind=i4b)                            :: i
    integer(kind=i4b)                            :: k
    integer(kind=i4b), dimension(:), allocatable :: factors, multiplicities

    a = 1
    k = nb_distinct_prime_factors(p-1)
    allocate(factors(k), multiplicities(k))
    call factor(p-1, factors, multiplicities)
    do

      a = a + 1
      i = 1

      do

        e = power_mod(a, (p-1)/factors(i), p)
        if (e == 1) then
          exit
        else
          i = i + 1
          if (i > k) then
            r = a
            return
          end if
        end if

      end do

    end do
    deallocate(factors, multiplicities)

  end subroutine primitive_root_prime


  ! Find *a* primitive root modulo n.
  !
  subroutine primitive_root(n, r)
    integer(kind=i4b), intent(in)  :: n
    integer(kind=i4b), intent(out) :: r

    integer(kind=i4b)                            :: g, g1, p, a
    integer(kind=i4b)                            :: nb_factors
    integer(kind=i4b), dimension(:), allocatable :: factors, multiplicities

    if (n == 2 .or. n == 4) then
      r = n-1
    else

      ! Factor n
      nb_factors = nb_distinct_prime_factors(n)
      allocate(factors(nb_factors), multiplicities(nb_factors))
      call factor(n, factors, multiplicities)


      ! n is of the form 2^a with a >= 3.

      if ( (size(factors) == 1) .and. (factors(1) == 2) &
             .and. (multiplicities(1) > 2) ) then


        r = 5
        return


      ! n is of the form p^a with a >=2.

      else if ( (size(factors) == 1) .and. (multiplicities(1) > 1) ) then


        p = factors(1)

        ! Compute a primitive root modulo p
        call primitive_root_prime(p, g)
        g1 = power_mod(g, p-1, p*p)
        if (g1 /= 1) then
          r = g
        else
          r = g + p ! modulo necessary???
        end if
        return


      ! n is of the form 2p^a with p an odd prime.
        
      else if ( (size(factors) == 2)                        &
                   .and. (factors(1) == 2)                  &
                     .and. (multiplicities(1) == 1) ) then


        p = factors(2)
        a = multiplicities(2)

        ! Compute g as a primitive root modulo p, it will also be a primitive
        ! root modulo p^a for a >= 2, because:
        !
        !   euler_totient(p^a) = (p-1)*p^{a-1}
        !
        !   g^{euler_totient(p^a)} = g^{(p-1)*p^{a-1}} = (g^{(p-1)})^{p^{a-1}}
        !                                              = 1^{p^{a-1}}
        !
        call primitive_root_prime(p, g)
        if (modulo(g, 2) == 1) then
          r = g
        else
          !r = g + power_mod(p, a, 2*p**a) ! necessary???
          r = g + p**a
        end if
        return

      else

        write(unit=*, fmt="(A)") &
          "ERROR: n is not of the correct form in primitive_root(n)!"
        stop

      end if

    end if

  end subroutine primitive_root


  ! Find the multiplicative inverse of b modulo m using the Extended Euclid
  ! Algorithm 1.3.6 as described on page 16-18 of [1].
  !
  ! Notes:
  !
  !   * This algorithm computes (u, v, d) such that b*u + m*v = d = gcd(b,m).
  !     The computations for v are commented out because we do not actually
  !     need it.
  !   * If d>1, then b is not invertible modulo m.
  !   * Zero does not have a multiplicative inverse, so z must
  !     be larger than zero!
  !
  elemental function inverse_mod(b, m) result(res)
    integer(kind=i4b), intent(in) :: b
    integer(kind=i4b), intent(in) :: m
    integer(kind=i4b)             :: res

    integer(kind=i4b) :: u, d
    !integer(kind=i4b) :: v
    integer(kind=i4b) :: q, v1, v3, t1, t3


    ! Initialize
    u = 1
    d = b
    if (m == 0) then
      !v = 0
      res = u
    else
      v1 = 0
      v3 = m
      do
        ! Finished?
        if (v3 == 0) then
          !v = (d-b*u)/m
          res = modulo(u, m) ! TODO: check if modulo necessary here
          exit
        else ! Euclidian step
          q = d/v3
          t3 = modulo(d, v3)
          t1 = u - q*v1
          u = v1
          d = v3
          v1 = t1
          v3 = t3
        end if
      end do
    end if

  end function inverse_mod


  ! Compute g^n modulo m using the naive but slow method that requires n-1
  ! group multiplications.  This method however does not suffer as much
  ! from the integer overflow problem as the fast method.
  !
  function power_mod_simple(g, n, m) result(res)
    integer(kind=i4b), intent(in) :: g, n, m
    integer(kind=i4b)             :: res

    integer(kind=i4b) :: i
    integer(kind=i4b) :: z

    if (n == 0) then
      res = 1
      return
    else if (n < 0) then
      ! g^{-n} mod m = (g^{-1})^n mod m
      z = inverse_mod(g, m)
    else
      z = g
    end if
    
    res = z
    do i = 2, abs(n)
      res = modulo(res*z, m)
    end do

  end function power_mod_simple


  ! Compute g^n modulo m using the Powering Algorithm 1.2.1 (Right-Left Binary)
  ! as described on page 8 in [1].  See also Knuth section 4.6.3???
  !
  ! WARNING: Watch out for integer overflow with this algorithm!
  !
  function power_mod_fast(g, n, m) result(res)
    integer(kind=i4b), intent(in) :: g, n, m
    integer(kind=i4b)             :: res

    integer(kind=i4b) :: y
    integer(kind=i4b) :: Nn
    integer(kind=i4b) :: z

    ! Some initializations
    y = 1
    if (n == 0) then
      res = y
      return
    else if (n < 0) then
      Nn = -n
      z = inverse_mod(g, m)
    else
      Nn = n
      z = g
    end if

    do

      ! Check the current bit of n to see if we need to multiply our
      ! accumulated result with the appropriate power of g.
      if (modulo(Nn, 2) == 1) then
        y = modulo(z*y, m) ! WARNING: watch out for integer overflow here!
      end if

      ! Halve N.
      Nn = Nn/2

      ! Have we checked all bits?
      if (Nn == 0) then

        res = y
        return

      else

        ! Update the current power of g that we should multiply our accumulated
        ! result with.
        z = modulo(z*z, m) ! WARNING: watch out for integer overflow here!

      end if

    end do

  end function power_mod_fast
 

  ! http://mathworld.wolfram.com/DistinctPrimeFactors.html
  ! http://mathworld.wolfram.com/PrimeFactorizationAlgorithms.html
  !
  ! Prime factorization.
  !
  subroutine factor(n, factors, multiplicities)
    integer(kind=i4b), intent(in)                :: n
    integer(kind=i4b), dimension(:), intent(out) :: factors
    integer(kind=i4b), dimension(:), intent(out) :: multiplicities

    integer(kind=i4b)                            :: temp
    integer(kind=i4b)                            :: i
    integer(kind=i4b), dimension(:), allocatable :: p
    integer(kind=i4b)                            :: factorindex

    temp = n

    allocate(p(nb_primes_up_to(n)))
    p = primes(nb_primes_up_to(n))
    factorindex = 0

    do i = 1, size(p)

      if (modulo(temp, p(i)) == 0) then

        factorindex = factorindex + 1
        factors(factorindex) = p(i)
        multiplicities(factorindex) = 1
        do
          temp = temp/p(i)
          if (modulo(temp, p(i)) /= 0) then
            exit
          end if
          multiplicities(factorindex) = multiplicities(factorindex) + 1
        end do

      end if

    end do

  end subroutine factor


  ! Return the number of prime factors of n.
  !
  function nb_distinct_prime_factors(n) result (res)
    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b)             :: res

    integer(kind=i4b) :: nb_primes
    integer(kind=i4b) :: i, x
    integer(kind=i4b), dimension(:), allocatable :: p

    x = n
    res = 0

    ! We have to check if each prime up to n is a possible divisor...
    nb_primes = nb_primes_up_to(n)
    allocate(p(nb_primes))
    p = primes(nb_primes)

    do i = 1, nb_primes

      ! If we have found a factor...
      if (modulo(x, p(i)) == 0) then

        ! ... increase the factor count...
        res = res + 1
        do

          ! ... and then factor it away
          x = x/p(i)
          if (modulo(x, p(i)) /= 0) then
            exit
          end if

        end do

      end if

    end do

    deallocate(p)

  end function nb_distinct_prime_factors


  ! Compute the continued fraction of numerator/denominator, where
  ! numerator > denominator >= 0.  This will always be finite, so
  ! we can also compute the sum of the partial quotients.
  !
  ! TODO: find a way to return the partial quotients instead of their sum.
  !
  subroutine cfrac_rational(numerator, denominator, sum_of_partial_quotients)
    integer(kind=i4b), intent(in)  :: numerator, denominator
    integer(kind=i4b), intent(out) :: sum_of_partial_quotients

    integer(kind=i4b) :: u, v
    integer(kind=i4b) :: a, old_u

    sum_of_partial_quotients = 0
    v = numerator
    u = denominator
    do
      if (v == 0) then
        exit
      end if

      ! Compute the next partial quotient.  Note that we rely
      ! on integer division here.
      a = u/v

      sum_of_partial_quotients = sum_of_partial_quotients + a

      old_u = u

      u = v
      v = modulo(old_u, v)

    end do

  end subroutine cfrac_rational


  ! Compute the continued fraction expansion of TOP/BOTTOM using
  ! a variant of Euclid's GCD algorithm.  See also [3], page 359.
  ! Note that we use the shorter representation that does not set the
  ! final term to 1 (See http://en.wikipedia.org/wiki/Continued_fraction)
  !
  ! NOTE: this was a subroutine sent to us by Tony Warnock.  Credits
  !       should go to him!
  !
  pure subroutine warnock_cfrac(numerator, denominator, maxquo, total, gcd)
    integer(kind=i4b), intent(in)  :: numerator   ! input numerator
    integer(kind=i4b), intent(in)  :: denominator ! input denominator
    integer(kind=i4b), intent(out) :: maxquo      ! maximum partial quotient
    integer(kind=i4b), intent(out) :: total       ! total partial quotients
    integer(kind=i4b), intent(out) :: gcd         ! greatest common divisor

    integer(kind=i4b) :: top    ! running numerator
    integer(kind=i4b) :: bottom ! running divisor
    integer(kind=i4b) :: rem    ! running remainder
    integer(kind=i4b) :: quo    ! running quotient

    total = 0
    maxquo = 0
    top = max(numerator, denominator)
    bottom = min(numerator, denominator)

    do
      if (bottom == 0) then
        exit
      end if
      quo = top/bottom
      rem = top-bottom*quo

      ! Keep track of the sum of the partial quotients.
      total = total + quo

      ! Check what is the maximum partial quotient.
      maxquo = max(maxquo, quo)

      top = bottom
      bottom = rem
    end do
    gcd = top

  end subroutine warnock_cfrac

end module mod_number_theory
