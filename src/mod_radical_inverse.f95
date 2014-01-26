! Module to implement different kinds of radical inverse functions.
!
! References:
!
!   [1] `An improved low-discrepancy sequence for multidimensional quasi-Monte
!       Carlo integration', Braaten, Eric and Weller, George, Journal of
!       Computational Physics, vol. 33, 1979, pages 249--258.
!
!   [2] `On the optimal Halton sequence', H. Chi, M. Mascagni, T. Warnock.
!
!   [3] `Generating and Testing the Modified Halton Sequences', Emanouil I.
!       Atanassov and Mariya K. Durchova in `Revised Papers from the 5th
!       International Conference on Numerical Methods and Applications',
!       Lecture Notes in Computer Science 2542, pp. 91-98, 2003.
!
!   [4] Warnock, Tony T., `Computational investigations of low-discrepancy point
!       sets.' in S. K. Zaremba, editor, `Applications of Number Theory to
!       Numerical Analysis', (Proc. Sympos., Univ. Montreal, Que. 1971), pages
!       319-343.  Academic Press, New York, 1972.
!
!   [5] Warnock, Tony T., `Computational investigations of
!       low-discrepancy point sets II', Monte Carlo and quasi-Monte
!       Carlo methods in scientific computing, Springer, 1995. p. 354-361.
!
!   [6] http://www.mcqmc.org/listarchive/mcqmc/discussion/98/Nov/msg00034.html
!
!   [7] http://public.lanl.gov/kmh/uncertainty/meetings/warnvgr.pdf
!
!   [8] `A survey of quadratic and inversive congruential pseudorandom
!        numbers', Juergen Eichenauer-Herrmann, Eva Herrmann and Stefan
!        Wegenkittl, in `Proceedings of the Second International Conference on
!        Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing',
!        Salzburg, July 9--12, 1996 of Lecture Notes in Statistics.
!        Springer-Verlag, New York, 1997, pages 67-97.
!
! TODO:
!   * it would be more natural if permutations had indexes p(0:b-1) instead
!     of p(1:b).  Try to fit this in...

module mod_radical_inverse

  use numeric_kinds
  use mod_utilities
  use mod_number_theory

  private

  private :: get_spin_hardcoded
  public :: warnock_spin_get

  public :: get_digits
  public :: radical_inverse
  public :: radical_inverse_scrambled
  public :: radical_inverse_lin_scrambled
  public :: radical_inverse_modified
  public :: radical_inverse_phicf
  public :: radical_inverse_folded
  public :: radical_inverse_inversive
  public :: radical_inverse_quadratic

  contains

    
    ! Get the digits of n in base b and put them in mydigits.
    ! The least significant digit is stored in the first
    ! element of mydigits.
    !
    pure subroutine get_digits(n, b, mydigits)

      integer(kind=i4b), intent(in)                :: n
      integer(kind=i4b), intent(in)                :: b
      integer(kind=i4b), dimension(:), intent(out) :: mydigits

      integer(kind=i4b) :: i
      integer(kind=i4b) :: k

      mydigits = 0

      i = 1
      k = n
      do 
        mydigits(i) = modulo(k, b)
        k = k/b  ! we rely on integer division here!
        if (k == 0) then
          exit
        end if
        i = i+1
      end do

    end subroutine get_digits


    ! The standard radical inverse function.
    !
    elemental function radical_inverse(n, b) result (radinv)

      integer(kind=i4b), intent(in) :: n 
      integer(kind=i4b), intent(in) :: b 
      real(kind=qp)                 :: radinv

      integer(kind=i4b) :: k
      real(kind=qp)     :: b_factor

      radinv = 0.0_qp
      k = n
      b_factor = 1.0_qp
      do
        b_factor = b_factor/b
        radinv = radinv + modulo(k, b)*b_factor
        k = k/b  ! we rely on integer division here!
        if (k == 0) then
          exit
        end if
      end do

    end function radical_inverse


    ! The scrambled radical inverse function, using the permutation p.
    ! (See [1])
    !
    function radical_inverse_scrambled(n, b, p) result (radinv)

      integer(kind=i4b), intent(in)               :: n
      integer(kind=i4b), intent(in)               :: b
      integer(kind=i4b), dimension(0:), intent(in) :: p
      real(kind=qp)                               :: radinv

      integer(kind=i4b) :: k
      real(kind=qp)     :: b_factor

      radinv = 0.0_qp
      k = n
      b_factor = 1.0_qp
      do
        b_factor = b_factor/b
        radinv = radinv + p(modulo(k, b))*b_factor
        k = k/b ! we rely on integer division here!
        if (k == 0) then
          exit
        end if
      end do

    end function radical_inverse_scrambled


    ! Calculate the linear scrambled radical inverse function of n in base b.
    ! The scrambling of digitvalue x is defined by the linear function
    !
    !                   f(x) = w*x + c (mod b)
    !
    ! (See [2])
    !
    elemental function radical_inverse_lin_scrambled(n, b, w, c) result (radinv)

      integer(kind=i4b), intent(in) :: n
      integer(kind=i4b), intent(in) :: b
      integer(kind=i4b), intent(in) :: w
      integer(kind=i4b), intent(in) :: c
      real(kind=qp)                 :: radinv

      integer(kind=i4b) :: k
      real(kind=qp)     :: b_factor

      radinv = 0.0_qp
      k = n
      b_factor = 1.0_qp
      do
        b_factor = b_factor/b
        radinv = radinv + modulo(modulo(k,b)*w+c, b)*b_factor
        k = k/b ! we rely on integer division here!
        if (k == 0) then
          exit
        end if
      end do
 
    end function radical_inverse_lin_scrambled


    ! Atanassov's modified radical inverse function. The scrambling of
    ! digitvalue x is defined by:
    !
    !                   f(x) = (x*k^(j)   + scr)  mod b    (not increased)
    !
    !                   f(x) = (x*k^(j+1) + scr)  mod b    (increased)
    !
    ! If scrambled is set to .true., then scr will be a random digit from the
    ! range 0..b-1.
    !
    ! (See [3])
    !
    subroutine radical_inverse_modified(n, b, k, increased, scrambled, radinv)

      integer(kind=i4b), intent(in) :: n
      integer(kind=i4b), intent(in) :: b
      integer(kind=i4b), intent(in) :: k
      logical, intent(in), optional :: increased
      logical, intent(in), optional :: scrambled
      real(kind=qp), intent(out)    :: radinv

      integer(kind=i4b) :: temp
      real(kind=qp)     :: b_factor
      integer(kind=i4b) :: k_factor
      integer(kind=i4b) :: scr
      logical           :: l_increased
      logical           :: l_scrambled

      if (present(increased)) then
        l_increased = increased
      else
        l_increased = .true.
      end if

      if (present(scrambled)) then
        l_scrambled = scrambled
      else
        l_scrambled = .false.
      end if

      radinv = 0.0_qp

      if (l_increased) then
        k_factor = k
      else
        k_factor = 1
      end if

      if (l_scrambled) then
        call random_integer(b, scr)
      else
        scr = 0
      end if

      temp = n
      b_factor = 1.0_qp
      do
        b_factor = b_factor/b
        radinv = radinv + modulo(modulo(temp,b)*k_factor+scr, b)*b_factor
        temp = temp/b  ! we rely on integer division here!
        if (temp == 0) then
          exit
        end if
        k_factor = modulo(k*k_factor, b)
      end do

    end subroutine radical_inverse_modified


    ! Warnock's PhiCf radical inverse function.
    !
    ! (See [5], [6] and [7])
    elemental function radical_inverse_phicf(n, b) result (radinv)

      integer(kind=i4b), intent(in) :: n
      integer(kind=i4b), intent(in) :: b
      real(kind=qp)                 :: radinv

      integer(kind=i4b) :: spin

      ! Calculate it ourselves or using a table lookup...
      !spin = warnock_spin_get(b)
      spin = get_spin_hardcoded(b)

      radinv = radical_inverse_lin_scrambled(n, b, spin, 0)

    end function radical_inverse_phicf


    ! Warnock's folded radical inverse function.
    ! (See [4])
    !
    elemental function radical_inverse_folded(n, b) result (radinv)

      integer(kind=i4b), intent(in) :: n
      integer(kind=i4b), intent(in) :: b
      real(kind=qp)                 :: radinv

      integer(kind=i4b) :: i
      integer(kind=i4b) :: temp
      real(kind=qp)     :: b_factor

      temp = n
      radinv = 0.0_qp
      b_factor = 1.0_qp
      i = 0

      ! Transform the known digits before the all-zero end-digits.
      do
        b_factor = b_factor/b
        radinv = radinv + modulo(modulo(temp,b)+i, b)*b_factor
        temp = temp/b
        if (temp == 0) then
          exit
        end if
        i = i+1
      end do

      ! Transform the all-zero digits at the end until there's no noticable
      ! difference anymore.
      do
        i = i+1
        b_factor = b_factor/b
        if (modulo(i, b) /= 0  & ! (necessary for spacing intrinsic)
              .and. b_factor*modulo(i, b) < spacing(radinv)) then
          exit
        end if
        radinv = radinv + b_factor*modulo(i, b)
      end do

    end function radical_inverse_folded


    ! Calculate the inversive scrambled radical inverse function of n in base b.
    ! The scrambling of digitvalue x is defined by the linear function
    !
    !                   f(x) = w*x^{-1} + c (mod b)
    !
    ! (See [8])
    !
    elemental function radical_inverse_inversive(n, b, w, c) result (radinv)

      integer(kind=i4b), intent(in) :: n
      integer(kind=i4b), intent(in) :: b
      integer(kind=i4b), intent(in) :: w
      integer(kind=i4b), intent(in) :: c
      real(kind=qp)                 :: radinv

      integer(kind=i4b) :: k
      integer(kind=i4b) :: digit
      real(kind=qp)     :: b_factor

      radinv = 0.0_qp
      k = n
      b_factor = 1.0_qp
      do
        b_factor = b_factor/b
        digit = modulo(k,b)

        ! Only compute the multiplicative inverse if the digit is non-zero!
        if (digit /= 0) then
          radinv = radinv + modulo( w*inverse_mod(digit, b)+c, b )*b_factor
        end if

        k = k/b ! we rely on integer division here!
        if (k == 0) then
          exit
        end if
      end do

    end function radical_inverse_inversive


    ! Calculate the quadratic scrambled radical inverse function of n in base b.
    ! The scrambling of digitvalue x is defined by the linear function
    !
    !                   f(x) = u*x^2 + v*x + w (mod b)
    !
    ! (See [8])
    !
    elemental function radical_inverse_quadratic(n, b, u, v, w) result (radinv)

      integer(kind=i4b), intent(in) :: n
      integer(kind=i4b), intent(in) :: b
      integer(kind=i4b), intent(in) :: u
      integer(kind=i4b), intent(in) :: v
      integer(kind=i4b), intent(in) :: w
      real(kind=qp)                 :: radinv

      integer(kind=i4b) :: k
      integer(kind=i4b) :: digit, digit2
      real(kind=qp)     :: b_factor

      radinv = 0.0_qp
      k = n
      b_factor = 1.0_qp
      do
        b_factor = b_factor/b

        digit = modulo(k, b)
        digit2 = modulo(digit*digit, b)

        ! The lots of modulo here is slower, but more robust against overflow...
        radinv = radinv + modulo( modulo(u*digit2, b) &
                                    + modulo(v*digit, b) + w, b )*b_factor

        k = k/b ! we rely on integer division here!
        if (k == 0) then
          exit
        end if
      end do

    end function radical_inverse_quadratic


    ! Get Warnock's `Spin'-value for base b (in a hardcoded table lookup)
    !
    pure function get_spin_hardcoded(b) result (spin)

      integer(kind=i4b), intent(in) :: b
      integer(kind=i4b)             :: spin

      select case (b)
      ! Only spin-values for the first 31 primes are listed here...
      case (2)          ! prime 1
         spin = 1
      case (3)          ! prime 2
         spin = 2
      case (5)          ! prime 3
         spin = 2
      case (7)          ! prime 4
         spin = 5
      case (11)         ! prime 5
         spin = 3
      case (13)         ! prime 6
         spin = 8
      case (17)         ! prime 7
         spin = 3
      case (19)         ! prime 8
         spin = 7
      case (23)         ! prime 9
         spin = 18
      case (29)         ! prime 10
         spin = 12
      case (31)         ! prime 11
         spin = 18
      case (37)         ! prime 12
         spin = 4
      case (41)         ! prime 13
         spin = 17 ! old version
         !spin = 18 ! new version...
      case (43)         ! prime 14
         spin = 24
      case (47)         ! prime 15
         spin = 40
      case (53)         ! prime 16
         spin = 14
      case (59)         ! prime 17
         spin = 41
      case (61)         ! prime 18
         spin = 50
      case (67)         ! prime 19
         spin = 12
      case (71)         ! prime 20
         spin = 30
      case (73)         ! prime 21
         spin = 40
      case (79)         ! prime 22
         spin = 70
      case (83)         ! prime 23
         spin = 10
      case (89)         ! prime 24
         spin = 39
      case (97)         ! prime 25
         spin = 82
      case (101)        ! prime 26
         spin = 6
      case (103)        ! prime 27
         spin = 16
      case (107)        ! prime 28
         spin = 37
      case (109)        ! prime 29
         spin = 48
      case (113)        ! prime 30
         spin = 71
      case (127)        ! prime 31
         spin = 34
      case default
         spin = 1

      end select

    end function get_spin_hardcoded


    ! Compute the best continued fraction expansion around.
    ! First: get the smallest lowertotal
    ! Second: get smallest partial quotient (if tie)
    ! Third: get closer approximation
    !
    ! Note: this is a slightly adapted version of the FORTRAN 77 code Tony
    ! Warnock sent us.
    ! 
    function warnock_spin_get(b) result (res)
      integer(kind=i4b), intent(in) :: b   ! number being checked
      integer(kind=i4b)             :: res ! the spin value found

      integer(kind=i4b) :: gcd        ! greatest common divisor (not used here)
      integer(kind=i4b) :: lower      ! floor of root
      integer(kind=i4b) :: lowerquo   ! lower's maximum partial quotient
      integer(kind=i4b) :: lowertotal ! lower total of partial quotients
      integer(kind=i4b) :: upper      ! ceiling of root
      integer(kind=i4b) :: upperquo   ! upper's maximum partial quotient
      integer(kind=i4b) :: uppertotal ! upper's total of partial quotients
      real(kind=qp) :: zzspin         ! proximate spin


      ! Base 2 and 3 are special cases: if Spin becomes 2 or 3, then
      ! (Spin*a)mod3 is always 0 for these cases.  We solve this by setting
      ! Spin to 1 (b=2) and 2 (b=3)
      if (b == 2) then
        res = 1
        return
      end if
      if (b == 3) then
        res = 2
        return
      end if

      ! From personal communication with Tony Warnock, we know that it might be
      ! better to take the multiplier for the prime 41 to be 18 and not 17. This
      ! removes the closeness of 17/41 to 12/29.  According to Tony, no other
      ! such relationship exists up to base 99991.
      if ((b == 41)) then
        res = 18
        return
      end if

      ! Once we get here, it means we are not in a special case
      ! and we simply follow the algorithm to compute the spin values.
      zzspin = b*frac_part(sqrt(real(b, kind=qp)))
      lower = floor(zzspin, kind=i4b)
      upper = lower + 1

      call warnock_cfrac(lower, b, lowerquo, lowertotal, gcd)
      call warnock_cfrac(upper, b, upperquo, uppertotal, gcd)

      if (uppertotal < lowertotal) then
        res = upper
        return
      else if (uppertotal > lowertotal) then
        res = lower
        return
      else if (upperquo < lowerquo) then
        res = upper
        return
      else if (upperquo > lowerquo) then
        res = lower
        return
      else if (abs(upper-zzspin) > abs(lower-zzspin)) then
        res = lower
        return
      else
        res = upper
      end if 

    end function warnock_spin_get

end module mod_radical_inverse
