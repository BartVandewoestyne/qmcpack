program test_number_theory

  use qmcpack

  integer(kind=i4b)                            :: i, j
  integer                                      :: my_unit
  integer(kind=i4b), dimension(:), allocatable :: p
  integer(kind=i4b)                            :: temp
  integer(kind=i4b)                            :: mysum
  integer(kind=i4b)                            :: primroot
  integer(kind=i4b), dimension(:), allocatable :: factors
  integer(kind=i4b), dimension(:), allocatable :: multiplicities
  integer(kind=i4b)                            :: nb_factors

  call show_test_header("MOD_NUMBER_THEORY")
  call init_assert()


  call start_test("Testing integer_factorial_recursive(x)...")
  call assert(integer_factorial_recursive(-1), 1)
  call assert(integer_factorial_recursive(0), 1)
  call assert(integer_factorial_recursive(1), 1)
  call assert(integer_factorial_recursive(2), 2)
  call assert(integer_factorial_recursive(3), 6)
  call assert(integer_factorial_recursive(4), 24)
  call assert(integer_factorial_recursive(5), 120)
  call assert(integer_factorial_recursive(6), 720)
  call assert(integer_factorial_recursive(7), 5040)
  call assert(integer_factorial_recursive(8), 40320)
  call assert(integer_factorial_recursive(9), 362880)
  call assert(integer_factorial_recursive(10), 3628800)
  call assert(integer_factorial_recursive(11), 39916800)
  call assert(integer_factorial_recursive(12), 479001600)
!  call assert(integer_factorial_recursive(13), 6227020800) ! OVERFLOW?
  call stop_test()

  call start_test("Testing integer_factorial_nonrecursive(x)...")
  call assert(integer_factorial_nonrecursive(-1), 1)
  call assert(integer_factorial_nonrecursive(0), 1)
  call assert(integer_factorial_nonrecursive(1), 1)
  call assert(integer_factorial_nonrecursive(2), 2)
  call assert(integer_factorial_nonrecursive(3), 6)
  call assert(integer_factorial_nonrecursive(4), 24)
  call assert(integer_factorial_nonrecursive(5), 120)
  call assert(integer_factorial_nonrecursive(6), 720)
  call assert(integer_factorial_nonrecursive(7), 5040)
  call assert(integer_factorial_nonrecursive(8), 40320)
  call assert(integer_factorial_nonrecursive(9), 362880)
  call assert(integer_factorial_nonrecursive(10), 3628800)
  call assert(integer_factorial_nonrecursive(11), 39916800)
  call assert(integer_factorial_nonrecursive(12), 479001600)
!  call assert(integer_factorial_nonrecursive(13), 6227020800) ! OVERFLOW?
  call stop_test()

  call start_test("Testing double_factorial(x)...")
  call assert(double_factorial(-1), 1)
  call assert(double_factorial(0), 1)
  call assert(double_factorial(1), 1)
  call assert(double_factorial(2), 2)
  call assert(double_factorial(3), 3)
  call assert(double_factorial(4), 8)
  call assert(double_factorial(5), 15)
  call assert(double_factorial(6), 48)
  call assert(double_factorial(7), 105)
  call assert(double_factorial(8), 384)
  call assert(double_factorial(9), 945)
  call assert(double_factorial(10), 3840)
  call assert(double_factorial(11), 10395)
  call assert(double_factorial(12), 46080)
  call assert(double_factorial(13), 135135)
  call assert(double_factorial(14), 645120)
  call assert(double_factorial(15), 2027025)
  call assert(double_factorial(16), 10321920)
  call stop_test()


  call start_test("Testing gcd(a,b)...")
  call assert(gcd(12, 18), 6)
  call assert(gcd(-4, 14), 2)
  call assert(gcd(5, 0), 5)
  call assert(gcd(9, 28), 1)
  call assert(gcd(42, 56), 14)
  call assert(gcd(1071, 1029), 21)
  call stop_test()

  
  call start_test("Testing euler_totient(n)...")
  call assert(euler_totient(1), 1)
  call assert(euler_totient(2), 1)
  call assert(euler_totient(3), 2)
  call assert(euler_totient(4), 2)
  call assert(euler_totient(5), 4)
  call assert(euler_totient(6), 2)
  call assert(euler_totient(7), 6)
  call assert(euler_totient(8), 4)
  call assert(euler_totient(9), 6)
  call assert(euler_totient((/ 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 /)), &
                            (/  4, 10,  4, 12,  6,  8,  8, 16,  6, 18 /))
  call assert(euler_totient((/ 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 /)), &
                            (/  8, 12, 10, 22,  8, 20, 12, 18, 12, 28 /))
  call assert(euler_totient((/ 30, 31, 32, 33, 34, 35, 36, 37, 38, 39 /)), &
                            (/  8, 30, 16, 20, 16, 24, 12, 36, 18, 24 /))
  call assert(euler_totient((/ 40, 41, 42, 43, 44, 45, 46, 47, 48, 49 /)), &
                            (/ 16, 40, 12, 42, 20, 24, 22, 46, 16, 42 /))
  call assert(euler_totient((/ 50, 51, 52, 53, 54, 55, 56, 57, 58, 59 /)), &
                            (/ 20, 32, 24, 52, 18, 40, 24, 36, 28, 58 /))
  call assert(euler_totient((/ 60, 61, 62, 63, 64, 65, 66, 67, 68, 69 /)), &
                            (/ 16, 60, 30, 36, 32, 48, 20, 66, 32, 44 /))
  call assert(euler_totient((/ 70, 71, 72, 73, 74, 75, 76, 77, 78, 79 /)), &
                            (/ 24, 70, 24, 72, 36, 40, 36, 60, 24, 78 /))
  call assert(euler_totient((/ 80, 81, 82, 83, 84, 85, 86, 87, 88, 89 /)), &
                            (/ 32, 54, 40, 82, 24, 64, 42, 56, 40, 88 /))
  call assert(euler_totient((/ 90, 91, 92, 93, 94, 95, 96, 97, 98, 99 /)), &
                            (/ 24, 72, 44, 60, 46, 72, 32, 96, 42, 60 /))
  call assert(euler_totient(2*11**5), 146410)
  call get_unit(my_unit)
  call checked_open(my_unit, "euler_totient.dat", "write")
  do i=1,1000
    write(unit=my_unit, fmt="(I4, I5)") i, euler_totient(i)
  end do
  call checked_close(my_unit)
  call stop_test()


  call start_test("Testing primitive_root(n) for n = 2 and 4...")
  call primitive_root(2, primroot)
  call assert(primroot, 1)
  call primitive_root(4, primroot)
  call assert(primroot, 3)
  call stop_test()


  call start_test("Testing primitive_root(n) for n=p^k with 2<=k<=5 and for the first 5 primes...")
  allocate(p(5))
  p = primes(5)
  do i = 2, size(p) ! Skip 2 as a prime, because the primes have to be odd
                    ! for this subroutine.

    do j = 2, 5

      ! Get the primitive root.
      call primitive_root(p(i)**j, primroot)

      ! The multiplicative order of the primitive root is Euler's
      ! totient function.
      temp = power_mod(primroot, euler_totient(p(i)**j), p(i)**j)
      call assert(temp, 1)

    end do

  end do
  deallocate(p)
  call stop_test()


  call start_test("Testing primitive_root(n) for n=2p^k with 2<=k<=4 and for the first 5 primes...")
  allocate(p(5))
  p = primes(5)
  do i = 2, size(p)
    do j = 2, 4

      ! Get the primitive root
      call primitive_root(2*p(i)**j, primroot)

      ! The multiplicative order of the primitive root is Euler's
      ! totient function
      temp = power_mod(primroot, euler_totient(2*p(i)**j), 2*p(i)**j)
      call assert(temp, 1)

    end do
  end do
  deallocate(p)
  call stop_test()


!  call start_test("Testing primitive_root_prime(n)...")
!  allocate(p(1500))
!  p = primes(1500)
!  do i = 2, size(p) ! Skip 2 as a prime, because the primes have to be odd
!                    ! for this subroutine.
!
!    ! Get the primitive root.
!    call primitive_root_prime(p(i), primroot)
!
!    ! The multiplicative order of the primitive root is Euler's
!    ! totient function.
!    temp = power_mod(primroot, euler_totient(p(i)), p(i))
!
!    ! We must have arrived at 1.
!    call assert(temp, 1)
!
!  end do
!  deallocate(p)
!  call stop_test()


  call start_test("Testing inverse_mod()...")
  allocate(p(1500))
  p = primes(1500)
  do i = 1, size(p)
    do j = 1, p(i)-1
      call assert(modulo(j*inverse_mod(j,p(i)), p(i)), 1)
    end do
  end do
  deallocate(p)
  call assert(inverse_mod(40, 41), 40)
  call stop_test()


  call start_test("Testing power_mod()...")
  call assert(power_mod(1, 0, 5), 1)
  call assert(power_mod(3, 0, 10), 1)
  call assert(power_mod(2, 4, 5), 1)
  call assert(power_mod(2, 5, 5), 2)
  call assert(power_mod(3, 1500, 7), 1)
  call assert(power_mod(2, 39366, 3**10), 1)
  call assert(power_mod(783, 437289, 4637), 617)
  ! These ones overflow
  ! call assert(power_mod(161053, 146410, 2*11**5), 1)
  call stop_test()


  call start_test("Testing nb_distinct_prime_factors()...")
  call assert(nb_distinct_prime_factors(1), 0)
  call assert(nb_distinct_prime_factors(2), 1)
  call assert(nb_distinct_prime_factors(2), 1)
  call assert(nb_distinct_prime_factors(3), 1)
  call assert(nb_distinct_prime_factors(4), 1)
  call assert(nb_distinct_prime_factors(5), 1)
  call assert(nb_distinct_prime_factors(6), 2)
  call assert(nb_distinct_prime_factors(7), 1)
  call assert(nb_distinct_prime_factors(8), 1)
  call assert(nb_distinct_prime_factors(9), 1)
  call assert(nb_distinct_prime_factors(10), 2)
  call assert(nb_distinct_prime_factors(11), 1)
  call assert(nb_distinct_prime_factors(12), 2)
  call assert(nb_distinct_prime_factors(13), 1)
  call assert(nb_distinct_prime_factors(14), 2)
  call assert(nb_distinct_prime_factors(15), 2)
  call assert(nb_distinct_prime_factors(16), 1)
  call stop_test()


  call start_test("Testing factor(...)...")
  
  !nb_factors = nb_distinct_prime_factors(1)
  !allocate(factors(nb_factors), multiplicities(nb_factors))
  !call factor(1, factors, multiplicities)
  !call assert(factors, (/ 1 /) )
  !call assert(multiplicities, (/ 0 /) )
  !deallocate(factors, multiplicities)

  nb_factors = nb_distinct_prime_factors(2)
  allocate(factors(nb_factors), multiplicities(nb_factors))
  call factor(2, factors, multiplicities)
  call assert(factors, (/ 2 /) )
  call assert(multiplicities, (/ 1 /) )
  deallocate(factors, multiplicities)

  nb_factors = nb_distinct_prime_factors(3)
  allocate(factors(nb_factors), multiplicities(nb_factors))
  call factor(3, factors, multiplicities)
  call assert(factors, (/ 3 /) )
  call assert(multiplicities, (/ 1 /) )
  deallocate(factors, multiplicities)

  nb_factors = nb_distinct_prime_factors(123)
  allocate(factors(nb_factors), multiplicities(nb_factors))
  call factor(123, factors, multiplicities)
  call assert(factors, (/ 3, 41/) )
  call assert(multiplicities, (/ 1, 1 /) )
  deallocate(factors, multiplicities)

  nb_factors = nb_distinct_prime_factors(1234)
  allocate(factors(nb_factors), multiplicities(nb_factors))
  call factor(1234, factors, multiplicities)
  call assert(factors, (/ 2, 617 /) )
  call assert(multiplicities, (/ 1, 1 /) )
  deallocate(factors, multiplicities)

  nb_factors = nb_distinct_prime_factors(68600)
  allocate(factors(nb_factors), multiplicities(nb_factors))
  call factor(68600, factors, multiplicities)
  call assert(factors, (/ 2, 5, 7 /) )
  call assert(multiplicities, (/ 3, 2, 3 /) )
  deallocate(factors, multiplicities)

  nb_factors = nb_distinct_prime_factors(32767)
  allocate(factors(nb_factors), multiplicities(nb_factors))
  call factor(32767, factors, multiplicities)
  call assert(factors, (/ 7, 31, 151 /) )
  call assert(multiplicities, (/ 1, 1, 1 /) )
  deallocate(factors, multiplicities)

  call stop_test()


  ! Test against the values in Knuth's book (Volume 2) on page 368.
  call start_test("Testing cfrac_rational()...")
  call cfrac_rational(1, 29, mysum)
  call assert(mysum, 29)
  call cfrac_rational(2, 29, mysum)
  call assert(mysum, 16)
  call cfrac_rational(3, 29, mysum)
  call assert(mysum, 12)
  call cfrac_rational(4, 29, mysum)
  call assert(mysum, 11)
  call cfrac_rational(5, 29, mysum)
  call assert(mysum, 10)
  call cfrac_rational(6, 29, mysum)
  call assert(mysum, 10)
  call cfrac_rational(7, 29, mysum)
  call assert(mysum, 11)
  call cfrac_rational(8, 29, mysum)
  call assert(mysum, 8)
  call cfrac_rational(9, 29, mysum)
  call assert(mysum, 9)
  call cfrac_rational(10, 29, mysum)
  call assert(mysum, 12)
  call cfrac_rational(11, 29, mysum)
  call assert(mysum, 8)
  call cfrac_rational(12, 29, mysum)
  call assert(mysum, 8)
  call cfrac_rational(13, 29, mysum)
  call assert(mysum, 9)
  call cfrac_rational(14, 29, mysum)
  call assert(mysum, 16)
  call cfrac_rational(15, 29, mysum)
  call assert(mysum, 16)
  call cfrac_rational(16, 29, mysum)
  call assert(mysum, 9)
  call cfrac_rational(17, 29, mysum)
  call assert(mysum, 8)
  call cfrac_rational(18, 29, mysum)
  call assert(mysum, 8)
  call cfrac_rational(19, 29, mysum)
  call assert(mysum, 12)
  call cfrac_rational(20, 29, mysum)
  call assert(mysum, 9)
  call cfrac_rational(21, 29, mysum)
  call assert(mysum, 8)
  call cfrac_rational(22, 29, mysum)
  call assert(mysum, 11)
  call cfrac_rational(23, 29, mysum)
  call assert(mysum, 10)
  call cfrac_rational(24, 29, mysum)
  call assert(mysum, 10)
  call cfrac_rational(25, 29, mysum)
  call assert(mysum, 11)
  call cfrac_rational(26, 29, mysum)
  call assert(mysum, 12)
  call cfrac_rational(27, 29, mysum)
  call assert(mysum, 16)
  call cfrac_rational(28, 29, mysum)
  call assert(mysum, 29)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_number_theory
