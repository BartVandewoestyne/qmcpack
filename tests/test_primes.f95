! Tests various procedures used for and with primes.

program test_primes

  use qmcpack

  integer(kind=i4b)                            :: i
  integer(kind=i4b), dimension(:), allocatable :: myprimes
  integer(kind=i4b), dimension(60), parameter  :: prime_ge_data = &
    (/  2,  2,  3,  5,  5,  7,  7, 11, 11, 11,                    &
       11, 13, 13, 17, 17, 17, 17, 19, 19, 23,                    &
       23, 23, 23, 29, 29, 29, 29, 29, 29, 31,                    &
       31, 37, 37, 37, 37, 37, 37, 41, 41, 41,                    &
       41, 43, 43, 47, 47, 47, 47, 53, 53, 53,                    &
       53, 53, 53, 59, 59, 59, 59, 59, 59, 61 /)


  call show_test_header("MOD_PRIMES")
  call init_assert()

  call start_test("Testing generation of first 1500 primes...")
  call assert(primes(1500), allprimes(1:1500))
  call stop_test()

  call start_test("Testing primes(1)...")
  call assert(primes(1), (/ 2 /))
  call stop_test()

  call start_test("Testing primes(5)...")
  call assert(primes(5), (/ 2, 3, 5, 7, 11 /))
  call stop_test()

  call start_test("Testing primes(1, skip=500)...")
  call assert(primes(1, skip=500), (/ 3581 /))
  call stop_test()

  call start_test("Testing primes(2, skip=500)...")
  call assert(primes(2, skip=500), (/3581, 3583/))
  call stop_test()

  call start_test("Testing primes(1, skip=1)...")
  call assert(primes(1, skip=1), (/ 3 /))
  call stop_test()

  call start_test("Testing primes(1, skip=2)...")
  call assert(primes(1, skip=2), (/ 5 /))
  call stop_test()

  call start_test("Checking the how-maniest prime 99991 is...")
  allocate(myprimes(9592))
  myprimes = primes(9592)
  write(unit=*, fmt="(A, I0)") "  prime 9592 = ", myprimes(9592)
  deallocate(myprimes)
  call stop_test()

  call start_test("Testing isprime(n) against the first 1500 primes...")
  allocate(myprimes(1500))
  myprimes = primes(1500)
  do i=1,1500
    if (.not. isprime(myprimes(i))) then
      print *, ""
      write(unit=*, fmt="(A, I0, A)") &
        "  isprime(n) said ", myprimes(i), " is not prime."
      call increase_nb_assert_errors(1)
    end if
  end do
  deallocate(myprimes)
  call stop_test()

  call start_test("Checking prime_ge(n) for n=1,60...")
  do i=1,size(prime_ge_data)
    call assert(prime_ge_data(i), prime_ge(i))
  end do
  call stop_test()

  call start_test("Testing nb_primes_up_to()...")
  call assert(nb_primes_up_to(2), 1)
  call assert(nb_primes_up_to(3), 2)
  call assert(nb_primes_up_to(4), 2)
  call assert(nb_primes_up_to(5), 3)
  call assert(nb_primes_up_to(6), 3)
  call assert(nb_primes_up_to(7), 4)
  call stop_test()

  call show_test_summary(get_nb_assert_errors()) 

end program test_primes
