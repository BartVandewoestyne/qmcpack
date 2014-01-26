program test_radical_inverse

  use qmcpack

  integer, dimension(:), allocatable           :: mydigits
  integer(kind=i4b), dimension(:), allocatable :: p
  integer(kind=i4b), dimension(:), allocatable :: myprimes
  integer(kind=i4b)                            :: i
  integer(kind=i4b), dimension(31), parameter  :: warnock_spinvalues = &
    (/  1,  2,  2,  5,  3,  8,  3,  7, 18, 12, 18,  4, 18, 24, 40, &
       14, 41, 50, 12, 30, 40, 70, 10, 39, 82,  6, 16, 37, 48, 71, &
       34 /)
  
  call show_test_header("MOD_RADICAL_INVERSE")
  call init_assert()


  call start_test("Testing get_nb_digits() and get_digits()...")

  call assert(get_nb_digits(154754, 10), 6)
  allocate(mydigits(get_nb_digits(154754, 10)))
  call get_digits(154754, 10, mydigits)
  call assert(mydigits, (/ 4,5,7,4,5,1 /))
  deallocate(mydigits)

  call assert(get_nb_digits(2, 10), 1)
  allocate(mydigits(get_nb_digits(2, 10)))
  call get_digits(2, 10, mydigits)
  call assert(mydigits, (/ 2 /))
  deallocate(mydigits)

  call assert(get_nb_digits(0, 2), 1)
  allocate(mydigits(get_nb_digits(0, 2)))
  call get_digits(0, 2, mydigits)
  call assert(mydigits, (/ 0 /))
  deallocate(mydigits)

  call assert(get_nb_digits(0, 3), 1)
  allocate(mydigits(get_nb_digits(0, 3)))
  call get_digits(0, 3, mydigits)
  call assert(mydigits, (/ 0 /))
  deallocate(mydigits)

  call stop_test()


  call start_test("Testing radical_inverse()...")
  call assert(radical_inverse(2, 11), 2.0_qp/11)
  call assert(radical_inverse(154754, 10), 0.457451_qp)
  call assert(radical_inverse(2, 10), 0.2_qp)
  call assert(radical_inverse(0, 10), 0.0_qp)
  call assert(radical_inverse(1037, 2), 0.68798828125_qp)
  call assert(radical_inverse(586, 5), 194.0_qp/625)
  call assert(radical_inverse((/ 1, 2 /), 2), (/ 0.5_qp, 0.25_qp /))
  call stop_test()


  call start_test("Testing radical_inverse_scrambled()...")
  allocate(p(0:10))
  p = (/ 0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 10 /)
  call assert(radical_inverse_scrambled(2, 11, p), p(2)/11.0_qp)
  call assert(radical_inverse_scrambled(3, 11, p), p(3)/11.0_qp)
  deallocate(p)
  call stop_test()


  call start_test("Testing radical_inverse_lin_scrambled()...")
  call assert(radical_inverse_lin_scrambled(2, 11, 3,  4), 10/11.0_qp)
  call assert(radical_inverse_lin_scrambled(15, 10, 11, 3), 0.84_qp)
  call stop_test()


  call start_test("Testing warnock_spin_get()...")
  allocate(myprimes(31))
  myprimes = primes(31)    
  do i=1,size(myprimes)
    call assert(warnock_spin_get(myprimes(i)), &
                    warnock_spinvalues(i))
  end do
  deallocate(myprimes)
  call stop_test()


  call start_test("Testing radical_inverse_folded()...")
  call assert(radical_inverse_folded(10, 2), 1.0_qp/48)
  call assert(radical_inverse_folded(11, 3), 551.0_qp/702)
  call stop_test()

  
  call start_test("Testing radical_inverse_inversive()...")
  call assert(radical_inverse_inversive(1, 2, 5, 3), 0.0_qp)
  call stop_test()


  call start_test("Testing radical_inverse_quadratic()...")
  call assert(radical_inverse_quadratic(1, 2, 5, 3, 2), 0.0_qp)
  call assert(radical_inverse_quadratic(1, 2, 5, 3, 1), 0.5_qp)
  call assert(radical_inverse_quadratic(586, 5, 7, 3, 2), 296.0_qp/625)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_radical_inverse
