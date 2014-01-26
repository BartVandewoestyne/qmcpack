program test_halton

  use qmcpack

  integer(kind=i4b)                            :: n, s
  integer(kind=i4b)                            :: i
  real(kind=qp), dimension(:,:), allocatable   :: x2
  real(kind=qp), dimension(:), allocatable     :: x
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b), dimension(:), allocatable :: u, v, w
  character(len=20)                            :: fmt_string
  integer(kind=i4b), dimension(:), allocatable :: myprimes
  integer                                      :: my_unit

  call show_test_header("MOD_HALTON")
  call init_assert()


  call start_test("Testing some Halton numbers against their theoretical value...")
  n = 1000
  s = 50
  allocate(x2(n,s), startindex(s), myprimes(s))
  startindex = 1
  myprimes = primes(s)
  call halton(n, s, startindex, myprimes, x2)
  call assert(x2(998,1), 0.4052734375_qp)
  call assert(x2(998,50), 0.3581548788161935_qp)
  call assert(x2(1000,5), 0.9316303531179565_qp)
  deallocate(x2, startindex, myprimes)
  call stop_test()


  call start_test("Testing against Halton numbers from 'Statistics and Computing' by James E. Gentle...")
  n = 5
  s = 3
  allocate(x2(n,s), startindex(s), myprimes(s))
  startindex = 15
  myprimes = primes(s)
  call halton(n, s, startindex, myprimes, x2)
  call assert(x2(1,:), (/ 0.937500_qp,    &
                          0.259259_qp,    &
                          0.120000_qp /))
  call assert(x2(2,:), (/ 0.031250_qp,    &
                          0.592593_qp,    &
                          0.320000_qp /))
  call assert(x2(3,:), (/ 0.531250_qp,    &
                          0.925926_qp,    &
                          0.520000_qp /))
  call assert(x2(4,:), (/ 0.281250_qp,    &
                          0.074074_qp,    & ! This one was wrong in the book!
                          0.720000_qp /))
  call assert(x2(5,:), (/ 0.781250_qp,    &
                          0.407407_qp,    & ! This one was wrong in the book
                          0.920000_qp /))
  deallocate(x2, startindex, myprimes)
  call stop_test()


  call start_test("Writing a Halton sequence to file halton_np_sd.dat...")
  n = 1024
  s = 2
  allocate(x2(n,s), startindex(s), myprimes(s))
  startindex = 1
  myprimes = primes(s)
  call halton(n, s, startindex, myprimes, x2)
  call write_sequence(x2, "halton", n, s)
  !call random_number(x2)
  !call write_sequence(x2, "pseudorandom", n, s)
  deallocate(x2, startindex, myprimes)
  call stop_test()


  call start_test("Testing halton()...")
  allocate(x2(4,3), startindex(3), myprimes(3))
  startindex = 0
  myprimes = primes(3)
  call halton(4, 3, startindex, myprimes, x2)
  call assert(x2(1,:), (/ 0.0_qp, 0.0_qp, 0.0_qp /))
  call assert(x2(2,:), (/ 0.5_qp, 1.0_qp/3, 1.0_qp/5 /))
  call assert(x2(3,:), (/ 0.25_qp, 2.0_qp/3, 2.0_qp/5 /))
  call assert(x2(4,:), (/ 0.75_qp, 1.0_qp/9, 3.0_qp/5 /))
  deallocate(x2, startindex, myprimes)
  call stop_test()


  call start_test("Testing init_halton() and next_halton()...")
  allocate(x(3), startindex(3), myprimes(3))
  startindex = 0
  myprimes = primes(3)
  call init_halton(3, startindex, myprimes)
  call next_halton(x)
  call assert(x, (/ 0.0_qp, 0.0_qp, 0.0_qp /))
  call next_halton(x)
  call assert(x, (/ 0.5_qp, 1.0_qp/3, 1.0_qp/5 /))
  call next_halton(x)
  call assert(x, (/ 0.25_qp, 2.0_qp/3, 2.0_qp/5 /))
  call next_halton(x)
  call assert(x, (/ 0.75_qp, 1.0_qp/9, 3.0_qp/5 /))
  call free_halton()
  deallocate(x, startindex, myprimes)
  call stop_test()


  call start_test("Testing if next_halton() generates points for which abs(x)>1...")
  allocate(x(3), startindex(3), myprimes(3))
  startindex = 0
  myprimes = primes(3)
  call init_halton(3, startindex, myprimes)
  do i=1,1000
    call next_halton(x)
    if (maxval(abs(x)) > 1) then
      write(unit=*, fmt="(A)")
      write(unit=*, fmt="(A)") "  ERROR: found wrong Halton number > 1"
      call increase_nb_assert_errors(1)
    end if
  end do
  call free_halton()
  deallocate(x, startindex, myprimes)
  call stop_test()


  !print *, "Testing if we stop on overflow..."
  !allocate(x(2))
  !call init_halton(2, (/ huge(i4b)-1, huge(i4b)-1 /), primes(2))
  !call next_halton(x)
  !write(unit=*, fmt="(A)", advance="no") x
  !call next_halton(x)
  !write(unit=*, fmt="(A)", advance="no") x
  !deallocate(x)
  !write(unit=*, fmt="(A)", advance="no") "done."


  call start_test("Testing init_halton() and next_halton_chi_permuted()...")
  allocate(x(40))
  allocate(startindex(40), myprimes(40))
  myprimes = primes(40)
  startindex = 0
  call init_halton(40, startindex, myprimes)
  call next_halton_chi_permuted(x)
  call assert(x(1),  0.0_qp)
  call assert(x(40), 0.0_qp)
  call next_halton_chi_permuted(x)
  call assert(x(1), 0.5_qp)
  call assert(x(40), 66.0_qp/173)
  call next_halton_chi_permuted(x)
  call assert(x(20), 17.0_qp/71)
  deallocate(x, startindex, myprimes)
  call free_halton()
  call stop_test()


  call start_test("Testing next_halton_modified() against some numbers generated"//&
    "with the program from Atanassov's site...")
  allocate(x(10), startindex(10), myprimes(10))
  startindex = 1
  myprimes = primes(10)
  call init_halton(10, startindex, myprimes)
  do i=1,100000
    call next_halton_modified(x)
  end do
  call assert(x(1),  0.021019_qp)
  call assert(x(2),  0.424822_qp)
  call assert(x(3),  0.00023296_qp)
  call assert(x(4),  0.928287_qp)
  call assert(x(5),  0.865614_qp)
  call assert(x(6),  0.445804_qp)
  call assert(x(7),  0.942873_qp)
  call assert(x(8),  0.10689_qp)
  call assert(x(9),  0.740303_qp)
  call assert(x(10), 0.665782_qp)
  deallocate(x, startindex, myprimes)
  call free_halton()
  call stop_test()


  call start_test("Testing next_halton_phicf()...")
  s = 31
  allocate(x(s), startindex(s), myprimes(s))
  startindex = 1
  myprimes = primes(s)
  call init_halton(s, startindex, myprimes)
  do i=1,12
    call next_halton_phicf(x)
  end do
  deallocate(x, startindex, myprimes)
  call free_halton()
  call stop_test()


  call start_test("Testing next_halton_inversive_scrambled()...")
  s = 5
  n = 1000
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F18.15)"
  allocate(x(s), startindex(s), myprimes(s))
  startindex = 1
  myprimes = primes(s)
  call get_unit(my_unit)
  call checked_open(my_unit, "halton_inversive_scrambled.dat", "write")
  call init_halton(s, startindex, myprimes)
  do i=1,n
    call next_halton_inversive_scrambled(x)
    write(unit=my_unit, fmt=fmt_string) x
  end do
  call checked_close(my_unit)
  deallocate(x, startindex, myprimes)
  call free_halton()
  call stop_test()


  call start_test("Testing next_halton_quadratic_scrambled()...")
  s = 20
  n = 1000
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F18.15)"
  allocate(x(s), startindex(s), myprimes(s), u(s), v(s), w(s))
  myprimes = primes(s)
  startindex = 1
  call get_unit(my_unit)
  call checked_open(my_unit, "halton_quadratic_scrambled.dat", "write")
  call init_halton(s, startindex, myprimes)
  call random_integer(myprimes(1)-1, u(1))
  v(1) = 1
  call random_integer(myprimes(1)-1, w(1))
  do i = 2,s
    call random_integer(myprimes(i)-1, u(i))
    call primitive_root_prime(myprimes(i), v(i))
    call random_integer(myprimes(i)-1, w(i))
  end do
  call set_quadratic_scramble(u, v, w)
  write(unit=*, fmt="(A, 20I3)") "p = ", myprimes
  write(unit=*, fmt="(A, 20I3)") "u = ", u
  write(unit=*, fmt="(A, 20I3)") "v = ", v
  write(unit=*, fmt="(A, 20I3)") "w = ", w
  do i=1,n
    call next_halton_quadratic_scrambled(x)
    write(unit=my_unit, fmt=fmt_string) x
  end do
  call checked_close(my_unit)
  deallocate(x, startindex, myprimes, u, v, w)
  call free_halton()
  call stop_test()


  call start_test("Testing next_halton_scrambled using Faure scrambling...")
  allocate(x(5), startindex(5), myprimes(5))
  startindex = 1
  myprimes = primes(5)
  call init_halton(5, startindex, myprimes, "Faure")
  do i=1,10
    call next_halton_scrambled(x)
  end do
  call assert(x, (/ 0.3125000000000000_qp, &
                    0.3703703703703703_qp, &
                    2.0_qp/25.0_qp, &
                    0.4693877551020408_qp, &
                    0.9090909090909092_qp /), &
                    2.0E-16_qp)
  deallocate(x, startindex, myprimes)
  call free_halton()
  call stop_test()


  call start_test("Testing next_halton_scrambled using Braaten and Weller scrambling...")
  allocate(x(5), startindex(5), myprimes(5))
  startindex = 1
  myprimes = primes(5)
  call init_halton(5, startindex, myprimes, "Braaten_Weller")
  do i=1,10
    call next_halton_scrambled(x)
  end do
  call assert(x, (/ 0.3125000000000000_qp, &
                    0.7407407407407407_qp, &
                    0.0400000000000000_qp, &
                    0.9387755102040816_qp, &
                    0.3636363636363636_qp /), &
                    0.0000000000000002_qp)
  deallocate(x, startindex, myprimes)
  call free_halton()
  call stop_test()


  call start_test("Testing next_halton_scrambled using Tuffin's MCL scrambling...")
  allocate(x(5), startindex(5), myprimes(5))
  startindex = 1
  myprimes = primes(5)
  call init_halton(5, startindex, myprimes, "MCL")
  do i=1,100
    call next_halton_scrambled(x)
  end do
  call assert(x, (/ 0.1484375000000000_qp,    &
                    0.4115226337448559_qp,    &
                    0.0240000000000000_qp,    &
                    0.2915451895043731_qp,    &
                    0.1239669421487603_qp /))
  deallocate(x, startindex, myprimes)
  call free_halton()
  call stop_test()


  call start_test("Testing random_startindex()...")
  allocate(startindex(5))
  call random_startindex(startindex)
  write(unit=*, fmt="(A)") "  Random startindices with maximum capacity:"
  write(unit=*, fmt="(A, 5I15)") "    ", startindex
  call random_startindex(startindex, 100)
  write(unit=*, fmt="(A)") "  Random startindices between 0 and 100:"
  write(unit=*, fmt="(A, 5I4)") "    ", startindex
  deallocate(startindex)
  call stop_test()


  call show_test_summary(get_nb_assert_errors())

end program test_halton
