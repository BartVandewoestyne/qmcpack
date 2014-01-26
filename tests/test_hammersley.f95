! References:
!

program test_hammersley

  use qmcpack

  integer(kind=i4b)                            :: n, s
  real(kind=qp), dimension(:,:), allocatable   :: x
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b), dimension(:), allocatable :: myprimes
  character(len=20)                            :: fmt_string
  integer(kind=i4b)                            :: i
  integer                                      :: my_unit

  call show_test_header("MOD_HAMMERSLEY")
  call init_assert()

  call start_test("Testing hammersley() with Shaw/Niederreiter's definition...")
  n = 10
  s = 3
  allocate(x(n,s), startindex(s-1), myprimes(s-1))
  startindex = 0
  myprimes = primes(s-1)
  call hammersley(n, s, myprimes, startindex, def="Shaw", x=x)
  call assert(x(10,:), (/ 0.9_qp, &
                          0.5625_qp, &
                          0.03703703703703703_qp /), &
                          1.0e-15_qp)
  deallocate(x, startindex, myprimes)
  call stop_test()


  call start_test("Testing hammersley() with Halton's definition...")
  n = 10
  s = 3
  allocate(x(n,s), startindex(s-1), myprimes(s-1))
  startindex = 1
  myprimes = primes(s-1)
  call hammersley(n, s, myprimes, startindex, def="Halton", x=x)
  call assert(x(10,:), (/ 1.0_qp, &
                          0.3125_qp, &
                          0.3703703703703703_qp /), &
                          1.0e-15_qp)
  deallocate(x, startindex, myprimes)
  call stop_test()


  call start_test("Testing hammersley() with Fang and Wang's definition...")
  n = 10
  s = 3
  allocate(x(n,s), startindex(s-1), myprimes(s-1))
  startindex = 1
  myprimes = primes(s-1)
  call hammersley(n, s, myprimes, startindex, def="Fang_Wang", x=x)
  call assert(x(10,:), (/ 0.95_qp, &
                          0.3125_qp, &
                          0.3703703703703703_qp /), &
                          1.0e-15_qp)
  deallocate(x, startindex, myprimes)
  call stop_test()


  call start_test("Writing some Hammersley points (FangWang) to file...")
  s = 2
  n = 1000000
  write(unit=fmt_string, fmt="(A)") "(F10.7, F18.15)"
  allocate(x(n,s), startindex(s-1), myprimes(s-1))
  startindex = 1
  myprimes = primes(s-1)
  call get_unit(my_unit)
  call checked_open(my_unit, "hammersley.dat", "write")
  call hammersley(n, s, myprimes, startindex, startindex, "Fang_Wang", x)
  do i=1,n
    write(unit=my_unit, fmt=fmt_string) x(i,:)
  end do 
  call checked_close(my_unit)
  deallocate(x, startindex, myprimes)
  call free_hammersley()
  call stop_test()


  call show_test_summary(get_nb_assert_errors())

end program test_hammersley
