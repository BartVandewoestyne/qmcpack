! References:
!
!   [1] 'Number-theoretic Methods in Statistics', Wang, Y. and Fang, K.-T.,
!       1994, page 27.

program test_weyl

  use qmcpack

  integer(kind=i4b)                            :: n, s
  real(kind=qp), dimension(:,:), allocatable   :: x
  real(kind=qp), dimension(:), allocatable     :: x2
  integer(kind=i4b), dimension(:), allocatable :: startindex, step
  integer(kind=i4b)                            :: i, j
  integer(kind=i4b), dimension(:), allocatable :: myprimes
  integer                                      :: ios

  call show_test_header("MOD_WEYL")
  call init_assert()

  call start_test("Testing block assignment with square_root_seq()...")
  n = 1000
  s = 3
  allocate(x(n,s), startindex(s), myprimes(s))
  startindex = 1
  myprimes = primes(s)
  call square_root_seq(n, s, startindex, myprimes, x)
  call assert(x(1,:), (/ 0.4142135623730951_qp,    &
                         0.7320508075688772_qp,    &
                         0.2360679774997898_qp /))
  call assert(x(2,:), (/ 0.8284271247461903_qp,    &
                         0.4641016151377544_qp,    &
                         0.4721359549995796_qp /))
  call assert(x(3,:), (/ 0.2426406871192857_qp,    &
                         0.1961524227066320_qp,    &
                         0.7082039324993694_qp /))
  call assert(x(4,:), (/ 0.6568542494923806_qp,    &
                         0.9282032302755088_qp,    &
                         0.9442719099991592_qp /))
  call write_sequence(x, "weyl", n, s)
  deallocate(x, startindex, myprimes)
  call stop_test()


  call start_test("Testing next_square_root() function...")
  allocate(x2(3), myprimes(3), startindex(3))
  startindex = 1
  myprimes = primes(3)
  call init_square_root_seq(3, myprimes, startindex)
  call next_square_root(x2)
  call assert(x2, (/ 0.4142135623730951_qp,    &
                     0.7320508075688772_qp,    &
                     0.2360679774997898_qp /))

  call next_square_root(x2)
  call assert(x2, (/ 0.8284271247461903_qp,    &
                     0.4641016151377544_qp,    &
                     0.4721359549995796_qp /))
  call next_square_root(x2)
  call assert(x2, (/ 0.2426406871192857_qp,    &
                     0.1961524227066320_qp,    &
                     0.7082039324993694_qp /))
  call next_square_root(x2)
  call assert(x2, (/ 0.6568542494923806_qp,    &
                     0.9282032302755088_qp,    &
                     0.9442719099991592_qp /))
  deallocate(x2, myprimes, startindex)
  call stop_test()

  call start_test("Testing next_gp_set_b() function...")
  call init_gp_set_b(s=5, p=2)
  allocate(x2(5))
  do i=1,10
    call next_gp_set_b(x2)
  end do
  call assert(x2, (/ 0.2246204830937302_qp, &
                     0.5992104989487324_qp, &
                     0.1421356237309510_qp, &
                     0.8740105196819936_qp, &
                     0.8179743628067868_qp /))
  deallocate(x2)
  call free_weyl()
  call stop_test()


  call start_test("Testing gp_set_b()...")
  allocate(x(10, 5))
  call gp_set_b(n=10, s=5, p=2, x=x)
  call assert(x(10,:), (/ 0.2246204830937302_qp, &
                          0.5992104989487324_qp, &
                          0.1421356237309510_qp, &
                          0.8740105196819936_qp, &
                          0.8179743628067868_qp /)) 
  call stop_test()


  ! Write out a file in Takhtamyshev format
  ! to see if we use the same points
  call start_test("Writing out some points to richtmyer_takhtamyshev.dat...")
  n = 20
  s = 100
  allocate(x2(s), startindex(s), step(s), myprimes(s))
  myprimes = primes(s, 501)
  startindex = 100002
  step = 50001
  call init_square_root_seq(s, myprimes, startindex, step)
  open(unit=unit_number, file="richtmyer_takhtamyshev.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  do i=1,n
    write(unit=unit_number, fmt="(A18, I2)") "Square-root point ", i
    call next_square_root(x2)
    do j=1,10
      write(unit=unit_number, fmt="(10F20.15)") x2(10*(j-1)+1:10*j)
    end do
  end do
  deallocate(x2, startindex, step, myprimes)
  close(unit=unit_number)
  call free_weyl()
  call stop_test()


  call show_test_summary(get_nb_assert_errors())

end program test_weyl
