program test_distance

  use qmcpack

  real(kind=qp), dimension(:), allocatable :: x, y
  real(kind=qp)                            :: d

  call show_test_header("MOD_DISTANCE")
  call init_assert()

  call start_test("Testing euclidian_distance(x,y)...")
  allocate(x(3), y(3))
  x = (/ 1.0_qp, 2.0_qp, 3.0_qp /)
  y = (/ 8.0_qp, 7.0_qp, 9.0_qp /)
  d = euclidian_distance(x, y)
  call assert(d, sqrt(110.0_qp))
  deallocate(x,y)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_distance
