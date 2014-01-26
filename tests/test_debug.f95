program test_debug

  use qmcpack

  integer(kind=i4b), dimension(6, 6) :: y

  call show_test_header("MOD_DEBUG")
  call init_assert()

  y(1,:) = (/ 4,   0,  0,  0,  0,  0 /)
  y(2,:) = (/ 24, 16,  0,  0,  0,  0 /)
  y(3,:) = (/ 40, 24, 29,  0,  0,  0 /)
  y(4,:) = (/ 11, 40, 24, 23,  0,  0 /)
  y(5,:) = (/ 15, 11, 40, 24, 30,  0 /)
  y(6,:) = (/ 29, 15, 11, 40, 24, 28 /)
  call print_matrix(y)

  call show_test_summary(get_nb_assert_errors())

end program test_debug
