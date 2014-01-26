program test_owen_scrambling

  use qmcpack

  integer(kind=i4b)                              :: b
  integer(kind=i4b), dimension(:,:), allocatable :: x
  integer(kind=i4b), dimension(:), allocatable   :: shift

  call show_test_header("MOD_OWEN_SCRAMBLING")
  call init_assert()
  call random_seed()
  
  call start_test("Testing random_linear_scrambling(...)...")
  b = 3
  allocate(x(10,10), shift(10))
  call random_linear_scrambling(b, x, shift)
  call print_matrix(x)
  call print_vector(shift, "tranposed")
  deallocate(x, shift)
  call stop_test()


  call start_test("Testing random_linear_digit_scrambling(...)...")
  b = 3
  allocate(x(10,10), shift(10))
  call random_linear_digit_scrambling(b, x, shift)
  call print_matrix(x)
  call print_vector(shift, "tranposed")
  deallocate(x, shift)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_owen_scrambling
