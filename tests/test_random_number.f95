program test_random_number

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

  call show_test_header("RANDOM NUMBERS")
  call init_assert()


  call start_test("Writing random numbers to random_number_np_sd.dat...")
  n = 1024
  s = 2
  allocate(x2(n,s))
  call random_number(x2)
  call write_sequence(x2, "pseudorandom", n, s)
  deallocate(x2)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_random_number
