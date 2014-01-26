program pub_jaeckel

  use qmcpack

  integer(kind=i4b)                            :: n, s
  real(kind=qp), dimension(:), allocatable     :: x_1p
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b)                            :: i

  write(unit=*, fmt="(A)") "#################### REPRODUCING PUBLISHED RESULTS ################"
  write(unit=*, fmt="(A)") "Book: Monte Carlo methods in finance by Peter Jaeckel."
  call init_assert()

  call start_test("Reproducing Halton numbers from page 80...")
  n = 37
  s = 4
  allocate(x_1p(s), startindex(s))
  startindex = 1
  call init_halton(s, startindex)
  do i=1, n
    call next_halton(x_1p)
  end do
  call assert(x_1p, (/0.640625_qp, 0.382716_qp,   &
                      0.488000_qp, 0.387755_qp/), &
                      0.000001_qp)
  deallocate(x_1p, startindex)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program pub_jaeckel
