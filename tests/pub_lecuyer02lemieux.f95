program pub_lecuyer02lemieux

  use qmcpack

  integer(kind=i4b)                            :: n, s
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b), dimension(:), allocatable :: step
  integer(kind=i4b)                            :: i
  real(kind=qp), dimension(:), allocatable     :: x_1p
  integer                                      :: my_unit

  print *, "#################### REPRODUCING ARTICLE RESULTS #################"
  print *, "Article: Recent Advances in Randomized Quasi-Monte Carlo Methods"
  call init_assert()

  call start_test("Reproducing points from Figure 1.2 on page 12.")
  n = 81
  s = 2
  allocate(x_1p(s), step(s), startindex(s))
  step = 1
  startindex = 0
  call get_unit(my_unit)
  call checked_open(my_unit, "pub_lecuyer02lemieux.dat", "write")
  call init_faure(n, s, "None", init_base=3, init_startindex=startindex, init_step=step)
  do i=1,n
    call next_faure(x_1p)
    write(unit=my_unit, fmt="(2F18.15)") x_1p
  end do
  deallocate(x_1p, step , startindex)
  call free_faure()
  call stop_test()

end program pub_lecuyer02lemieux
