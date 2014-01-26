! References:
!
!   [1] Press, William H. and Teukolsky, Saul A., `Quasi- (that is, Sub-)
!       Random Numbers', Computers in Physics, Vol. 3, nov/dec 1989,
!       pages 76-79.
!
! Notes:
!
!   * Another test for correctness could be to compare a plot of the points
!     with the plots on page 77 of [1].

program test_sobol

  use qmcpack

  integer(kind=i4b)                            :: n, s
  real(kind=qp), dimension(:), allocatable     :: x_1p
  real(kind=qp), dimension(:,:), allocatable   :: x_np
  type(soboltype)                              :: mysoboltype
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b)                            :: i
  character(len=10)                            :: fmt_string
  integer                                      :: my_unit

  call show_test_header("MOD_SOBOL")
  call init_assert()

  call start_test("Testing sobol(n, s, ...) by writing numbers to a datafile, using")
  n = 10
  s = 4
  ! Set the type of Sobol sequence we want
  mysoboltype%poly_order = "JoeKuo"
  mysoboltype%dirnumbers = "JoeKuo"
  mysoboltype%use_ones = .true.
  mysoboltype%use_antonov_saleev = .true.
  call print_sobol_type(mysoboltype)
  allocate(x_np(n,s), startindex(s))
  startindex = 1
  call sobol(n, s, startindex, x_np, mysoboltype)
  call write_sequence(x_np, "sobol", n, s)
  deallocate(x_np, startindex)
  write(unit=*, fmt="(A)", advance="no") "... "
  call stop_test()


  ! Check for the numbers on page 792 in [1].
  call start_test("Testing sobol(n, s, ...) against Sobol numbers from a paper by Sobol...")
  n = 9
  s = 4
  allocate(x_np(n,s), startindex(s))
  mysoboltype%poly_order = "BratleyFox"
  mysoboltype%dirnumbers = "BratleyFox"
  mysoboltype%use_ones = .true.
  mysoboltype%use_antonov_saleev = .false.
  startindex = 0
  call sobol(n, s, startindex, x_np, mysoboltype)
  call assert(x_np(1,:), &
        (/ 0.0_qp, 0.0_qp, 0.0_qp, 0.0_qp /))
  call assert(x_np(2,:), &
        (/ 0.5_qp, 0.5_qp, 0.5_qp, 0.5_qp /))
  call assert(x_np(3,:), &
        (/ 0.25_qp, 0.75_qp, 0.25_qp, 0.75_qp /))
  call assert(x_np(4,:), &
        (/ 0.75_qp, 0.25_qp, 0.75_qp, 0.25_qp /))
  call assert(x_np(5,:), &
        (/ 0.125_qp, 0.625_qp, 0.875_qp, 0.875_qp /))
  call assert(x_np(6,:), &
        (/ 0.625_qp, 0.125_qp, 0.375_qp, 0.375_qp /))
  call assert(x_np(7,:), &
        (/ 0.375_qp, 0.375_qp, 0.625_qp, 0.125_qp /))
  call assert(x_np(8,:), &
        (/ 0.875_qp, 0.875_qp, 0.125_qp, 0.625_qp /))
  call assert(x_np(9,:), &
        (/ 0.0625_qp, 0.9375_qp, 0.6875_qp, 0.3125_qp /))
  deallocate(x_np, startindex)
  call stop_test()


  call start_test("Testing sobol(n, s, ...) for the Numerical Recipes version...")
  n = 1000
  s = 6
  allocate(x_np(n,s), startindex(s))
  mysoboltype%poly_order = "NumericalRecipes"
  mysoboltype%dirnumbers = "NumericalRecipes"
  mysoboltype%use_ones = .false.
  mysoboltype%use_antonov_saleev = .true.
  startindex = 0
  call sobol(n, s, startindex, x_np, mysoboltype)
  if (any(x_np > 1) .or. any(x_np < 0) ) then
    write(unit=*, fmt="(A)")
    write(unit=*, fmt="(A)") &
      "ERROR: found sobol point which is larger than one!!!"
  end if
  deallocate(x_np, startindex)
  call stop_test()


  call start_test("Testing sobol(n, s, ...) against some Sobol numbers from " // &
    "Bratley and Fox's paper... ")
  ! Check for the numbers on page 91-92 in [2].
  ! Since in Bratley and Fox their implementation the primpoly for the
  ! 4th dimension is x^3+x+1 = 11, we can check against this fourth dimension to
  ! see if we get the values from the paper.
  !
  n = 24
  s = 4
  allocate(x_np(n,s), startindex(s))
  mysoboltype%poly_order = "BratleyFox"
  mysoboltype%dirnumbers = "BratleyFox"
  mysoboltype%use_ones = .true.
  mysoboltype%use_antonov_saleev = .true.
  startindex = 0
  call sobol(n, s, startindex, x_np, mysoboltype)
  call assert(x_np(1,4), 0.0_qp)
  call assert(x_np(2,4), 0.5_qp)
  call assert(x_np(3,4), 0.25_qp)
  call assert(x_np(4,4), 0.75_qp)
  call assert(x_np(24,4), 0.53125_qp)
  deallocate(x_np, startindex)
  call stop_test()


  call start_test("Testing init_sobol() and next_sobol() by writing points to " // &
    "sobol_joe_kuo.dat, using")
  mysoboltype%poly_order = "JoeKuo"
  mysoboltype%dirnumbers = "JoeKuo"
  mysoboltype%use_ones = .true.
  mysoboltype%use_antonov_saleev = .true.
  call print_sobol_type(mysoboltype)
  n = 1000
  s = 40
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F19.16)"
  allocate(startindex(s), x_1p(s))
  startindex = 1
  call init_sobol(n, s, startindex, mysoboltype)
  call get_unit(my_unit)
  call checked_open(my_unit, "sobol_joe_kuo.dat", "write")
  do i=1,n
    call next_sobol(x_1p)
    write(unit=my_unit, fmt=fmt_string) x_1p
  end do
  call checked_close(my_unit)
  deallocate(startindex, x_1p)
  write(unit=*, fmt="(A)", advance="no") "... "
  call stop_test()


  call start_test("Testing init_sobol() and next_sobol() for starting at a random starting " // &
    "point... ")
  mysoboltype%poly_order = "JoeKuo"
  mysoboltype%dirnumbers = "JoeKuo"
  mysoboltype%use_ones = .true.
  mysoboltype%use_antonov_saleev = .true.
  n = 10
  s = 100
  allocate(startindex(s), x_1p(s))
  startindex = 1000
  call init_sobol(n, s, startindex, mysoboltype)
  call next_sobol(x_1p)
  call assert(x_1p(1), 0.2197265625_qp)
  call next_sobol(x_1p)
  call assert(x_1p(1), 0.7197265625_qp)
  call next_sobol(x_1p)
  call assert(x_1p(1), 0.9697265625_qp)
  call next_sobol(x_1p)
  call assert(x_1p(1), 0.4697265625_qp)
  deallocate(startindex, x_1p)
  call stop_test()


  ! We use the 4th dimension, because this is the dimension that has
  ! primpoly 11.
  call start_test("Testing init_sobol() and " // &
    "next_sobol() against some Sobol numbers from Bratley and Fox's paper... ")
  mysoboltype%poly_order = "BratleyFox"
  mysoboltype%dirnumbers = "BratleyFox"
  mysoboltype%use_ones = .true.
  mysoboltype%use_antonov_saleev = .true.
  n = 24
  s = 4 
  allocate(startindex(s), x_1p(s))
  startindex = 1
  call init_sobol(n, s, startindex, mysoboltype)
  call next_sobol(x_1p)
  call assert(x_1p(4), 0.5_qp)
  call next_sobol(x_1p)
  call assert(x_1p(4), 0.25_qp)
  call next_sobol(x_1p)
  call assert(x_1p(4), 0.75_qp)
  do i=1,20
    call next_sobol(x_1p)
  end do
  call assert(x_1p(4), 0.53125_qp)
  deallocate(startindex, x_1p)
  call stop_test()


!  write(unit=*, fmt="(A)", advance="no") &
!    "Testing init_sobol() and next_sobol() by generating 2^22 points &
!    &in 1000 dimensions... "
!  mysoboltype%poly_order = "JoeKuo"
!  mysoboltype%dirnumbers = "JoeKuo"
!  mysoboltype%use_ones = .true.
!  mysoboltype%use_antonov_saleev = .true.
!  ! 2**22 = 4194304
!  n = 2**22
!  s = 1000
!  allocate(startindex(s), x_1p(s))
!  startindex = 1
!  call init_sobol(n, s, startindex, mysoboltype)
!  do i=1,n
!    call next_sobol(x_1p)
!  end do
!  deallocate(startindex, x_1p)
!  write(unit=*, fmt="(A)") "done."


  call start_test("Testing init_sobol() and next_sobol() with NumericalRecipes " // &
    "and writing to sobol_numerical_recipes.dat...")
  mysoboltype%poly_order = "NumericalRecipes"
  mysoboltype%dirnumbers = "NumericalRecipes"
  mysoboltype%use_ones = .false.
  mysoboltype%use_antonov_saleev = .true.
  n = 100
  s = 6
  allocate(startindex(s), x_1p(s))
  startindex = 1
  call init_sobol(n, s, startindex, mysoboltype)
  call get_unit(my_unit)
  call checked_open(my_unit, "sobol_numerical_recipes.dat", "write")
  do i=1,n
    call next_sobol(x_1p)
    write(unit=my_unit, fmt="(6F19.16)") x_1p
    if (any(x_1p > 1) .or. any(x_1p < 0) ) then
      write(unit=*, fmt="(A)")
      write(unit=*, fmt="(A)") &
        "ERROR: found sobol point which is larger than one!!!"
    end if
  end do
  deallocate(startindex, x_1p)
  call stop_test()


  call start_test("Testing init_sobol() and next_sobol() on 1-dimensional sequences " // &
    "in 1000 dimensions... ")
  mysoboltype%poly_order = "JoeKuo"
  mysoboltype%dirnumbers = "JoeKuo"
  mysoboltype%use_ones = .true.
  mysoboltype%use_antonov_saleev = .true.
  n = 1000
  s = 1
  allocate(startindex(s), x_1p(s))
  startindex = 1
  call init_sobol(n, s, startindex, mysoboltype)
  do i=1,n
    call next_sobol(x_1p)
  end do
  deallocate(startindex, x_1p)
  call stop_test()

  call free_sobol()

  call show_test_summary(get_nb_assert_errors())

end program test_sobol
