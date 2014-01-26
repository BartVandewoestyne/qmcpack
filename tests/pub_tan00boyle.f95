program pub_tan00boyle

  use qmcpack

  integer(kind=i4b)                            :: n, s
  integer(kind=i4b)                            :: i
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b), dimension(:), allocatable :: step
  real(kind=qp), dimension(:), allocatable     :: x_1p
  character(len=20)                            :: fmt_string
  integer                                      :: my_unit

  call show_test_header("PUB TAN00BOYLE")
  call init_assert()

  ! Generate the (0,3,2)-net in base 2 with 2^3=8 elements from Fig 1.
  ! Note: the startindex is 64, but this is not mentioned in the paper!
  call start_test("Generating (0,3,2)-net in base 2 with 2^3=8 elements...")
  n = 8
  s = 2
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F18.15)"
  allocate(x_1p(s), step(s), startindex(s))
  step = 1
  startindex = 2**6
  call init_faure(n, s, init_scrambletype="None", init_startindex=startindex, init_step=step)
  call get_unit(my_unit)
  call checked_open(my_unit, "pub_tan00boyle_fig1.dat", "write")
  do i=1,n
    call next_faure(x_1p)
    write(unit=my_unit, fmt=fmt_string) x_1p
  end do
  call checked_close(my_unit)
  deallocate(x_1p, startindex, step)
  call free_faure()
  call stop_test()


  ! Generate the (0,3,2)-net in base 2 with 2^4=16 elements from Fig 2.
  ! Note: the startindex is 64, but this is not mentioned in the paper!
  call start_test("Generating (0,3,2)-net in base 2 with 2^4=16 elements...")
  n = 16
  s = 2
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F18.15)"
  allocate(x_1p(s), step(s), startindex(s))
  step = 1
  startindex = 2**6
  call init_faure(n, s, init_scrambletype="None", init_startindex=startindex, init_step=step)
  call get_unit(my_unit)
  call checked_open(my_unit, "pub_tan00boyle_fig2.dat", "write")
  do i=1,n
    call next_faure(x_1p)
    write(unit=my_unit, fmt=fmt_string) x_1p
  end do
  call checked_close(my_unit)
  deallocate(x_1p, startindex, step)
  call free_faure()
  call stop_test()


  ! Generate the (0,4,7)-net in base 7.
  ! Note: startindex and number of points are not mentioned in the paper,
  ! so we made a guess...
  call start_test("Generating (0,4,7)-net in base 7...")
  n = 2000
  s = 7
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F18.15)"
  allocate(x_1p(s), step(s), startindex(s))
  step = 1
  startindex = 0
  call init_faure(n, s, init_scrambletype="None", init_startindex=startindex, init_step=step)
  call get_unit(my_unit)
  call checked_open(my_unit, "pub_tan00boyle_fig3.dat", "write")
  do i=1,n
    call next_faure(x_1p)
    write(unit=my_unit, fmt=fmt_string) x_1p
  end do
  call checked_close(my_unit)
  deallocate(x_1p, step, startindex)
  call free_faure()
  call stop_test()

  ! Generate the (0,3,19)-net in base 19.
  ! Note: startindex and number of points are not mentioned in the paper,
  ! so we made a guess...
  call start_test("Generating (0,3,19)-net in base 19...")
  n = 4000
  s = 19
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F18.15)"
  allocate(x_1p(s), step(s), startindex(s))
  step = 1
  startindex = 0
  call init_faure(n, s, init_scrambletype="None", init_startindex=startindex, init_step=step)
  call get_unit(my_unit)
  call checked_open(my_unit, "pub_tan00boyle_fig5.dat", "write")
  do i=1,n
    call next_faure(x_1p)
    write(unit=my_unit, fmt=fmt_string) x_1p
  end do
  call checked_close(my_unit)
  deallocate(x_1p, step, startindex)
  call free_faure()
  call stop_test()


  call show_test_summary(get_nb_assert_errors())

end program pub_tan00boyle
