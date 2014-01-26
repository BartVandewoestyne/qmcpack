program test_pdf

  use qmcpack

  real(kind=qp)     :: x
  integer(kind=i4b) :: i
  integer(kind=i4b) :: my_unit
  integer           :: ios

  call show_test_header("MOD_PDF")
  call init_assert()

  call random_seed()

  call start_test("Testing randn_box_muller(x) by writing 10000 samples " // &
    "to file randn_box_muller.dat...")

  call get_unit(my_unit)
  open(unit=my_unit, file="randn_box_muller.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  if (ios /= 0) then
    write(unit=*, fmt="(A)") "ERROR while opening randn_box_muller.dat!"
    call increase_nb_assert_errors(1)
  else
    do i=1,10000
      call randn_box_muller(x, 5.0_qp, 4.0_qp)
      write(unit=my_unit, fmt="(ES25.15)") x
    end do
  end if

  close(unit=my_unit)
  call stop_test()


  call start_test("Testing norm_pdf()...")
  call assert(norm_pdf(4.0_qp, 0.0_qp, 1.0_qp), 0.0001338302257648854_qp)
  call assert(norm_pdf(0.0_qp, 1.0_qp, 4.0_qp), 0.1760326633821498_qp)
  call stop_test()


  call start_test("Testing norm_pdf(...) by writing samples to file norm_pdf.dat...")
  call get_unit(my_unit)
  open(unit=my_unit, file="norm_pdf.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  if (ios /= 0) then
    write(unit=*, fmt="(A)") "ERROR while opening norm_pdf.dat!"
    call increase_nb_assert_errors(1)
  else
    do i=-1,100
      write(unit=my_unit, fmt="(ES25.15)") norm_pdf(0.1_qp*i, 5.0_qp, 4.0_qp)
    end do
  end if

  close(unit=my_unit)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_pdf
