program test_periodize

  use qmcpack

  integer(kind=i4b)           :: nb_errors
  integer(kind=i4b)           :: N, i
  real(kind=qp), dimension(1) :: x, fi, si
  real(kind=qp)               :: h
  real(kind=qp)               :: y
  integer                     :: my_unit
  integer(kind=i4b)           :: showdim
  type(functionparams)        :: params

  call show_test_header("MOD_PERIODIZE")

  nb_errors = 0

  showdim = 1

  write(unit=*, fmt="(A)", advance="no") "Testing f_imt() by writing out data to file f_imt_1d.dat..."
  N = 1000
  h = 4.0_qp/N
  call get_unit(my_unit)
  call checked_open(my_unit, "f_imt_1d.dat", "write")
  do i = 1,N
    x = -2 + i*h
    call f_imt(x, fi, si)
    write(unit=my_unit, fmt="(3(es24.15e4))") x(showdim), fi(showdim), si(showdim)
  end do
  call checked_close(my_unit)
  write(unit=*, fmt="(A)") "done."


  write(unit=*, fmt="(A)", advance="no") "Testing f_mori() by writing out data to file f_mori_1d.dat..."
  N = 3000
  h = 1.0_qp/N
  call get_unit(my_unit)
  call checked_open(my_unit, "f_mori_1d.dat", "write")
  do i = 0,N
    x = 0 + i*h
    call f_mori(x, fi, si)
    write(unit=my_unit, fmt="(3(es24.15e4))") x(showdim), fi(showdim), si(showdim)
  end do
  call checked_close(my_unit)
  write(unit=*, fmt="(A)") "done."


  write(unit=*, fmt="(A)", advance="no") "Testing f_tanh() by writing out data to file f_tanh_1d.dat..."
  N = 1000
  h = 1.0_qp/N
  call get_unit(my_unit)
  call checked_open(my_unit, "f_tanh_1d.dat", "write")
  do i = 0,N
    x = 0 + i*h
    call f_tanh(x, fi, si)
    write(unit=my_unit, fmt="(3(es24.15e4))") x(showdim), fi(showdim), si(showdim)
  end do
  call checked_close(my_unit)
  write(unit=*, fmt="(A)") "done."


  write(unit=*, fmt="(A)", advance="no") "Testing f_de() by writing out data to file f_de_1d.dat..."
  N = 1000
  h = 1.0_qp/N
  call get_unit(my_unit)
  call checked_open(my_unit, "f_de_1d.dat", "write")
  do i = 0,N
    x = 0 + i*h
    call f_de(x, fi, si)
    write(unit=my_unit, fmt="(3(es24.15e4))") x(showdim), fi(showdim), si(showdim)
  end do
  call checked_close(my_unit)
  write(unit=*, fmt="(A)") "done."


  write(unit=*, fmt="(A)", advance="no") "Testing f_imt() by writing out periodized version of some function..."
  N = 1000
  h = 1.0_qp/N
  call get_unit(my_unit)
  call checked_open(my_unit, "f_periodized_1d.dat", "write")
  do i = 0,N
    x = 0 + i*h
    call periform(f_warnock02_1p, params, x, "IMT", y)
    write(unit=my_unit, fmt="(2(es24.15e4))") x(showdim), y
  end do
  call checked_close(my_unit)
  write(unit=*, fmt="(A)") "done."



  call show_test_summary(nb_errors)

end program test_periodize
