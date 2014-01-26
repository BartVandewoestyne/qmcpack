! TODO: add the graycode(n) for the startindex when choosing the
! Antonov-Saleev variant
program write_pointset

  use qmcpack

  integer(kind=i4b)                             :: n, s
  real(kind=qp), dimension(:), allocatable      :: x_1d
  integer(kind=i4b)                             :: nb_errors
  type(soboltype)                               :: mysoboltype
  integer(kind=i4b), dimension(:), allocatable  :: startindex
  integer(kind=i4b)                             :: start_n
  integer(kind=i4b)                             :: i
  integer(kind=i4b)                             :: my_unit
  integer                                       :: ios

  write(unit=*, fmt="(A)") "#################### WRITING SEQUENCE #################"
  nb_errors = 0

  write(unit=*, fmt="(A)", advance="no") "Enter start_n: "
  read *, start_n

  write(unit=*, fmt="(A)") "Testing init_sobol() and next_sobol() by writing to a data file"
  mysoboltype%poly_order = "JoeKuo"
  mysoboltype%dirnumbers = "JoeKuo"
  mysoboltype%use_ones = .true.
  mysoboltype%use_antonov_saleev = .true.
  n = 6
  s = 100
  allocate(startindex(s), x_1d(s))
  startindex = start_n
  call init_sobol(n, s, startindex, mysoboltype)
  call get_unit(my_unit)
  open(unit=my_unit, file="sobol_compare_to_george.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  if (ios == 0) then
    do i=1,n
      write(unit=my_unit, fmt="(A12, I5)") "Sobol point ", i+start_n-1
      call next_sobol(x_1d)
      write(unit=my_unit, fmt="(10F9.6)") x_1d
    end do
    close(unit=my_unit)
  else
    write(unit=*, fmt="(A)") "ERROR: could not write sobol dataset to file."
    nb_errors = nb_errors + 1
  end if
  deallocate(startindex, x_1d)
  write(unit=*, fmt="(A)") "done."


  if (nb_errors > 0) then
    write(unit=*, fmt="(I0.0, A)") nb_errors, " errors found in mod_sobol!!!"
  else
    write(unit=*, fmt="(A)") "=> mod_sobol is ok."
  end if


end program write_pointset
