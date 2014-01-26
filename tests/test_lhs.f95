program test_lhs

  use qmcpack

!  implicit none

  integer(kind=i4b)                          :: n, s
  real(kind=qp), dimension(:,:), allocatable :: x
  integer(kind=i4b)                          :: i
  integer(kind=i4b)                          :: nb_errors

  call show_test_header("MOD_LHS")

  nb_errors = 0

  n = 10
  s = 3

  allocate(x(n,s))

  print *, "Generating ", n, "Latin Hypercube Samples (type MCB)..."
  call lhs(n, s, "MCB", x)
  do i=1,n
    write(unit=*, fmt="(3F20.15)") x(i,:)
  end do
  print *, "done."

  print *, "Generating ", n, "Latin Hypercube Samples (type Patterson)..."
  call lhs(n, s, "Patterson", x)
  do i=1,n
    write(unit=*, fmt="(3F20.15)") x(i,:)
  end do
  print *, "done."

  print *, "Generating ", n, "Latin Hypercube Samples (type Modified)..."
  call lhs(n, s, "Modified", x)
  do i=1,n
    write(unit=*, fmt="(3F20.15)") x(i,:)
  end do
  print *, "done."

  deallocate(x)


  call show_test_summary(nb_errors)

end program test_lhs
