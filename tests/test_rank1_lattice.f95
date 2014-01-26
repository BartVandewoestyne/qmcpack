program test_rank1_lattice

  use qmcpack

  integer(kind=i4b)                            :: N, s
  integer(kind=i4b), dimension(:), allocatable :: z
  integer(kind=i4b)                            :: i
  real(kind=qp), dimension(:,:), allocatable   :: x
  integer(kind=i4b)                            :: nb_errors
  character(len=10)                            :: fmt_string
  integer                                      :: ios

  call show_test_header("MOD_RANK1_LATTICE")
  nb_errors = 0


  write(unit=*, fmt="(A)", advance="no") "Testing block assignment with rank1_lattice()..."
  N = 34
  s = 2
  allocate(x(N,s), z(s))
  z = (/ 1, 21 /)
  call rank1_lattice(N, z, x)
  deallocate(x, z)
  write(unit=*, fmt="(A)") "done."


  ! Write out a rank1 Fibonacci lattice of order 34 to a file.
  !
  write(unit=*, fmt="(A)", advance="no") "Writing out Fibonacci lattice rule of order 34 to file fibonacci_34.dat... "
  n = 34
  s = 2
  allocate(x(N,s), z(s))
  z = (/ 1, 21 /)
  call get_unit(unit_number)
  open(unit=unit_number, file="fibonacci_34.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F20.15)"
  call rank1_lattice(N, z, x)
  do i=1,N
    write(unit=unit_number, fmt=fmt_string) x(i,:)
  end do
  deallocate(x, z)
  close(unit=unit_number)
  write(unit=*, fmt="(A)") "done."


  ! Write out a rank1 Fibonacci lattice of order 610 to a file.
  !
  write(unit=*, fmt="(A)", advance="no") "Writing out Fibonacci lattice rule of order 610 to file fibonacci_610.dat... "
  n = 610
  s = 2
  allocate(x(N,s), z(s))
  z = (/ 1, 377 /)
  call get_unit(unit_number)
  open(unit=unit_number, file="fibonacci_610.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F20.15)"
  call rank1_lattice(N, z, x)
  do i=1,N
    write(unit=unit_number, fmt=fmt_string) x(i,:)
  end do
  deallocate(x, z)
  close(unit=unit_number)
  write(unit=*, fmt="(A)") "done."


  ! Write out a 3D Korobov lattice to a file.
  ! See Maissonneuve paper page 158.
  !
  write(unit=*, fmt="(A)", advance="no") "Writing out Korobov lattice rule to file korobov.dat... "
  !n = 101
  !n = 199
  !n = 307
  !n = 523
  n = 701
  s = 3
  allocate(x(N,s), z(s))
  !z = (/ 1, 40,   85 /)
  !z = (/ 1, 30,  104 /)
  !z = (/ 1, 75,   99 /)
  !z = (/ 1, 78,  331 /)
  z = (/ 1, 215, 660 /)
  call get_unit(unit_number)
  open(unit=unit_number, file="korobov.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F18.15)"
  call rank1_lattice(N, z, x)
  do i=1,N
    write(unit=unit_number, fmt=fmt_string) x(i,:)
  end do
  deallocate(x, z)
  close(unit=unit_number)
  write(unit=*, fmt="(A)") "done."


  call show_test_summary(nb_errors)

end program test_rank1_lattice
