! References:
!
!   [1] `On computing the lattice rule criterion R', Joe, S. and Sloan, I. H.,
!       Mathematics of Computation, vol. 59, pages 557-568, 1992.

program pub_sloan92joe

  use qmcpack

  integer(kind=i4b)               :: nb_errors
  integer(kind=i4b), dimension(7) :: z
  integer(kind=i4b)               :: a
  integer(kind=i4b)               :: N

  call show_test_header("Article sloan92joe")

  nb_errors = 0

  N = 15019
  a = 12439
  z = create_korobov_vector(N, a, 7)
  print *, "z = ", z
  !write(unit=*, fmt="(A, F8.2)") "R(z, N) = ", R(z, N)
  write(unit=*, fmt="(A, F6.3)") "P_2(z, N) = ", P_2(z, N)

  N = 18101
  a = 17487
  z = create_korobov_vector(N, a, 7)
  !write(unit=*, fmt="(A, F8.2)") "R(z, N) = ", R(z, N)
  write(unit=*, fmt="(A, F6.3)") "P_2(z, N) = ", P_2(z, N)

  N = 24041
  a = 1833
  z = create_korobov_vector(N, a, 7)
  !write(unit=*, fmt="(A, F8.2)") "R(z, N) = ", R(z, N)
  write(unit=*, fmt="(A, F6.3)") "P_2(z, N) = ", P_2(z, N)

  N = 33139
  a = 7642
  z = create_korobov_vector(N, a, 7)
  !write(unit=*, fmt="(A, F8.2)") "R(z, N) = ", R(z, N)
  write(unit=*, fmt="(A, F6.3)") "P_2(z, N) = ", P_2(z, N)

  N = 46213
  a = 37900
  z = create_korobov_vector(N, a, 7)
  !write(unit=*, fmt="(A, F8.2)") "R(z, N) = ", R(z, N)
  write(unit=*, fmt="(A, F6.3)") "P_2(z, N) = ", P_2(z, N)

  N = 57091
  a = 35571
  z = create_korobov_vector(N, a, 7)
  !write(unit=*, fmt="(A, F8.2)") "R(z, N) = ", R(z, N)
  write(unit=*, fmt="(A, F6.3)") "P_2(z, N) = ", P_2(z, N)

  N = 71053
  a = 31874
  z = create_korobov_vector(N, a, 7)
  !write(unit=*, fmt="(A, F8.2)") "R(z, N) = ", R(z, N)
  write(unit=*, fmt="(A, F6.3)") "P_2(z, N) = ", P_2(z, N)

  N = 100063
  a = 39040
  z = create_korobov_vector(N, a, 7)
  !write(unit=*, fmt="(A, F8.2)") "R(z, N) = ", R(z, N)
  write(unit=*, fmt="(A, F6.3)") "P_2(z, N) = ", P_2(z, N)

  call show_test_summary(nb_errors)

end program pub_sloan92joe
