! References:
!
!   [1] `Recherche et utilisation des bons trellis.  Programmation et
!        resultats numeriques', Maisonneuve, D., Proceedings of the Symposium
!        on Applications of Number Theory to Numerical Analysis (Montreal
!        1971), eds. Zaremba, S. K., Academic Press, New York, pages 121-201.

program pub_maisonneuve

  use qmcpack

  integer(kind=i4b), dimension(3) :: z
  integer(kind=i4b)               ::  N

  call show_test_header("Article Maisonneuve")
  call init_assert()

  N = 101
  z = (/ 1, 40, 85 /)
  call assert(P_2(z, N), 2.6E-1_qp, 0.1E-1_qp)
  call assert(P_4(z, N), 1.7E-3_qp, 0.1E-3_qp)

  N = 199
  z = (/ 1, 30, 104 /)
  call assert(P_2(z, N), 8.2E-2_qp, 0.1E-2_qp)
  call assert(P_4(z, N), 1.3E-4_qp, 0.1E-4_qp)

  N = 307
  z = (/ 1, 75, 99 /)
  call assert(P_2(z, N), 4.2E-2_qp, 0.1E-2_qp)
  call assert(P_4(z, N), 3.0E-5_qp, 0.1E-5_qp)

  N = 523
  z = (/ 1, 78, 331 /)
  call assert(P_2(z, N), 1.8E-2_qp, 0.1E-2_qp)
  call assert(P_4(z, N), 6.6E-6_qp, 0.1E-6_qp)

  N = 701
  z = (/ 1, 215, 660 /)
  call assert(P_2(z, N), 1.2E-2_qp, 0.1E-2_qp)
  call assert(P_4(z, N), 3.2E-6_qp, 0.1E-6_qp)

  N = 1069
  z = (/ 1, 136, 323 /)
  call assert(P_2(z, N), 5.7E-3_qp, 0.1E-3_qp)
  call assert(P_4(z, N), 1.0E-6_qp, 0.1E-6_qp)

  N = 1543
  z = (/ 1, 355, 1042 /)
  call assert(P_2(z, N), 3.0E-3_qp, 0.1E-3_qp)
  call assert(P_4(z, N), 1.7E-7_qp, 0.1E-7_qp)

  N = 2129
  z = (/ 1, 359, 1141 /)
  call assert(P_2(z, N), 1.7E-3_qp, 0.1E-3_qp)
  call assert(P_4(z, N), 4.6E-8_qp, 0.1E-8_qp)

  N = 3001
  z = (/ 1, 276, 1151 /)
  call assert(P_2(z, N), 1.0E-4_qp, 0.1E-4_qp)
  call assert(P_4(z, N), 1.7E-8_qp, 0.1E-8_qp)

  N = 4001
  z = (/ 1, 722, 1154 /)
  call assert(P_2(z, N), 5.9E-4_qp, 0.1E-4_qp)
  call assert(P_4(z, N), 6.1E-9_qp, 0.1E-9_qp)

  N = 5003
  z = (/ 1, 1476, 2271 /)
  call assert(P_2(z, N), 4.4E-4_qp, 0.1E-4_qp)
  call assert(P_4(z, N), 4.0E-9_qp, 0.1E-9_qp)

  N = 6007
  z = (/ 1, 592, 2058 /)
  call assert(P_2(z, N), 2.8E-4_qp, 0.1E-4_qp)
  call assert(P_4(z, N), 1.2E-9_qp, 0.1E-9_qp)

  N = 8191
  z = (/ 1, 739, 5515 /)
  call assert(P_2(z, N), 1.7E-4_qp, 0.1E-4_qp)
  call assert(P_4(z, N), 4.0E-10_qp, 0.1E-10_qp)

  N = 10007
  z = (/ 1, 544, 5733 /)
  call assert(P_2(z, N), 1.3E-4_qp, 0.1E-4_qp)
  call assert(P_4(z, N), 2.5E-10_qp, 0.1E-10_qp)

  N = 20039
  z = (/ 1, 5704, 12319 /)
  call assert(P_2(z, N), 6.4E-5_qp, 0.1E-5_qp)

  N = 28117
  z = (/ 1, 19449, 5600 /)
  call assert(P_2(z, N), 3.0E-5_qp, 0.1E-5_qp)

  N = 39029
  z = (/ 1, 10607, 26871 /)
  call assert(P_2(z, N), 2.1E-5_qp, 0.1E-5_qp)

  N = 57091
  z = (/ 1, 48188, 21101 /)
  call assert(P_2(z, N), 9.8E-6_qp, 0.1E-6_qp)

  N = 82001
  z = (/ 1, 21252, 67997 /)
  call assert(P_2(z, N), 4.1E-6_qp, 0.1E-6_qp)

  N = 100063
  z = (/ 1, 53584, 37334 /)
  call assert(P_2(z, N), 4.9E-6_qp, 0.1E-6_qp)

  call show_test_summary(get_nb_assert_errors())

end program pub_maisonneuve
