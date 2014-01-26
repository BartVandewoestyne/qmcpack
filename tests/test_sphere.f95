program test_sphere

  use qmcpack

  real(kind=qp), dimension(2) :: c
  real(kind=qp), dimension(3) :: x
  real(kind=qp), dimension(3) :: g
  real(kind=qp)               :: mysum
  integer(kind=i4b)           :: i
  integer(kind=i4b)           :: N
  integer(kind=i4b)           :: nb_errors
  integer(kind=i4b)           :: my_unit

  call show_test_header("MOD_SPHERE")

  nb_errors = 0

  N = 10000

  write(unit=*, fmt="(A)", advance="no") &
    "Testing transform_to_sphere3D(...) by writing Halton-based samples to file transform_to_sphere3D_halton.dat..."
  call get_unit(my_unit)
  call checked_open(my_unit, "transform_to_sphere3D_Halton.dat", "write")
  call init_halton(2)
  do i=1,N
    call next_halton(c)
    call transform_to_sphere3D_fang_wang(c, x)
    write(unit=my_unit, fmt="(3ES25.15)") x
  end do
  call checked_close(my_unit)
  call free_halton()
  write(unit=*, fmt="(A)") "done."


  call random_seed()
  write(unit=*, fmt="(A)", advance="no") &
    "Testing transform_to_sphere3D(...) by writing pseudorandom samples to file transform_to_sphere3D_random.dat..."
  call get_unit(my_unit)
  call checked_open(my_unit, "transform_to_sphere3D_random.dat", "write")
  do i=1,N
    call random_number(c)
    call transform_to_sphere3D_fang_wang(c, x)
    write(unit=my_unit, fmt="(3ES25.15)") x
  end do
  call checked_close(my_unit)
  write(unit=*, fmt="(A)") "done."

  N = 1000000

  ! Example of an elliptic integral.
  call get_unit(my_unit)
  call checked_open(my_unit, "elliptic_integral_random.dat", "write")
  g = (/ 4.0_qp/3, 4.0_qp/3, 2.0_qp/3 /)
  mysum = 0
  do i=1,N
    call random_number(c)
    call transform_to_sphere3D_fang_wang(c, x)
    mysum = mysum + sqrt(sum(g*x*x))
    write(unit=my_unit, fmt="(ES25.15)") 1.049523-mysum/i
  end do
  print *, "INTEGRAL       RANDOM = ", mysum/N
  call checked_close(my_unit)

  ! Example of an elliptic integral.
  call get_unit(my_unit)
  call checked_open(my_unit, "elliptic_integral_quasirandom.dat", "write")
  g = (/ 4.0_qp/3, 4.0_qp/3, 2.0_qp/3 /)
  mysum = 0
  call init_halton(2)
  do i=1,N
    call next_halton(c)
    call transform_to_sphere3D_fang_wang(c, x)
    mysum = mysum + sqrt(sum(g*x*x))
    write(unit=my_unit, fmt="(ES25.15)") 1.049523-mysum/i
  end do
  print *, "INTEGRAL QUASI-RANDOM = ", mysum/N
  call free_halton()
  call checked_close(my_unit)

  call show_test_summary(nb_errors)

end program test_sphere
