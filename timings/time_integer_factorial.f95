! Time the integer_factorial routine from mod_number_theory.

program time_integer_factorial

  use qmcpack

  integer(kind=i4b) :: f, n_max, nb_it
  integer(kind=i4b) :: i, j
  real              :: t1, t2

  n_max = 12
  nb_it = 1000000

  call cpu_time(t1)
  do i=1,nb_it
    do j=1,n_max
      f = integer_factorial_recursive(j)
    end do
  end do
  call cpu_time(t2)
  write (unit=*, fmt="(A, F10.2, A)") "Time: ", t2-t1, " seconds."

  call cpu_time(t1)
  do i=1,nb_it
    do j=1,n_max
      f = integer_factorial_nonrecursive(j)
    end do
  end do
  call cpu_time(t2)
  write (unit=*, fmt="(A, F10.2, A)") "Time: ", t2-t1, " seconds."

end program time_integer_factorial
