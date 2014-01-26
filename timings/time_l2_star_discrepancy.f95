! Program to time the speed of L2 star discrepancy calculations.

program time_l2_star_discrepancy

  use qmcpack

  real(kind=qp), dimension(:,:), allocatable :: x
  real(kind=qp), dimension(:,:), allocatable :: y
  real(kind=qp), dimension(:), allocatable   :: tn2
  integer(kind=i4b)                          :: s, N
  real                                       :: t1, t2

  do

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the number of points you want to time for: "
    read *, N

    if (n==-1) then
      exit
    end if

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the dimension of your pointset: "
    read *, s

    allocate(x(s,N), y(N,s))
    allocate(tn2(N))

    ! We will calculate the L2 star discrepancy of a
    ! pseudorandom pointset.
    call random_number(x)

    call cpu_time(t1)
    call T_N_star_squared(x, tn2)
    call cpu_time(t2)
    write (unit=*, fmt="(A, ES20.10)") &
          "T_N_star_squared result = ", tn2(N)
    write (unit=*, fmt="(A, F7.2, A)") &
          "T_N_star_squared time: ", t2-t1, " seconds."

    y = transpose(x)
    call cpu_time(t1)
    call discd2(N, s, y, tn2)
    call cpu_time(t2)
    write (unit=*, fmt="(A, ES20.10)") &
          "discd2 result = ", tn2(N)
    write (unit=*, fmt="(A, F7.2, A)") &
          "discd2 time: ", t2-t1, " seconds."

    deallocate(x, y, tn2)

  end do

end program time_l2_star_discrepancy
