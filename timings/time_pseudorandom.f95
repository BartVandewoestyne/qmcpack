program time_pseudorandom

  use qmcpack

  integer(kind=i4b)                        :: n, s
  real(kind=qp), dimension(:), allocatable :: x_1p
  integer(kind=i4b)                        :: i
  real                                     :: t1, t2

  call get_unit(unit_number)
  open(unit=unit_number, file="pseudorandom_timing.dat", &
       status="replace", access="sequential", action="write")
  write(unit=unit_number, fmt="(A12, A10, A25)") "# nb of points", "dimension", "time"

  do

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the number of points you want to time for: "
    read *, n

    if (n==-1) then
      exit
    end if

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the dimension of your pointset: "
    read *, s

    allocate(x_1p(s))

    call random_seed()

    call cpu_time(t1)

    do i=1,n
      call random_number(x_1p)
    end do

    call cpu_time(t2)

    write (unit=*, fmt="(I12, A, I10, A, F14.6, A)") n, " ", s, &
           "-dimensional pseudorandom-points generated in ", t2-t1, " seconds."
    write (unit=unit_number,fmt="(I12, I10, ES25.6)") n, s, t2-t1

    deallocate(x_1p)

  end do

  close(unit=unit_number)

end program time_pseudorandom
