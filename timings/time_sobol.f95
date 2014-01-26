! Testprogram to time the generation of Sobol points.

program time_sobol

  use qmcpack

  integer(kind=i4b)                            :: n, s
  real(kind=qp), dimension(:), allocatable     :: x_1p
  type(soboltype)                              :: mysoboltype
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b)                            :: i
  real                                         :: t1, t2

  call get_unit(unit_number)
  open(unit=unit_number, file="sobol_timing.dat", &
       status="replace", access="sequential", action="write")
  write(unit=unit_number, fmt="(A12, A10, A25)") "# nb of points", "dimension", "time"

  do

    write(unit=*, fmt="(A)", advance="no") "Enter number of points: "
    read *, n

    if (n==-1) then
      exit
    end if

    write(unit=*, fmt="(A)", advance="no") "Enter dimension: "
    read *, s

    ! Set the type of Sobol sequence we want
    mysoboltype%poly_order = "JoeKuo"
    mysoboltype%dirnumbers = "JoeKuo"
    mysoboltype%use_ones = .true.
    mysoboltype%use_antonov_saleev = .true.

    allocate(startindex(s), x_1p(s))
    startindex = 1
    call init_sobol(n, s, startindex, mysoboltype)

    call cpu_time(t1)

    do i=1,n
      call next_sobol(x_1p)
    end do

    call cpu_time(t2)

    write (unit=*, fmt="(A, F10.2, A)") &
      "Time: ", t2-t1, " seconds."
    write (unit=unit_number,fmt="(I12, I10, ES25.6)") n, s, t2-t1

    deallocate(startindex, x_1p)

  end do

  close(unit=unit_number)

end program time_sobol
