! Testprogram to see how fast the radical inverse routines are.
!
program time_get_degree

  use qmcpack

  integer(kind=i4b) :: n
  integer(kind=i4b) :: i
  real              :: t1, t2
  integer(kind=i4b) :: deg

  do

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the number of points you want to time for: "
    read *, n
    if (n==-1) then
      exit
    end if

    call cpu_time(t1)
    do i=1,n
      deg = get_degree(i)
    end do
    call cpu_time(t2)
    write (unit=unit_number,fmt="(A, F25.6)") &
      "Time: ", t2-t1

  end do

end program time_get_degree
