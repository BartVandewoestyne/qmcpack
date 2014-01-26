! Program to see how fast the routine that finds the rightmost zero bit
! position is.
!
program time_find_rightmost_zero_bitpos

  use qmcpack

  integer(kind=i4b) :: i, n
  integer(kind=i4b) :: pos
  real              :: t1, t2

  do

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the number of points you want to time for: "
    read *, n

    if (n==-1) then
      exit
    end if

    call cpu_time(t1)
    do i=1,n
      pos = find_rightmost_zerobitpos(i)
    end do
    call cpu_time(t2)
    write (unit=unit_number,fmt="(A, F10.2)") &
      "Time: ", t2-t1

  end do

end program time_find_rightmost_zero_bitpos
