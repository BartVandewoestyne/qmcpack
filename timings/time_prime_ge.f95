! Testprogram to see how fast the prime_ge procedure is.
!
program time_prime_ge

  use qmcpack

  integer(kind=i4b) :: n
  integer(kind=i4b) :: myprime
  integer(kind=i4b) :: i
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
      myprime = prime_ge(i)
    end do
    call cpu_time(t2)
    write (unit=unit_number, fmt="(A, F25.6, A)") &
      "Time: ", t2-t1, " seconds."

  end do

end program time_prime_ge
