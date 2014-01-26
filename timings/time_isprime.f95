! Testprogram to see how fast the isprime procedure is.
!
program time_isprime

  use qmcpack

  integer(kind=i4b) :: n
  integer(kind=i4b) :: i, nb_found
  real              :: t1, t2

  do

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the number of points you want to time for: "
    read *, n
    if (n==-1) then
      exit
    end if

    nb_found = 0
    call cpu_time(t1)
    do i=1,n
      if (isprime(i)) then
        nb_found = nb_found + 1
      end if
    end do
    call cpu_time(t2)
    write (unit=unit_number, fmt="(A, I0.0, A, F25.6, A)") &
      "Found ", nb_found, " primes in ", t2-t1, " seconds."

  end do

end program time_isprime
