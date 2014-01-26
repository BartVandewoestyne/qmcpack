program test_sort

  use qmcpack

  real(kind=qp), dimension(10)                 :: x
  real(kind=qp), dimension(:,:), allocatable   :: x_np
  integer(kind=i4b), dimension(10)             :: indices
  integer(kind=i4b), dimension(:), allocatable :: startindex, myprimes
  integer(kind=i4b)                            :: i
  integer(kind=i4b)                            :: n, s
  integer(kind=i4b)                            :: nb_errors

  call show_test_header("MOD_SORT")

  nb_errors = 0

  call random_seed()

  write(unit=*, fmt="(A)") "Testing quicksort() with indices..."

  call random_number(x)
  indices = (/ (i, i=1,10) /)

  write(unit=*, fmt="(A)") "  Before sorting:"
  write(unit=*, fmt="(A, 10F10.5)") "  ", x
  write(unit=*, fmt="(A, 10I10)") "  ", indices

  call quicksort(x, indices)

  write(unit=*, fmt="(A)") "  After sorting:"
  write(unit=*, fmt="(A, 10F10.5)") "  ", x
  write(unit=*, fmt="(A, 10I10)") "  ", indices
  if (.not. is_increasing(x)) then
    write(unit=*, fmt="(A)") "ERROR: not correctly sorted!"
    nb_errors = nb_errors + 1
  end if

  write(unit=*, fmt="(A)") "done."


  write(unit=*, fmt="(A)") "Testing quicksort() without indices..."

  call random_number(x)

  write(unit=*, fmt="(A)") "  Before sorting:"
  write(unit=*, fmt="(A, 10F10.5)") "  ", x

  call quicksort(x)

  write(unit=*, fmt="(A)") "  After sorting:"
  write(unit=*, fmt="(A, 10F10.5)") "  ", x
  if (.not. is_increasing(x)) then
    write(unit=*, fmt="(A)") "ERROR: not correctly sorted!"
    nb_errors = nb_errors + 1
  end if

  write(unit=*, fmt="(A)") "done."

  
  write(unit=*, fmt="(A)", advance="no") "Testing quicksort() on a 2D square root sequence... "
  n = 1000
  s = 3
  allocate(x_np(n,s), startindex(s), myprimes(s))
  startindex = 1
  myprimes = primes(s)
  call square_root_seq(n, s, startindex, myprimes, x_np)
  call quicksort(x_np(:,1))
  if (.not. is_increasing(x_np(:,1))) then
    write(unit=*, fmt="(A)") "ERROR: not correctly sorted!"
    nb_errors = nb_errors + 1
  end if
  deallocate(x_np, startindex, myprimes)
  write(unit=*, fmt="(A)") "done."

  call show_test_summary(nb_errors)

end program test_sort
