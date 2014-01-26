program test_cranley_patterson

  use qmcpack

!  implicit none

  integer(kind=i4b)                          :: n, s
  integer(kind=i4b)                          :: i
  real(kind=qp), dimension(:), allocatable   :: x_1p
  real(kind=qp), dimension(:,:), allocatable :: x_np
  integer(kind=i4b)                          :: nb_errors

  call show_test_header("CRANLEY PATTERSON")

  nb_errors = 0

  n = 1000
  s = 200

  write(unit=*, fmt="(A)", advance="no") "Testing random shift on single points... "
  allocate(x_1p(s))
  call random_number(x_1p)
  call random_shift(x_1p)
  if (any(x_1p>1) .or. any(x_1p<0)) then
    write(unit=*, fmt="(A)")
    write(unit=*, fmt="(A)") "  ERROR: found number outside [0,1]!!!"
    nb_errors = nb_errors + 1
  end if
  deallocate(x_1p)
  write(unit=*, fmt="(A)") "done."


  write(unit=*, fmt="(A)", advance="no") "Testing random shift on pointsets... "
  allocate(x_np(n, s))
  call random_number(x_np)
  call random_shift(x_np)
  if (any(x_np>1) .or. any(x_np<0)) then
    write(unit=*, fmt="(A)")
    write(unit=*, fmt="(A)") "  ERROR: found number outside [0,1]!!!"
    nb_errors = nb_errors + 1
  end if
  deallocate(x_np)
  write(unit=*, fmt="(A)") "done."


  write(unit=*, fmt="(A)", advance="no") "Checking if we go outside [0,1] when shifting Halton numbers... "
  call init_halton(s)
  call init_random_shift(s)
  allocate(x_1p(s))
  do i=1,10000
    call next_halton(x_1p)
    call random_shift(x_1p)
    if (any(x_1p>1) .or. any(x_1p<0)) then
      write(unit=*, fmt="(A)")
      write(unit=*, fmt="(A)") "  ERROR: found number outside [0,1]!!!"
      nb_errors = nb_errors + 1
    end if
  end do
  deallocate(x_1p)
  call free_halton()
  call free_cranley_patterson()
  write(unit=*, fmt="(A)") "done."
  

  call show_test_summary(nb_errors)

end program test_cranley_patterson
