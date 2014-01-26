! TODO: adapt this to the new interface for the models...

program test_system_model

  use qmcpack

!  implicit none

  real(kind=qp), dimension(:), allocatable :: x
  integer(kind=i4b)                        :: i
  integer(kind=i4b)                        :: nb_errors
  integer(kind=i4b)                        :: my_unit
  integer                                  :: ios

  call show_test_header("MOD_SYSTEM_MODEL")

  call random_seed()

  nb_errors = 0

  write(unit=*, fmt="(A)", advance="no") &
    "Testing next_gordon93salmond(x) by writing 50 samples to file " // &
    "gordon93samlond.dat..."

  call get_unit(my_unit)
  open(unit=my_unit, file="gordon93salmond.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  if (ios /= 0) then
    write(unit=*, fmt="(A)") "ERROR while opening gordon93salmond.dat!"
    nb_errors = nb_errors + 1
  else
    allocate(x(1))
    x(1) = 0.1_qp
    write(unit=my_unit, fmt="(ES25.15)") x
    do i=1,50
!      call next_gordon93salmond(x)
      write(unit=my_unit, fmt="(ES25.15)") x
    end do
  end if

  deallocate(x)
  close(unit=my_unit)
  write(unit=*, fmt="(A)") "done."

  call show_test_summary(nb_errors)

end program test_system_model
