program test_file_utils

  use qmcpack

  integer           :: my_unit
  integer(kind=i4b) :: nb_errors

  call show_test_header("MOD_FILE_UTILS")
  
  nb_errors = 0

  write(unit=*, fmt="(A)") "Testing get_unit()..."
  call get_unit(my_unit)
  write(unit=*, fmt="(A, I0)") "  Found unit ", my_unit
  write(unit=*, fmt="(A)") "done."

  call show_test_summary(nb_errors)

end program test_file_utils
