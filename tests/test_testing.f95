! Testfile for mod_testing

program test_testing

  use qmcpack

  integer(kind=i4b)  :: nb_columns, nb_points
  logical            :: equal
  character(len=100) :: fmt_string

  call show_test_header("MOD_TESTING")
  call init_assert()


  call start_test("Testing get_format_specifier(...)...")
  call get_format_specifier("test_data/pointset_random1_10rows_3columns.dat", &
                            fmt_string)
  if (fmt_string /= "(3F7.4)") then
    call increase_nb_assert_errors(1)
  end if
  call stop_test()


  call start_test("Testing get_nb_columns(...)...")
  call get_nb_columns("test_data/pointset_random_13rows_5columns.dat", &
                      nb_columns)
  call assert(nb_columns, 5)
  call get_nb_columns("test_data/FAURE_BRATLEY_FOX_40D_1025P.DAT", &
                      nb_columns)
  call assert(nb_columns, 40)

  call get_nb_columns("test_data/pointset_random_13rows_1column.dat", &
                      nb_columns)
  call assert(nb_columns, 1)
  call stop_test()


  call start_test("Testing get_nb_points(...)...")
  call get_nb_points("test_data/pointset_random_13rows_5columns.dat", &
                     nb_points)
  if (nb_points /= 13) then
    call increase_nb_assert_errors(1)
  end if
  call stop_test()


  call start_test("Testing compare_pointsets(...)...")
  call compare_pointsets("test_data/pointset_random1_10rows_3columns.dat", &
                         "test_data/pointset_random2_10rows_3columns.dat", &
                         0.0001_qp, equal)
  if (.not. equal) then
    call increase_nb_assert_errors(1)
  end if
  call compare_pointsets("test_data/pointset_random_13rows_5columns.dat", &
                         "test_data/pointset_random_13rows_5columns_b.dat", &
                         0.0001_qp, equal)
  if (equal) then
    call increase_nb_assert_errors(1)
  end if
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_testing
