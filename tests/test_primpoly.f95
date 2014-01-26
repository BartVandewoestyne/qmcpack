program test_primpoly

  use qmcpack

  call show_test_header("MOD_PRIMPOLY")
  call init_assert()

  call start_test("Testing get_degree()...")
  call assert(get_degree(11), 3)
  call assert(get_degree(0), 0)
  call assert(get_degree(1), 0)
  call assert(get_degree(2), 1)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_primpoly
