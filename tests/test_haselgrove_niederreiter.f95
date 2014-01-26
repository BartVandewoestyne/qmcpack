program test_haselgrove_niederreiter

  use qmcpack

  integer(kind=i4b)                        :: q
  type(functionparams)                     :: params
  real(kind=qp), dimension(:), allocatable :: irrationals

  call show_test_header("MOD_HASELGROVE_NIEDERREITER")
  call init_assert()


  call start_test("Testing init_nied_weyl()...")
  q = 3
  allocate(irrationals(3))
  irrationals = sqrt(real(primes(3), kind=qp))
  call init_nied_weyl(q, irrationals, f_warnock02_1p, params, "none")
  !call init_nied_weyl(q, irrationals, f_simple_product_1p, params, "none")
  !call init_nied_q3_weyl(irrationals, f_warnock02_1p, params, "none")
  deallocate(irrationals)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_haselgrove_niederreiter
