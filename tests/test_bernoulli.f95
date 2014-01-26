program test_bernoulli

  use qmcpack

  real(kind=qp)                            :: res
  real(kind=qp), dimension(:), allocatable :: b

  call show_test_header("MOD_BERNOULLI")
  call init_assert()
  

  call start_test("Testing bernoulli2(n,b)")
  allocate(b(1))
  call bernoulli2(0, b)
  call assert(b, (/ 1.0_qp /))
  deallocate(b)
  allocate(b(2))
  call bernoulli2(1, b)
  call assert(b, (/ 1.0_qp, -0.5_qp /))
  deallocate(b)
  allocate(b(3))
  call bernoulli2(2, b)
  call assert(b, (/ 1.0_qp, -0.5_qp, 1.0_qp/6 /))
  deallocate(b)
  allocate(b(4))
  call bernoulli2(3, b)
  call assert(b, (/ 1.0_qp, -0.5_qp, 1.0_qp/6, 0.0_qp /))
  deallocate(b)
  allocate(b(5))
  call bernoulli2(4, b)
  call assert(b, (/ 1.0_qp, -0.5_qp, 1.0_qp/6, 0.0_qp, -1.0_qp/30 /))
  deallocate(b)
  call stop_test()


  call start_test("Testing bernoulli3(n)...")
  call assert(bernoulli3(0), 1.0_qp)
  call assert(bernoulli3(1), -1.0_qp/2)
  call assert(bernoulli3(2), 1.0_qp/6)
  call assert(bernoulli3(3), 0.0_qp)
  call assert(bernoulli3(4), -1.0_qp/30)
  call assert(bernoulli3(5), 0.0_qp)
  call assert(bernoulli3(6), 1.0_qp/42)
  call assert(bernoulli3(7), 0.0_qp)
  call assert(bernoulli3(8), -1.0_qp/30)
  call assert(bernoulli3(9), 0.0_qp)
  call assert(bernoulli3(10), 5.0_qp/66)
  call assert(bernoulli3(11), 0.0_qp)
  call stop_test()


  call start_test("Testing bernoulli_poly(n,x)...")
  call assert(bernoulli_poly(2, 0.5_qp), -1.0_qp/12)
  call assert(bernoulli_poly(2, 0.25_qp), -1.0_qp/48)
  call stop_test()


  call start_test("Testing bernoulli_poly2(n,x,res)...")
  call bernoulli_poly2(2, 0.5_qp, res)
  call assert(res, -1.0_qp/12)
  call bernoulli_poly2(2, 0.25_qp, res)
  call assert(res, -1.0_qp/48)
  call stop_test()


  call show_test_summary(get_nb_assert_errors())

end program test_bernoulli
