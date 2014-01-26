program test_special_functions

  use qmcpack

  real(kind=qp), dimension(:), allocatable :: w
  integer(kind=i4b)                        :: N
  real(kind=qp)                            :: u
  integer(kind=i4b)                        :: i

  N = 1000

  call show_test_header("MOD_SPECIAL_FUNCTIONS")
  call init_assert()


  call start_test("Testing normal cumulative distribution function phi(x)...")
  call assert(phi(0.0_qp), 0.5_qp)
  call assert(phi(1.0_qp), 0.8413447460685429_qp)
  call assert(phi(5.0_qp), 0.9999997133484281_qp)
  call write_function_data(phi, -5.0_qp, 5.0_qp, 10.0_qp/N, "phi.dat")
  call stop_test()


  call start_test("Testing erfc(x)...")
  call assert(erfc(0.0_qp), 1.0_qp)
  call assert(erfc(0.00001_qp), 0.99998871620832942100065005896470_qp)
  call assert(erfc(0.2_qp), 0.77729741078952154585986099319986_qp)
  call assert(erfc(0.4_qp), 0.57160764495333154489639615467983_qp)
  call assert(erfc(0.5_qp), 0.47950012218695346231725334610804_qp)
  call assert(erfc(0.75_qp), 0.28884436634648486840106216540859_qp)
  call assert(erfc(1.0_qp), 0.15729920705028513065877936491739_qp)
  call assert(erfc(5.0_qp), 0.15374597944280348501883434853834e-11_qp)
  call assert(erfc(10.0_qp), 0.20884875837625447570007862949578e-44_qp)
  call write_function_data(erfc, -4.0_qp, 4.0_qp, 8.0_qp/N, "erfc.dat")
  call stop_test()


  call start_test("Testing erf(x)...")
  call assert(erf(0.0_qp), 0.0_qp)
  call assert(erf(0.00001_qp), 0.11283791670578999349941035298254e-4_qp)
  call assert(erf(0.2_qp), 0.22270258921047845414013900680014_qp)
  call assert(erf(0.4_qp), 0.42839235504666845510360384532017_qp)
  call assert(erf(0.5_qp), 0.52049987781304653768274665389196_qp)
  call assert(erf(0.75_qp), 0.71115563365351513159893783459141_qp)
  call assert(erf(1.0_qp), 0.84270079294971486934122063508261_qp)
  call write_function_data(erf, -4.0_qp, 4.0_qp, 8.0_qp/N, "erf.dat")
  call stop_test()


  call start_test("Testing gamma(x)...")
  call assert(gamma(0.2_qp), 4.59084371199880_qp)
  call assert(gamma(8.5_qp), 14034.4072934834_qp)
  call write_function_data(gamma, 0.01_qp, 4.0_qp, 4.0_qp/N, "gamma.dat")
  call stop_test()

  call start_test("Testing norminv(x)...")
  call assert(norminv(0.25_qp), -0.67448975019608174320222701454131_qp)
  call assert(norminv(0.001_qp), -3.0902323061678135415403998301074_qp)
  call assert(norminv(1.0e-20_qp), -9.262340089798408_qp)
  call assert(norminv(0.1_qp), -1.2815515655446004669651033294487_qp)
  call write_function_data(norminv, 0.0_qp, 1.0_qp, 1.0_qp/N, "norminv.dat")
  call random_seed()
  allocate(w(N))
  do i = 1,N
    call random_number(u)
    w(i) = norminv(u)
  end do
  call write_vector(w, "norminv_samples.dat", "(F20.15)")
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_special_functions
