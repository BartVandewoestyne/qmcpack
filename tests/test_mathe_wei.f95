program test_mathe_wei

  use qmcpack

  integer(kind=i4b)                          :: n, d, s, max_size
  real(kind=qp), dimension(:,:), allocatable :: x
  real(kind=qp)                              :: res
  integer                                    :: nberrors
  type(functionparams)                       :: dummy_params

  call show_test_header("MATHE_WEI")

  nberrors = 0

  call random_seed()
  N = 1000000

  call init_mathe_wei(2.0_qp, 2.0_qp)

  write(unit=*, fmt="(A)") "==> First integral from Mathe and Wei's article:"
  s = 4
  d = 2
  call calculate_max_pointset_size(N, d, "ECD", s, max_size)
  allocate(x(d, max_size))
  call random_number(x)
  call mathe_wei(N, x, f_mathe_wei01_np, dummy_params, rho_mathe_wei01_np, dummy_params, "ECD", s, res)
  write(unit=*, fmt="(A16, F20.15)") "approximation = ", res
  write(unit=*, fmt="(A16, F20.15)") "exact value   = ", 2*log(2.0_qp)-1
  deallocate(x)
  print *


  write(unit=*, fmt="(A)") "==> Second integral from Mathe and Wei's article:"
  s = 4
  d = 2
  call calculate_max_pointset_size(N, d, "ECD", s, max_size)
  allocate(x(d, max_size))
  call random_number(x)
  call mathe_wei(N, x, f_mathe_wei02_np, dummy_params, rho_mathe_wei02_np, dummy_params, "ECD", s, res)
  write(unit=*, fmt="(A16, F20.15)") "approximation = ", res
  write(unit=*, fmt="(A16, F20.15)") "exact value   = ", 2*log(2.0_qp)-1
  deallocate(x)
  print *


  write(unit=*, fmt="(A)") "==> Integral with bivariate Normal from Tim Pillards his article:"
  s = 0 ! dummy value
  d = 2
  call calculate_max_pointset_size(N, d, "MVN", s, max_size)
  allocate(x(d, max_size))
  call random_number(x)
  call mathe_wei(N, x, f_sum_of_squares_np, dummy_params, f_bivariate_normal_np, dummy_params, "MVN", s, res)
  write(unit=*, fmt="(A16, F20.15)") "approximation = ", res
  deallocate(x)
  print *


  write(unit=*, fmt="(A)") "==> Integral with f=1 and rho the bivariate Normal:"
  s = 0 ! dummy value
  d = 2
  call calculate_max_pointset_size(N, d, "MVN", s, max_size)
  allocate(x(d, max_size))
  call random_number(x)
  call mathe_wei(N, x, f_one_np, dummy_params, f_bivariate_normal_np, dummy_params, "MVN", s, res)
  write(unit=*, fmt="(A16, F20.15)") "approximation = ", res
  write(unit=*, fmt="(A16, F20.15)") "exact value   = ", 1.0_qp
  deallocate(x)
  print *

  call show_test_summary(nberrors)

end program test_mathe_wei
