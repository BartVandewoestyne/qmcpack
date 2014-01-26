program pub_mascagni_chi

  use qmcpack

  integer(kind=i4b)                                   :: n, s
  real(kind=qp), dimension(:,:), allocatable          :: x
  type(functionparams)                                :: params
  type(integrationbounds)                             :: bounds
  type(integration_result), dimension(:), allocatable :: myres
  integer(kind=i4b), dimension(:), allocatable        :: startindex

  print *, "#################### REPRODUCING ARTICLE RESULTS #################"
  print *, "Article: Mascagni, M. and Chi, H."

  n = 50000

  call start_test("Reproducing Halton results...")
  s = 13
  print *, "Dimension = ", s
  allocate(x(n,s), myres(n), startindex(s))
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1
  call halton(n, s, startindex, primes(s), x)
  call sim_mc(f_fox86_np, transpose(x), params, f_fox86_exact(params, bounds), myres)
  call assert(myres(1000)%i, 1.171_qp, 0.001_qp)
  call assert(myres(2000)%i, 1.091_qp, 0.001_qp)
!  call assert(myres(3000)%i, 1.091_qp, 0.001_qp)
  call assert(myres(5000)%i, 0.978_qp, 0.001_qp)
  call assert(myres(7000)%i, 0.922_qp, 0.001_qp)
  call assert(myres(30000)%i, 0.979_qp, 0.001_qp)
  call assert(myres(40000)%i, 0.974_qp, 0.001_qp)
  call assert(myres(50000)%i, 0.984_qp, 0.001_qp)
  deallocate(x, myres, startindex)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)


  s = 20
  print *, "Dimension = ", s
  allocate(x(n,s), myres(n), startindex(s))
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1
  call halton(n, s, startindex, primes(s), x)
  call sim_mc(f_fox86_np, transpose(x), params, f_fox86_exact(params, bounds), myres)
  call assert(myres(1000)%i, 2.324_qp, 0.001_qp)
  call assert(myres(2000)%i, 1.444_qp, 0.001_qp)
  call assert(myres(3000)%i, 1.362_qp, 0.001_qp)
  call assert(myres(5000)%i, 1.140_qp, 0.001_qp)
  call assert(myres(7000)%i, 0.998_qp, 0.001_qp)
  call assert(myres(30000)%i, 0.888_qp, 0.001_qp)
  call assert(myres(40000)%i, 0.889_qp, 0.001_qp)
  call assert(myres(50000)%i, 0.903_qp, 0.001_qp)
  deallocate(x, myres, startindex)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)


  s = 25
  print *, "Dimension = ", s
  allocate(x(n,s), myres(n), startindex(s))
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1
  call halton(n, s, startindex, primes(s), x)
  call sim_mc(f_fox86_np, transpose(x), params, f_fox86_exact(params, bounds), myres)
  call assert(myres(1000)%i, 34.513_qp, 0.001_qp)
  call assert(myres(2000)%i, 17.450_qp, 0.001_qp)
  call assert(myres(3000)%i, 12.178_qp, 0.001_qp)
  call assert(myres(5000)%i, 7.811_qp, 0.001_qp)
  call assert(myres(7000)%i, 5.782_qp, 0.001_qp)
  call assert(myres(30000)%i, 2.13_qp, 0.01_qp)
  call assert(myres(40000)%i, 1.796_qp, 0.001_qp)
  call assert(myres(50000)%i, 1.568_qp, 0.001_qp)
  deallocate(x, myres, startindex)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)

  s = 40
  print *, "Dimension = ", s
  allocate(x(n,s), myres(n), startindex(s))
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1
  call halton(n, s, startindex, primes(s), x)
  call sim_mc(f_fox86_np, transpose(x), params, f_fox86_exact(params, bounds), myres)
  call assert(myres(1000)%i, 681382.379_qp, 0.001_qp)
  call assert(myres(2000)%i, 340691.207_qp, 0.001_qp)
  call assert(myres(3000)%i, 227127.541_qp, 0.001_qp)
  call assert(myres(5000)%i, 136276.627_qp, 0.001_qp)
  call assert(myres(7000)%i, 97340.706_qp, 0.001_qp)
  call assert(myres(30000)%i, 22713.137_qp, 0.01_qp)
  call assert(myres(40000)%i, 17035.076_qp, 0.001_qp)
  call assert(myres(50000)%i, 13628.735_qp, 0.001_qp)
  deallocate(x, myres, startindex)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)


  print *, "Reproducing DHalton results..."
  s = 13
  print *, "Dimension = ", s
  allocate(x(n,s), myres(n), startindex(s))
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1
  call halton_chi_permuted(n, s, startindex, x)
  call sim_mc(f_fox86_np, transpose(x), params, f_fox86_exact(params, bounds), myres)
  call assert(myres(1000)%i, 0.875_qp, 0.001_qp)
!  call assert(myres(2000)%i, 0.922_qp, 0.001_qp)
  call assert(myres(3000)%i, 0.908_qp, 0.001_qp)
  call assert(myres(5000)%i, 0.942_qp, 0.001_qp)
  call assert(myres(7000)%i, 0.942_qp, 0.001_qp)
  call assert(myres(30000)%i, 0.988_qp, 0.001_qp)
  call assert(myres(40000)%i, 1.014_qp, 0.001_qp)
  call assert(myres(50000)%i, 1.006_qp, 0.001_qp)
  deallocate(x, myres, startindex)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)


  s = 20
  print *, "Dimension = ", s
  allocate(x(n,s), myres(n), startindex(s))
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1
  call halton_chi_permuted(n, s, startindex, x)
  call sim_mc(f_fox86_np, transpose(x), params, f_fox86_exact(params, bounds), myres)
  call assert(myres(1000)%i, 0.601_qp, 0.001_qp)
!  call assert(myres(2000)%i, 0.952_qp, 0.001_qp)
  call assert(myres(3000)%i, 0.869_qp, 0.001_qp)
  call assert(myres(5000)%i, 0.985_qp, 0.001_qp)
  call assert(myres(7000)%i, 1.216_qp, 0.001_qp)
  call assert(myres(30000)%i, 1.097_qp, 0.001_qp)
!  call assert(myres(40000)%i, 1.118_qp, 0.001_qp)
  call assert(myres(50000)%i, 1.116_qp, 0.001_qp)
  deallocate(x, myres, startindex)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)


  s = 25
  print *, "Dimension = ", s
  allocate(x(n,s), myres(n), startindex(s))
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1
  call halton_chi_permuted(n, s, startindex, x)
  call sim_mc(f_fox86_np, transpose(x), params, f_fox86_exact(params, bounds), myres)
  call assert(myres(1000)%i, 0.612_qp, 0.001_qp)
!  call assert(myres(2000)%i, 0.846_qp, 0.001_qp)
  call assert(myres(3000)%i, 0.769_qp, 0.001_qp)
  call assert(myres(5000)%i, 1.979_qp, 0.001_qp)
  call assert(myres(7000)%i, 1.742_qp, 0.001_qp)
  call assert(myres(30000)%i, 1.171_qp, 0.01_qp)
  call assert(myres(40000)%i, 1.381_qp, 0.001_qp)
  call assert(myres(50000)%i, 1.289_qp, 0.001_qp)
  deallocate(x, myres, startindex)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)

  s = 40
  print *, "Dimension = ", s
  allocate(x(n,s), myres(n), startindex(s))
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1 ! Note: changing this to 0 gives totally different results...
  call halton_chi_permuted(n, s, startindex, x)
  call sim_mc(f_fox86_np, transpose(x), params, f_fox86_exact(params, bounds), myres)
  call assert(myres(1000)%i, 0.311_qp, 0.001_qp)
!  call assert(myres(2000)%i, 0.255_qp, 0.001_qp)
  call assert(myres(3000)%i, 0.515_qp, 0.001_qp)
  call assert(myres(5000)%i, 0.419_qp, 0.001_qp)
  call assert(myres(7000)%i, 0.489_qp, 0.001_qp)
  call assert(myres(30000)%i, 1.276_qp, 0.01_qp)
!  call assert(myres(40000)%i, 1.118_qp, 0.001_qp)
  call assert(myres(50000)%i, 1.034_qp, 0.001_qp)
  deallocate(x, myres, startindex)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)

  call show_test_summary(get_nb_assert_errors())

end program pub_mascagni_chi
