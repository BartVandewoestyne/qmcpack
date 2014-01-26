program test_monte_carlo

  use qmcpack

  integer(kind=i4b)                                   :: n, s
  real(kind=qp), dimension(:,:), allocatable          :: x
  type(functionparams)                                :: params
  type(integrationbounds)                             :: bounds
  type(integration_result), dimension(:), allocatable :: myres
  integer(kind=i4b)                                   :: i
  integer(kind=i4b), dimension(:), allocatable        :: startindex
  integer(kind=i4b), dimension(:), allocatable        :: myprimes
  integer                                             :: my_unit

  call show_test_header("MOD_MONTE_CARLO")
  call init_assert()

  call start_test("Running a simulation and writing it to mc_result.dat...")
  n = 500
  s = 100
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 0.21_qp
  params%u = 0.90_qp
  call set_unit_hypercube(bounds)

  allocate(x(n,s), myres(n), startindex(s), myprimes(s))
  !call random_number(x)
  startindex = 1
  myprimes = primes(s)
  call square_root_seq(n, s, startindex, myprimes, x)

  call sim_mc(f_c0_np, transpose(x), params, f_c0_exact(params, bounds), myres)

  call get_unit(my_unit)
  call checked_open(my_unit, "mc_result.dat", "write")
  do i=1,n
    write(unit=my_unit, fmt="(I10,3(F20.15))") &
      i, myres(i)%i, myres(i)%abs_err, myres(i)%rel_err
  end do
  call checked_close(my_unit)

  deallocate(x, myres, startindex, myprimes)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)
  call stop_test()


  call show_test_summary(get_nb_assert_errors())

end program test_monte_carlo
