program test_testfunctions

  use qmcpack

  real(kind=qp), dimension(:), allocatable   :: x_1p
  real(kind=qp), dimension(:,:), allocatable :: x_np
  real(kind=qp), dimension(:), allocatable   :: y
  type(functionparams)                       :: params
  type(integrationbounds)                    :: bounds
  integer(kind=i4b)                          :: N, s

  call show_test_header("TESTFUNCTIONS")
  call init_assert()

  ! Test evaluation in more than one point and dimension
  N = 3
  s = 2

  ! Setup the parameters...
  allocate(params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 0.2_qp
  params%u = 0.1_qp
  call set_unit_hypercube(bounds)

  ! The point we're evaluating
  allocate(x_1p(s), x_np(s,N), y(N))
  x_np(1,:) = (/ 0.1_qp, 0.2_qp, 0.3_qp /)
  x_np(2,:) = (/ 0.4_qp, 0.5_qp, 0.6_qp /)
  x_1p = (/ 0.1_qp, 0.4_qp /)

  call start_test("Testing the Gaussian from Braaten and Weller's paper...")
  y =  f_gaussian_bw(x_np)
  call assert(y(:), (/ 0.7117703227626098_qp, &
                                0.8352702114112720_qp, &
                                0.9048374180359596_qp /))
  call stop_test()


  call start_test("Testing f_c0 from Genz...")
  call assert(f_c0_1p(x_1p, params),  &
                 0.9417645335842487_qp)
  y = f_c0_np(x_np, params)
  call assert(y(:), (/ 0.9417645335842487_qp, &
                                0.9048374180359595_qp, &
                                0.8693582353988057_qp /))
  call assert( f_c0_np(x_np(:,1:1), params), &
                                           (/0.9417645335842487_qp/))
  call stop_test()


  call start_test("Testing Fox's testfunction...")
  params%a = 4.0_qp
  params%u = 2.0_qp
  y = f_fox86_np(x_np, params)
  call assert(y(:), (/  0.64_qp, &
                                                    0.0_qp, &
                                                    0.32_qp /))
  call stop_test()


  call start_test("Testing f_atanassov03_f2_1p...")
  call assert(f_atanassov03_f2_1p(x_1p, params), &
                                    0.611314_qp, &
                                    spacing(0.611314_qp))
  call stop_test()


  call start_test("Testing f_atanassov03_f3_1p...")
  call assert(f_atanassov03_f3_1p(x_1p, params), &
                                    0.06_qp,                               &
                                    spacing(0.06_qp))
  call stop_test()


  call start_test("Testing Wang and Hickernell's testfunction...")
  params%a = 0.01_qp
  y = f_wang_hickernell(x_np, params)
  call assert(y(:), (/  0.6470934222135084_qp, &
                                                    0.01186158219782374_qp, &
                                                    0.3255563180080384_qp /))
  call stop_test()


  call start_test("Testing exact value for integral over unit-cube of Haselgrove's exponential function...")
  call assert( f_haselgrove61_exact(params, bounds), &
                                     0.7965995993_qp, &
                                     0.0000000001_qp)
  call stop_test()

  deallocate(x_np, y)
  deallocate(params%a, params%u, bounds%lb, bounds%ub)

  call start_test("Testing f_sine_power()...")
  allocate(params%k(2))
  params%k = 6
  call assert(f_sine_power_1p((/ 0.2_qp, 0.3_qp /), params), &
              1.2831435264880281E-7_qp)
  call stop_test()

  call start_test("Testing f_sine_power_prod_deriv(...)...")
  call assert(f_sine_power_deriv(0.4_qp, 3, 2, 6), 28.89208093268047_qp)
  call assert(f_sine_power_deriv(0.4_qp, 3, 2, 5), 18.21437945150671_qp)
  call assert(f_sine_power_deriv(0.6_qp, 4, 3, 6), 2150.402495485229_qp)
  call assert(f_sine_power_deriv(0.6_qp, 4, 3, 5), 1565.46058355222_qp)
  call assert(f_sine_power_deriv(0.0_qp, 4, 3, 5), 0.0_qp)
  call stop_test()

  call start_test("Testing f_sine_power_1d_exact()...")
  call assert(f_sine_power_1d_exact(1, 2, 0.0_qp, 1.0_qp), 0.5_qp)
  call assert(f_sine_power_1d_exact(2, 3, 0.0_qp, 1.0_qp), 0.0_qp)
  call assert(f_sine_power_1d_exact(3, 4, 0.0_qp, 1.0_qp), 3.0_qp/8)
  call assert(f_sine_power_1d_exact(4, 5, 0.0_qp, 1.0_qp), 0.0_qp)
  call assert(f_sine_power_1d_exact(5, 6, 0.0_qp, 1.0_qp), 5.0_qp/16)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_testfunctions
