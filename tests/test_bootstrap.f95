program test_bootstrap

  use qmcpack

  real(kind=qp)                            :: x_real, y_measured
  type(system_model_parameters)            :: sys_param
  type(measurement_model_parameters)       :: meas_param
  real(kind=qp), dimension(:), allocatable :: x, w
  integer(kind=i4b)                        :: i
  integer(kind=i4b)                        :: my_unit1, my_unit2
  integer                                  :: ios1, ios2
  character(len=100)                       :: fmt_string

  integer(kind=i4b), parameter :: nb_samples=100, nb_timesteps=50
  real(kind=qp), dimension(nb_samples, nb_timesteps) :: samples, weights
  
  call show_test_header("MOD_BOOTSTRAP")
  call init_assert()

  call start_test("Testing bootstrap filter by writing 50 real and measured positions to file bootstrap.dat...")

  call get_unit(my_unit1)
  open(unit=my_unit1, file="bootstrap.dat", iostat=ios1, &
       status="replace", access="sequential", action="write")
  if (ios1 /= 0) then
    write(unit=*, fmt="(A)") "ERROR while opening bootstrap data file!"
    call increase_nb_assert_errors(1)
  else

    x_real = 0.1_qp
    sys_param%k = 1

    ! TODO: check if this should be 0
    meas_param%k = 0

    allocate(x(nb_samples), w(nb_samples))
    call initialize_bootstrap(nb_samples)
    do i=1,nb_timesteps
      sys_param%k = sys_param%k + 1
      call step_bootstrap(x_real, y_measured,                  &
                          x, w,                                &
                          next_gordon93salmond, sys_param,     &
                          measure_gordon93salmond, meas_param)

      write(unit=my_unit1, fmt="(I10, 3(ES25.15e4))") i, x_real, y_measured, sum(x*w)
      samples(:,i) = x
      weights(:,i) = w

    end do
    deallocate(x, w)
  end if

  close(unit=my_unit1)

  call free_bootstrap()

  if (any(weights>1.0_qp)) then
    print *, "ERROR: some weights are larger than 1!!!!"
  end if

  ! Write the samples and weights to a datafile
  call get_unit(my_unit2)
  open(unit=my_unit2, file="bootstrap_samples_weights.dat", iostat=ios2, recl=2*nb_timesteps*25, &
       status="replace", access="sequential", action="write")
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", 2*nb_timesteps, "ES25.15e4)"
  do i=1,nb_samples
    write(unit=my_unit2, fmt=fmt_string) samples(i,:), weights(i,:)
  end do
  close(unit=my_unit2)

  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_bootstrap
