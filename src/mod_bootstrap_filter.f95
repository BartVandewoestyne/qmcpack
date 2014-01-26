! Module that implements the bootstrap-filer.
!
! References:
!
!  [1] Gordon, N. J., Salmond, D. J., Smith, A. F. M, `Novel approach to
!      nonlinear/non-Gaussian Bayesian state estimation', IEE Proceedings-F,
!      Vol. 140, No. 2, April 1993.

module mod_bootstrap_filter

  use numeric_kinds
  use mod_system_model
  use mod_measurement_model
  use mod_pdf
  use mod_sort

  private

  public :: initialize_bootstrap
  public :: step_bootstrap
  public :: free_bootstrap

  private :: predict
  private :: update
  private :: resample
  private :: get_measurement

  ! The number of samples
  integer(kind=i4b), private                        :: N

  ! The actual real value of x
!  real(kind=qp), private                            :: x_real

  ! An array containing the samples for the bootstrap filter
  real(kind=qp), dimension(:), allocatable, private :: x

  ! The weights accompanying each sample
  real(kind=qp), dimension(:), allocatable, private :: w

contains


  ! Initialize the bootstrap filter:
  !
  !   * set the current real value of x.
  !
  !   * set the samples of the bootstrap filter by sampling N times from
  !     the prior and assigning uniform weights.
  !
  !   * assign uniform weights.
  ! 
  ! Notes:
  !
  !   * currently, only sampling from a normal distribution is supported.  It
  !     would be nice if we could support any kind of distribution.
  !
  !   * this is also the place where we would use quasi-Monte Carlo points
  !     instead of pseudorandom points to generate the samples.
  !
  ! TODO:
  !   * try to see if we can get rid of this hardcoded stuff...
  !
  !   * make this more generic so we can use other priors than simply a normal
  !     distribution in 1 dimension.
  !
  subroutine initialize_bootstrap(num_samples)
    integer(kind=i4b), intent(in) :: num_samples

    integer(kind=i4b) :: i

    N = num_samples

    allocate(x(1:N))

    ! Get N samples from a normal distribution with mean 0 and variance 2
    ! TODO: check of matlab implementatie 2 of 5 gebruikt...
    do i=1,N
      call randn_box_muller(x(i), 0.0_qp, sqrt(2.0_qp))
    end do

    ! Set the weights to uniform weights
    allocate(w(1:N))
    w = 1.0_qp/N

  end subroutine initialize_bootstrap


  ! Update the internals of the bootstrap filter using the specified
  ! system model and measurement model.
  !
  subroutine step_bootstrap(x_real,                                  &
                            x_measured,                              &
                            x, w,                                    &
                            f_system_model, system_params,           &
                            f_measurement_model, measurement_params)
    real(kind=qp), intent(inout)             :: x_real
    real(kind=qp), intent(out)               :: x_measured
    real(kind=qp), dimension(:), intent(out) :: x
    real(kind=qp), dimension(:), intent(out) :: w
    interface
      subroutine f_system_model(fparams, x, noisy)
        use numeric_kinds
        use mod_system_model
        type(system_model_parameters), intent(in) :: fparams
        real(kind=qp), intent(inout)              :: x
        logical, intent(in)                       :: noisy
      end subroutine f_system_model
    end interface
    type(system_model_parameters), intent(in) :: system_params
    interface
      subroutine f_measurement_model(params, x, y, noisy)
        use numeric_kinds
        use mod_measurement_model
        type(measurement_model_parameters), intent(in) :: params
        real(kind=qp), intent(in)                      :: x
        real(kind=qp), intent(out)                     :: y
        logical, intent(in)                            :: noisy
      end subroutine f_measurement_model
    end interface
    type(measurement_model_parameters), intent(in) :: measurement_params


    ! Predict step, based on the system model (with noise)
    call predict(f_system_model, system_params, x)

    ! Update the real position and get the next measurement, based on this
    ! new real position.
    call f_system_model(system_params, x_real, .true.)
    call get_measurement(f_measurement_model, measurement_params, &
                         x_real, x_measured)

    ! Update step
    ! TODO: check this, because in some steps there appear weights that are
    !       much larger then one!!!
    call update(x_measured, f_measurement_model, measurement_params, x, w)

    ! Resample step: change the particles and the weights to their
    ! resampled values
    ! TODO:
    !   * it might be advantageous not to do the resampling every time,
    !     but only at certain moments in time, depending on the value
    !     of the weights (See PhD Klaas).
    call resample(x, w)

  end subroutine step_bootstrap


  ! Use the system model to execute the prediction step and do an initial
  ! update of the samples.
  !
  subroutine predict(f_system_model, params, x)
    interface
      subroutine f_system_model(fparams, x, noisy)
        use numeric_kinds
        use mod_system_model
        type(system_model_parameters), intent(in) :: fparams
        real(kind=qp), intent(inout)              :: x
        logical, intent(in)                       :: noisy
      end subroutine f_system_model
    end interface
    type(system_model_parameters), intent(in)  :: params
    real(kind=qp), dimension(:), intent(inout) :: x

    integer(kind=i4b) :: i

    do i=1,size(x)
      call f_system_model(params, x(i), .true.)
    end do

  end subroutine predict


  ! Given the current real position, return a measurement using the
  ! specified measurement model.
  !
  subroutine get_measurement(f_measurement_model, params, x, y)
    interface
      subroutine f_measurement_model(params, x, y, noisy)
        use numeric_kinds
        use mod_measurement_model
        type(measurement_model_parameters), intent(in) :: params
        real(kind=qp), intent(in)                      :: x
        real(kind=qp), intent(out)                     :: y
        logical, intent(in)                            :: noisy
      end subroutine f_measurement_model
    end interface
    type(measurement_model_parameters), intent(in) :: params
    real(kind=qp), intent(in)  :: x
    real(kind=qp), intent(out) :: y

    call f_measurement_model(params, x, y, .true.)

  end subroutine get_measurement


  ! On receipt of the measurement y, evaluate the likelihood of each prior
  ! sample and obtain a normalized weight for each sample.
  !
  ! TODO: remove the hardcoded stuff here!!!
  !
  subroutine update(y, f_measurement_model, measurement_params, x, q)
    real(kind=qp), intent(in) :: y
    interface
      subroutine f_measurement_model(params, x, y, noisy)
        use numeric_kinds
        use mod_measurement_model
        type(measurement_model_parameters), intent(in) :: params
        real(kind=qp), intent(in)                      :: x
        real(kind=qp), intent(out)                     :: y
        logical, intent(in)                            :: noisy
      end subroutine f_measurement_model
    end interface
    type(measurement_model_parameters), intent(in) :: measurement_params
    real(kind=qp), dimension(:), intent(in)  :: x
    real(kind=qp), dimension(:), intent(out) :: q

    integer(kind=i4b)                 :: i
    real(kind=qp), dimension(size(x)) :: m
    real(kind=qp)                     :: mysum
    real(kind=qp)                     :: variance

    variance = 1.0_qp

    ! Apply the measurement model to the x-values (without noise!).
    do i=1,size(x)
      call f_measurement_model(measurement_params, x(i), m(i), .false.)
    end do

    ! Update the weights and make them normalized.
    mysum = 0.0_qp
    do i=1,size(m)
      mysum = mysum + norm_pdf(y, m(i), variance)
    end do
    do i=1,size(m)
      q(i) = norm_pdf(y, m(i), variance)/mysum
    end do

    if (any(q>1.0_qp)) then
        print *, "ERROR: found a weight > 1!!!!"
    end if

    !mysum = 0.0_qp
    !do i=1,size(m)
    !  mysum = mysum + exp(-0.5_qp*(y-m(i))**2)
    !end do
    !do i=1,size(m)
    !  q(i) = exp(-0.5_qp*(y-m(i))**2)/mysum
    !end do

  end subroutine update


  ! Resample from the discrete distribution to generate new samples.
  ! (See [1], page 109 for the description of the method).
  !
  subroutine resample(x, w)
    real(kind=qp), dimension(:), intent(inout) :: x, w

    real(kind=qp)                         :: cumsum
    real(kind=qp), dimension(size(x))     :: u
    real(kind=qp), dimension(size(x))     :: x_sorted, w_sorted
    integer(kind=i4b)                     :: i, j
    integer(kind=i4b), dimension(size(x)) :: indices
    integer(kind=i4b), dimension(size(x)) :: dummy

    ! Sort the particles (and re-organize the weights accordingly)
    x_sorted = x
    w_sorted = w
    indices = (/ ( i, i=1,size(x) ) /)
    call quicksort(x_sorted, indices)
    w_sorted = w(indices)
    !print *, "x_sorted = ", x_sorted
    !print *, "w_sorted = ", w_sorted

    ! Generate N uniform random numbers and sort them
    call random_number(u)
    call quicksort(u, dummy)
    !print *, "u = ", u

    ! Now determine the new particles and samples, using
    ! the cumulative distribution (See [1], page 109).
    cumsum = 0.0_qp
    j = 0
    do i=1,size(x)
      do
        if (cumsum > u(i)) then
          exit
        end if
        j = j + 1
        cumsum = cumsum + w_sorted(j)
      end do

      ! The new particles
      x(i) = x_sorted(j)

      ! The new weights (by my own invented method...)
      if (i==1) then
        w(i) = u(i)
      else
        w(i) = u(i) - u(i-1)
      end if

    end do

  end subroutine resample


  ! Free up any memory used by the bootstrap filter
  !
  subroutine free_bootstrap()
    
    if (allocated(w)) then
      deallocate(w)
    end if
    if (allocated(x)) then
      deallocate(x)
    end if

  end subroutine free_bootstrap

end module mod_bootstrap_filter
