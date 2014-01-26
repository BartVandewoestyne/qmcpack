! Module containing several measurement models.
!
! References:
!
!  [1] Gordon, N. J., Salmond, D. J., Smith, A. F. M, `Novel approach to
!      nonlinear/non-Gaussian Bayesian state estimation', IEE Proceedings-F,
!      Vol. 140, No. 2, April 1993.
!
! TODO:
!   * Check what the best approach is to make this module as generic as
!     possible.  Should we create elemental next_xxx(x) methods that can
!     operate on arrays of different dimensions, or should we see the
!     next_xxx(x) methods as basic building blocks and only allow real
!     values x as their arguments and let the caller do the looping if we
!     have an array of a certain dimension.

module mod_measurement_model

  use numeric_kinds
  use mod_pdf

!  implicit none

  private

  public :: measure_gordon93salmond

  type, public :: measurement_model_parameters
    integer(kind=i4b) :: k
  end type measurement_model_parameters

contains

  
  ! A nonlinear one-dimensional example.  See [1].
  !
  subroutine measure_gordon93salmond(params, x, y, noisy)
    type(measurement_model_parameters), intent(in) :: params
    real(kind=qp), intent(in)                      :: x
    real(kind=qp), intent(out)                     :: y
    logical, intent(in)                            :: noisy

    real(kind=qp)     :: v
    real(kind=qp)     :: variance

    variance = 1.0_qp

    ! Update the state
    y = x*x/20

    if (noisy) then

      call randn_box_muller(v, 0.0_qp, sqrt(variance))
      y = y + v

    end if

  end subroutine measure_gordon93salmond


end module mod_measurement_model
