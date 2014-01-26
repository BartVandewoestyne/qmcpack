! Module containing several system models.
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

module mod_system_model

  use numeric_kinds
  use mod_pdf

  private

  public :: next_gordon93salmond

  type, public :: system_model_parameters
    integer(kind=i4b) :: k
  end type system_model_parameters

contains

  
  ! A nonlinear one-dimensional example.  See [1].
  !
  subroutine next_gordon93salmond(params, x, noisy)
    type(system_model_parameters), intent(in) :: params
    real(kind=qp), intent(inout)              :: x
    logical, intent(in)                       :: noisy

    real(kind=qp) :: w
    real(kind=qp) :: variance

    variance = 10.0_qp

    ! Update the state
    x = 0.5_qp*x + 25*x/(1+x**2) + 8*cos(1.2_qp*(params%k-1))

    ! Add noise if necessary
    if (noisy) then

      call randn_box_muller(w, 0.0_qp, sqrt(variance))
      x = x + w

    end if

  end subroutine next_gordon93salmond


end module mod_system_model
