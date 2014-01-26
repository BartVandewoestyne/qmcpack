module mod_integration

  use numeric_kinds

  private

  public :: create_integration_result
  public :: set_unit_hypercube
  public :: allocate_integration_bounds


  type, public :: integration_result
    real(kind=qp)     :: i       ! the computed value of the integral
    real(kind=qp)     :: abs_err ! the absolute error
    real(kind=qp)     :: rel_err ! the relative error
    integer(kind=i4b) :: n       ! the value of n for which these results apply
  end type integration_result


  type, public :: integrationbounds
    real(kind=qp), dimension(:), pointer :: lb ! vector with lower bounds
    real(kind=qp), dimension(:), pointer :: ub ! vector with upper bounds
  end type integrationbounds


  contains

    ! Using the information of the estimate and the exact value, calculate
    ! an integration result containing
    !
    !   * the estimate for integral value
    !   * the absolute error for that estimate
    !   * the relative error for that estimate
    !   * the value of n for which this result applies
    !
    subroutine create_integration_result(myresult, i_estimate, i_exact, n)
      type(integration_result), intent(out) :: myresult
      real(kind=qp), intent(in)             :: i_estimate
      real(kind=qp), intent(in)             :: i_exact
      integer(kind=i4b), intent(in)         :: n

      myresult%i = i_estimate
      myresult%abs_err = i_exact-i_estimate
      myresult%rel_err = (i_exact-i_estimate)/i_exact
      myresult%n = n

    end subroutine create_integration_result


    ! Set the bounds to be bounds for the unit hypercube
    !
    subroutine set_unit_hypercube(bounds)
      type(integrationbounds), intent(inout) :: bounds

      bounds%lb = 0.0_qp
      bounds%ub = 1.0_qp

    end subroutine set_unit_hypercube


    ! Allocate memory for the integration bounds
    !
    subroutine allocate_integration_bounds(bounds, s)
      type(integrationbounds), intent(inout) :: bounds
      integer(kind=i4b), intent(in)          :: s

      allocate(bounds%lb(s))
      allocate(bounds%ub(s))

    end subroutine allocate_integration_bounds

end module mod_integration
