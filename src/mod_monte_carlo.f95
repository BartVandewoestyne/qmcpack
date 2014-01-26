! Module to do a simple Monte Carlo simulation.
!
! Note that doing it this way is limited, because the x-abscissa should already
! contain the points to evaluate the function at... and this might use a lot
! of memory if you want to do it for a lot of points...

module mod_monte_carlo

  use numeric_kinds
  use mod_function
  use mod_integration

  private

  public :: sim_mc

  contains

    ! Notes:
    !   * the rows of x are the dimensions!!!
    !   * func should be a function that can act on arrays of rank 2!!!
    !
    subroutine sim_mc(func, x, params, exact, myresult)
      interface
        function func(x, fparams) result (y)
          use numeric_kinds
          use mod_function
          real(kind=qp), dimension(:,:), intent(in) :: x
          type(functionparams), intent(in)          :: fparams
          real(kind=qp), dimension(size(x,dim=2))   :: y
        end function func
      end interface
      real(kind=qp), dimension(:,:), intent(in)           :: x
      type(functionparams), intent(in)                    :: params
      real(kind=qp), intent(in)                           :: exact
      type(integration_result), dimension(:), intent(out) :: myresult

      real(kind=qp), dimension(size(x, dim=2)) :: fvalues
      real(kind=qp)                            :: currentsum
      integer(kind=i4b)                        :: N
      integer(kind=i4b)                        :: i


      N = size(x, dim=2)

      fvalues = func(x, params)

      currentsum = 0.0_qp
      do i=1,N
        currentsum = currentsum + fvalues(i)
        call create_integration_result(myresult(i), currentsum/i, exact, i)
      end do

  end subroutine sim_mc

end module mod_monte_carlo
