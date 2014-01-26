! Module that implements testfunctions of unbounded variation.
!
module f_unbounded_variation

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_constants

  private
 
  public :: f_circle_1p
  public :: f_circle_exact

contains


  function f_circle_1p(x, params) result(res)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: res

    if (size(x) /= 2) then
      print *, "ERROR: f_circle is only to be used in 2 dimensions!"
    else 
      if ( (x(1)-0.5_qp)**2 + (x(2)-0.5_qp)**2 <= 0.25 ) then
        res = 1
      else
        res = 0
      end if
    end if

  end function f_circle_1p

  ! TODO:
  !   Currently, the exact value is only exact when integration is
  !   done over the unit cube!  WATCH OUT!
  !
  function f_circle_exact(params, bounds) result (i)
    type(functionparams), intent(in)    :: params
    type(integrationbounds), intent(in) :: bounds

    real(kind=qp) :: i

    i = 0.25*pi

  end function f_circle_exact


end module f_unbounded_variation
