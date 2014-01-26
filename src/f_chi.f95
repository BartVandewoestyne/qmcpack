! Module that implements Hongmei's testfunctions.
!
module f_chi

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_constants

  private
 
  public :: f_sin_1p
  public :: f_sin_exact

contains


  ! Simple sinus function:
  !              s
  !           --------'
  !          '  |  |
  !  f(x) :=    |  |    (1/2 Pi sin(pi x[i]))
  !             |  |
  !             |  |
  !            i = 1
  !
  function f_sin_1p(x, params) result(res)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: res

    res = product( 0.5_qp*pi*sin(pi*x) )

  end function f_sin_1p

  function f_sin_exact(params, bounds) result(i)
    type(functionparams), intent(in)    :: params
    type(integrationbounds), intent(in) :: bounds
    real(kind=qp)                       :: i

    real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

    lb = bounds%lb
    ub = bounds%ub

    i = product( 0.5_qp*(cos(pi*lb)-cos(pi*ub)) )

  end function f_sin_exact

end module f_chi
