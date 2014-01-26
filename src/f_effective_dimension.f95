! Module that implements testfunctions for which the effective
! dimension can be changed by setting certain values for the parameters.
!
! References:
!   [1] `The effective dimension and quasi-Monte Carlo integration', Wang,
!   Xiaoqun, Fang, Kai-Tai, Journal of Complexity, Vol. 19, pages 101-124, 2003.
!
module f_effective_dimension

  use numeric_kinds
  use mod_function
  use mod_integration

  private
 
  public :: f_wang03fang_1p
  public :: f_wang03fang_exact

contains

 !             d
 !          --------'
 !         '  |  |    | 4 x[k] - 2 | + a[k]
 ! f(x) =     |  |    ---------------------
 !            |  |          1 + a[k]
 !            |  |
 !           k = 1
 !      
  function f_wang03fang_1p(x, params) result(res)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: res

    res = product((abs(4*x-2) + params%a)/(1+params%a))

  end function f_wang03fang_1p

  ! TODO:
  !   Currently, the exact value is only exact when integration is
  !   done over the unit cube!  WATCH OUT!
  !
  function f_wang03fang_exact(params, bounds) result (i)
    type(functionparams), intent(in)    :: params
    type(integrationbounds), intent(in) :: bounds

    real(kind=qp) :: i

    i = 1.0_qp

  end function f_wang03fang_exact


end module f_effective_dimension
