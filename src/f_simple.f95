! Module that implements some simple testfunctions.
!
module f_simple

  use numeric_kinds
  use mod_function

  private
 
  public :: f_sum_of_squares_np
  public :: f_one_np

contains


  ! Sum of the squares of the coordinates.
  !
  function f_sum_of_squares_np(x, params) result(res)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: res

    res = sum(x*x, dim=1)

  end function f_sum_of_squares_np


  ! Function that is 1 everywhere.
  !
  function f_one_np(x, params) result (res)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: res

    res = 1

  end function f_one_np


end module f_simple
