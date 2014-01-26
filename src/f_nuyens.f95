! Module that implements some of Dirk's simple and quick testfunctions
! (obtained via private communication).

module f_nuyens

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_constants

  private

  public :: f_nuyens01_1p
  public :: f_nuyens01_exact

  contains

    !                s
    !             --------'
    !            '  |  |
    !    f(x) :=    |  |    (1 + sin(x[i]))
    !               |  |
    !               |  |
    !              i = 1
    !
    function f_nuyens01_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(1+sin(x))

    end function f_nuyens01_1p

    function f_nuyens01_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( ub - lb + cos(lb) - cos(ub) )

    end function f_nuyens01_exact


end module f_nuyens
