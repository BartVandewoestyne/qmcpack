! Module that implements some of Tony Warnock's testfunctions (obtained
! via private communication).
!

module f_warnock

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_constants

  private

  public :: f_warnock01_1p
  public :: f_warnock01_exact
  public :: f_warnock02_1p
  public :: f_warnock02_exact
  public :: f_warnock03_1p
  public :: f_warnock03_exact
  public :: f_warnock02_bernoulli3_1p
  public :: f_warnock02_bernoulli5_1p

  contains

    !                   s
    !                --------'
    !               '  |  |    2 exp(2 x[i]) - 1
    !       f(x) :=    |  |    -----------------
    !                  |  |       exp(2) - 1
    !                  |  |
    !                 i = 1
    !
    function f_warnock01_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product((2*exp(2*x)-1)/(exp(2.0_qp)-1))

    end function f_warnock01_1p

    function f_warnock01_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( (exp(2*ub) - exp(2*lb) + lb - ub)/(exp(2.0_qp)-1) )

    end function f_warnock01_exact


    !                 s
    !              --------'
    !             '  |  |           5                 2
    !     f(x) :=    |  |    (3 x[i]  + sin(6 pi x[i]) )
    !                |  |
    !                |  |
    !               i = 1
    !
    function f_warnock02_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( 3.0_qp*x**5 + (sin(6*pi*x))**2 )

    end function f_warnock02_1p

    function f_warnock02_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( 0.5*(ub**6-lb**6+ub-lb) &
                    - sin(12*pi*ub)/(24*pi) + sin(12*pi*lb)/(24*pi) )

    end function f_warnock02_exact


    !                      s
    !                   --------'
    !                  '  |  |            29
    !          f(x) :=    |  |    (30 x[i]  )
    !                     |  |
    !                     |  |
    !                    i = 1
    !
    function f_warnock03_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( 30*x**29 )

    end function f_warnock03_1p

    function f_warnock03_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product(ub**30-lb**30)

    end function f_warnock03_exact


    ! The `Bernoulli-method'-periodized version of
    !
    !                 s
    !              --------'
    !             '  |  |           5                 2
    !     f(x) :=    |  |    (3 x[i]  + sin(6 pi x[i]) )
    !                |  |
    !                |  |
    !               i = 1
    !
    ! for alpha = 3.
    !
    function f_warnock02_bernoulli3_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( (sin(6*pi*x))**2 &
                        + (x*(x*(x*(3*x*x-10)+7.5_qp)-0.5_qp)+0.25_qp) )

    end function f_warnock02_bernoulli3_1p


    ! The `Bernoulli-method'-periodized version of
    !
    !                 s
    !              --------'
    !             '  |  |           5                 2
    !     f(x) :=    |  |    (3 x[i]  + sin(6 pi x[i]) )
    !                |  |
    !                |  |
    !               i = 1
    !
    ! for alpha >= 5.
    !
    function f_warnock02_bernoulli5_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( 0.5_qp*(2*sin(6*pi*x)**2+1) )

    end function f_warnock02_bernoulli5_1p


end module f_warnock
