! Module that implements testfunctions from:
!
! References:
!
!  [1] Gajda, Piotr, Li, Youming, Plaskota, Leszek, Wasilkowski, Grzegorz W.,
!      `A Monte Carlo Algorithm for Weighted Integration over R^d', Mathematics
!      of Computation, Vol. 73, Nb. 246, pages 813-825, 2003.

module f_gajda

  use numeric_kinds
  use mod_function
  use mod_constants

  private
 
  public :: rho_gajda01_np
  public :: rho_gajda02_np
  public :: rho_gajda03_np
  public :: f_gajda01_np
  public :: f_gajda02_np
  public :: f_gajda03_np

contains


  ! rho_1 = exp(-||x||^2)
  !
  function rho_gajda01_np(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    y = exp(-sum(x**2, dim=1))

  end function rho_gajda01_np


  ! rho_2 = (1-||x||)/(1-||x||^{d+3})
  !
  function rho_gajda02_np(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    integer(kind=i4b) :: d

    d = size(x, dim=1)

    y = (1-sqrt(sum(x**2, dim=1)))/(1-(sqrt(sum(x**2, dim=1)))**(d+3.0_qp))

  end function rho_gajda02_np


  ! Normalized Gaussian weight (for mortgage-backed securities)
  !
  ! rho_3 = (2*pi)^(-d/2)*exp(-||x||^2/2)
  !
  function rho_gajda03_np(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    integer(kind=i4b) :: d

    d = size(x, dim=1)

    y = (2*pi)**(-0.5_qp*d)*exp(-0.5*sum(x*x, dim=1))

  end function rho_gajda03_np


  ! f_1 = cos(||x||)
  !
  function f_gajda01_np(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    y = cos(sqrt(sum(x**2, dim=1)))
    
  end function f_gajda01_np


  ! f_2 = sum_{k=1:d} 1/(1+sqrt(|x_k|))
  !
  function f_gajda02_np(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    y = sum(1/(1+sqrt(abs(x))), dim=1)
    
  end function f_gajda02_np


  ! f_3 = sum_{k=1:d}(|x_k|)
  !
  function f_gajda03_np(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    y = sum(abs(x), dim=1)
    
  end function f_gajda03_np


end module f_gajda
