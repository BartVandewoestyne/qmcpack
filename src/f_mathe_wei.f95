! Module that implements the testfunctions from the article of Mathe and Wei.
!
! References:
!
!   [1] `Quasi-Monte Carlo integration over R^d', Peter Mathe and Gang Wei,
!        Mathematics of Computation, Volume 73, Number 246, Pages 827-841, 2004.

module f_mathe_wei

  use numeric_kinds
  use mod_function
  use mod_constants

  private

  public :: f_mathe_wei01_np
  public :: f_mathe_wei02_np
  public :: rho_mathe_wei01_np
  public :: rho_mathe_wei02_np

contains

  ! The first function from Mathe and Wey's paper.
  ! (See page 837-838 in [1]).
  !
  function f_mathe_wei01_np(x, params) result(res)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: res

    real(kind=qp), dimension(2,2) :: inv_sigma
    integer(kind=i4b)             :: i

    inv_sigma(1,:) = (/ 0.25_qp, 0.0_qp /)
    inv_sigma(2,:) = (/ 0.0_qp,  1.0_qp /)

    do i=1,size(x,dim=2)
      res(i) = 1.0_qp/(1+0.25_qp*x(1,i)*x(1,i)+x(2,i)*x(2,i))
    end do

  end function f_mathe_wei01_np


  ! The first probability density function from Mathe and Wey's paper
  ! (See page 837-838 in [1]).
  !
  function rho_mathe_wei01_np(x, params) result(rho)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: rho

    real(kind=qp), dimension(2,2) :: inv_sigma
    integer(kind=i4b)             :: i

    inv_sigma(1,:) = (/ 0.25_qp, 0.0_qp /)
    inv_sigma(2,:) = (/ 0.0_qp,  1.0_qp /)

    do i=1,size(x,2)
      rho(i) = 1.0_qp/(4*pi)*(1+0.5_qp*(0.25_qp*x(1,i)*x(1,i)+x(2,i)*x(2,i)))**(-2.0_qp)
    end do

  end function rho_mathe_wei01_np


  ! The second function from Mathe and Wey's paper.
  ! (See page 837-838 in [1])
  !
  function f_mathe_wei02_np(x, params) result(res)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: res

    real(kind=qp), dimension(2,2) :: sigma
    real(kind=qp), dimension(2,2) :: inv_sigma
    real(kind=qp)                 :: det_sigma
    real(kind=qp), dimension(2)   :: rs
    integer(kind=i4b)             :: i

    sigma(1,:) = (/ 4.0_qp, 1.9_qp /)
    sigma(2,:) = (/ 1.9_qp, 1.0_qp /)

    ! Calculate determinant of sigma
    det_sigma = sigma(1,1)*sigma(2,2)-sigma(2,1)*sigma(1,2)

    ! Calculate inverse of sigma
    inv_sigma(1,:) = (/  sigma(2,2), -sigma(1,2) /)/det_sigma
    inv_sigma(2,:) = (/ -sigma(2,1),  sigma(1,1) /)/det_sigma

    do i=1,size(x,2)
      rs = matmul(x(1:2,i), inv_sigma)
      res(i) = 1/( 1+dot_product(rs, x(1:2,i)) )
    end do

  end function f_mathe_wei02_np


  ! The second probability density function from Mathe and Wey's paper.
  ! (See page 837-838 in [1])
  !
  function rho_mathe_wei02_np(x, params) result(rho)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: rho

    real(kind=qp), dimension(2,2) :: sigma
    real(kind=qp), dimension(2,2) :: inv_sigma
    real(kind=qp)                 :: det_sigma
    integer(kind=i4b)             :: s
    real(kind=qp), dimension(2)   :: rs
    integer(kind=i4b)             :: i

    s = 4

    sigma(1,:) = (/ 4.0_qp, 1.9_qp /)
    sigma(2,:) = (/ 1.9_qp, 1.0_qp /)

    ! Calculate determinant of sigma
    det_sigma = sigma(1,1)*sigma(2,2)-sigma(2,1)*sigma(1,2)

    ! Calculate inverse of sigma
    inv_sigma(1,:) = (/sigma(2,2), -sigma(1,2) /)/det_sigma
    inv_sigma(2,:) = (/-sigma(2,1), sigma(1,1) /)/det_sigma

    do i=1,size(x,2)
      rs = matmul(x(1:2,i), inv_sigma)
      rho(i) = 1/sqrt(det_sigma)*(s-2.0_qp)/(4*pi)*(1+0.5*dot_product(rs, x(1:2,i)))**(-0.5_qp*s)
    end do

  end function rho_mathe_wei02_np


end module f_mathe_wei
