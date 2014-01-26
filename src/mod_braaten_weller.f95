! Module that implements the testfunction that is used in the paper of
! Braaten and Weller.
!
! References:
!   [1] 'An improved low-discrepancy sequence for multidimensional quasi-Monte
!       Carlo integration', Braaten, Eric and Weller, George, Journal of
!       Computational Physics, vol. 33, 1979, pages 249--258.

module mod_braaten_weller

  use numeric_kinds

  private

  public :: f_gaussian_bw

  contains

    ! the rows for xt are the different dimensions!!!
    function f_gaussian_bw(xt) result (y)
      real(kind=qp), dimension(:,:), intent(in) :: xt
      real(kind=qp), dimension(size(xt, dim=2)) :: y

      y = product(exp(-2*(xt-0.5)**2), dim=1)
    
    end function f_gaussian_bw


    ! The exact value of the integral in [0,1]
    ! TODO: implement error function so we can get this exact value
    ! References:
    ! Masters thesis R. Schuerer
    !function f_gaussian_bw_exact(s) result (f)
    !  a = sqrt(2.0_qp)
    !  u = 0.5_qp
    !  t = 1
    !  r = 0
    !end function

end module mod_braaten_weller
