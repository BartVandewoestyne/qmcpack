! Module that implements weight-functions that can be used in numerical
! integration routines
!
! References:
!
!   [1] Sugihara, Masaaki and Murota, Kazuo, `A Note on Haselgrove's Method
!       for Numerical Integration', Mathematics of Computation, volume 39,
!       number 160, October 1982, pages 549-554.

module mod_weightfunction

  use numeric_kinds
  use mod_special_functions
  use mod_constants

  private

  public :: imt
  public :: sag_szekeres
  public :: sugihara_murota_gaussian


contains


  ! The IMT-transformation function, normalized so
  ! that it integrates to 1 in the interval [0,1].
  !
  function imt(x) result (y)
    real(kind=qp), intent(in) :: x
    real(kind=qp)             :: y

    real(kind=qp) :: c

    ! The normalizing constant
    c = 0.0070298584066096562392412705303540_qp

    y = exp(-1.0_qp/(x*(1.0_qp-x)))/c

  end function imt


  ! The TANH-transformation function from Sag and Szekeres,
  ! normalized so that it integrates to 1 in the interval [0,1].
  !
  function sag_szekeres(x) result (y)
    real(kind=qp), intent(in) :: x
    real(kind=qp)             :: y

    y = tanh(-0.5_qp*(1/x-1/(1-x))) + 1.0_qp

  end function sag_szekeres


  ! The Gaussian-type weightfunction from [1].  Note that in [1],
  ! the weightfunction does not integrate to 1 in [0,1] yet.  We
  ! normalize it here.
  !
  function sugihara_murota_gaussian(x, B, N, d) result (y)
    real(kind=qp), intent(in)     :: x
    real(kind=qp), intent(in)     :: B
    integer(kind=i4b), intent(in) :: N
    integer(kind=i4b), intent(in) :: d
    real(kind=qp)                 :: y

    real(kind=qp) :: C
    real(kind=qp) :: alpha

    alpha = 1.0_qp/(d+1)

    ! Calculate the normalizing constant that makes the
    ! weightfunction integrate to 1 over [0,1].
    C = erf(0.5_qp*sqrt(B*(N**alpha)))/N

    y = 1/C * sqrt(B/(Pi*(N**(2-alpha))))*exp(-B*(N**alpha)*((x-0.5_qp)**2))

  end function sugihara_murota_gaussian


end module mod_weightfunction
