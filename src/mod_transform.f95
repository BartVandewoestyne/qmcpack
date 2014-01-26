! Module that implements several transformations that can be applied to
! a pointset.
!
! References:
!
!   [1] `Obtaining O(N^{-2+\varepsilon}) Convergence for Lattice Quadrature
!       Rules', Fred J. Hickernell, Monte Carlo and Quasi-Monte Carlo
!       Methods 2000, pages 274--289, 2002.
!
! TODO:
!   * Add the transform from page 18 of Halve's technical report, which is
!     a generalization of the Baker's transform.

module mod_transform

  use numeric_kinds

  private

  public :: baker_transform

  interface baker_transform
    module procedure baker_transform_point, baker_transform_pointset
  end interface baker_transform


  private :: baker_transform_point
  private :: baker_transform_pointset

contains

  ! Apply the Baker's transform (also called `Baker's map') to each dimension
  ! of the point x.  The Baker's transform is defined as
  !
  !             phi(t) = 1 - abs(2*t-1)
  !
  subroutine baker_transform_point(x)
    real(kind=qp), dimension(:), intent(inout) :: x

    x = 1 - abs(2*x-1)

  end subroutine baker_transform_point


  ! Apply the Baker's transform (also called `Baker's map') to each dimension
  ! of the pointset x.  The Baker's transform is defined as
  !
  !             phi(t) = 1 - abs(2*t-1)
  !
  subroutine baker_transform_pointset(x)
    real(kind=qp), dimension(:,:), intent(inout) :: x

    x = 1 - abs(2*x-1)

  end subroutine baker_transform_pointset


end module mod_transform
