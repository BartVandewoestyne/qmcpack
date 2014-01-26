! Module that implements methods to calculate several types of distances between
! points.
!
! References:
!
!   [1] `Extremal systems of points and numerical integration on the sphere',
!       Ian H. Sloan and Robert S. Womersley, Advances in Computational
!       Mathematics, Vol. 21, pages 107--125, 2004.
!
module mod_distance

  use numeric_kinds

  private

  public :: euclidian_distance
  public :: geodesic_distance

contains

  ! Return the Euclidian distance between two multi-dimensional point x and y.
  !
  pure function euclidian_distance(x, y) result (d)
    real(kind=qp), dimension(:), intent(in) :: x, y
    real(kind=qp)                           :: d

    d = sqrt(sum((x-y)**2))

  end function euclidian_distance


  ! Return the geodesic/geodetic distance between points x and y on
  ! a multi-dimensional *unit* sphere (See [1], page 120).
  !
  ! Note: this code is only valid if the points x and y are on the unit sphere!
  !
  pure function geodesic_distance(x, y) result (d)
    real(kind=qp), dimension(:), intent(in) :: x, y
    real(kind=qp)                           :: d

    d = acos(sum(x*y))

  end function geodesic_distance

end module mod_distance
