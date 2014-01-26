! Module to transform points from the unit cube to the surface of a sphere.
!
! References:
!
!   [1] `Number-theoretic Methods in Statistics', Wang, Y. and Fang, K.-T.,
!       1994, ISBN 0-412-46520-5.
!
!   [2] `Sampling with Hammersley and Halton Points', Wong, Tien-Tsin and
!        Luk, Wai-Shing and Heng, Pheng-Ann.

module mod_sphere

  use numeric_kinds
  use mod_constants

  private

  public :: transform_to_sphere3D_fang_wang
  public :: transform_to_sphere3D_wong_luk

contains

  ! Transform the two-dimensional point c in the two-dimensional unit cube
  ! to a point x on the sphere, using Fang and Wang's transformation method.
  ! (See [1], page 50).
  !
  subroutine transform_to_sphere3D_fang_wang(c, x)
    real(kind=qp), dimension(:), intent(in)  :: c
    real(kind=qp), dimension(:), intent(out) :: x

    x(1) = 1-2*c(1)
    x(2) = 2*sqrt(c(1)*(1-c(1)))*cos(2*pi*c(2))
    x(3) = 2*sqrt(c(1)*(1-c(1)))*sin(2*pi*c(2))

  end subroutine transform_to_sphere3D_fang_wang

  
  ! Transform the two-dimensional point c in the two-dimensional unit cube
  ! to a point x on the sphere, using the transformation method as
  ! described in [2].  Note that this is about the same method as in
  ! Fang and Wang their book, with the following interchanges:
  !
  !     c(1) <--> c(2)
  !     x(1) <--> x(3)
  !     x(2) <--> x(1)
  !     x(3) <--> x(2)
  !
  subroutine transform_to_sphere3D_wong_luk(c, x)
    real(kind=qp), dimension(:), intent(in)  :: c
    real(kind=qp), dimension(:), intent(out) :: x

    x(1) = 2*sqrt(c(2)*(1-c(2)))*cos(2*pi*c(1))
    x(2) = 2*sqrt(c(2)*(1-c(2)))*sin(2*pi*c(1))
    x(3) = 2*c(2)-1
  end subroutine transform_to_sphere3D_wong_luk

end module mod_sphere
