! Module that implements rank-1 lattice rules.
!
! References:
!
!   [1] `Lattice Methods for Multiple Integration', Sloan, I. H. and Joe, S.,
!       Oxford Science Publications, 1994, ISBN 0 19 853472 8.

module mod_rank1_lattice

  use numeric_kinds
  use mod_utilities

  private

  public :: create_korobov_vector
  public :: rank1_lattice

contains


  ! Create a `reduced' Korobov vector that is suitable for working in the
  ! limited precision of a computer.  Use the property that
  !
  ! frac_part(i*a/N) = frac_part(i*(a1+a2)/N)
  !                  = frac_part(i*a1/N + i*a2/N)
  !                  ... and because a1 is a multiple of N we have ...
  !                  = frac_part(i*a2/N)
  !
  function create_korobov_vector(N, a, s) result (vec)

    integer(kind=i4b), intent(in)   :: N
    integer(kind=i4b), intent(in)   :: a
    integer(kind=i4b), intent(in)   :: s
    integer(kind=i4b), dimension(s) :: vec

    integer(kind=i4b) :: i

    vec(1) = 1
    vec(2) = modulo(a, N)

    do i=3,s
      vec(i) = modulo(vec(i-1)*vec(2), N)
    end do

  end function create_korobov_vector


  ! Construct a rank-1 lattice of order N and generating vector z.
  ! The columns of x represent the different dimensions.
  ! (See [1], page 69).
  !
  subroutine rank1_lattice(N, z, x)

    integer(kind=i4b), intent(in)               :: N
    integer(kind=i4b), dimension(:), intent(in) :: z
    real(kind=qp), dimension(:,:), intent(out)  :: x

    integer(kind=i4b) :: i

    ! TODO: Check if the vector z has no factor in common with N.

    ! Calculate the lattice
    x = frac_part( spread(z, 1, N) &
                   *spread( (/ (i/real(N, kind=qp), i=0,N-1) /), 2, size(z) ) )

  end subroutine rank1_lattice


end module mod_rank1_lattice
