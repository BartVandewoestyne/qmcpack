! Module that implements the random-shift technique introduced by Cranley and
! Patterson.
!
! References:
!
!   [1] `Randomization of number theoretic methods for multiple integration',
!       Cranley, R. and Patterson, T. N. L., SIAM Journal on Numberical
!       Analysis, vol. 13, nr. 6, december 1976, pages 904--914.

module cranley_patterson

  use numeric_kinds
  use mod_utilities

  private

  public :: init_random_shift
  public :: random_shift
  public :: free_cranley_patterson

  interface random_shift
    module procedure random_shift_point, random_shift_pointset
  end interface random_shift


  private :: random_shift_point
  private :: random_shift_pointset

  real(kind=qp), dimension(:), allocatable, private :: shiftvector

contains

  ! Initialize the random shift vector to be used in the
  ! randomization operation.
  !
  subroutine init_random_shift(s)
    integer(kind=i4b), intent(in) :: s

    if (allocated(shiftvector)) then
      if (size(shiftvector) /= s) then
        deallocate(shiftvector)
      end if
    end if
    if (.not. allocated(shiftvector)) then
      allocate(shiftvector(s))
    end if

    call random_number(shiftvector)
    
  end subroutine init_random_shift


  ! Randomly shift a single point x with the current shiftvector
  !
  subroutine random_shift_point(x)
    real(kind=qp), dimension(:), intent(inout) :: x

    if (.not. allocated(shiftvector)) then
      allocate(shiftvector(size(x)))
      call random_number(shiftvector)
    end if

    x = frac_part(x + shiftvector)

  end subroutine random_shift_point


  ! Randomly shift a pointset x with the current shiftvector
  ! NOTE: x should have it's dimensions as it's columns!!!
  !
  subroutine random_shift_pointset(x)
    real(kind=qp), dimension(:,:), intent(inout) :: x

    if (.not. allocated(shiftvector)) then
      allocate(shiftvector(size(x)))
      call random_number(shiftvector)
    end if

    x = frac_part(x + spread(shiftvector, 1, size(x,1)))

  end subroutine random_shift_pointset


  subroutine free_cranley_patterson()
    if (allocated(shiftvector)) then
      deallocate(shiftvector)
    end if
  end subroutine free_cranley_patterson 

    
end module cranley_patterson
