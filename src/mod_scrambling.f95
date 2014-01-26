! Module that implements several things we can use for simplified versions of
! Owen scrambling.
!
! References:
!
!   [1] `On the L_2-Discrepancy of Anchored Boxes', Jiri Matousek, Journal of
!       Complexity, Vol. 14, 1998, pages 527-556.
!
module mod_scrambling

  use numeric_kinds
  use mod_utilities

  private

  public :: random_linear_scrambling
  public :: random_linear_digit_scrambling
  public :: left_i_binomial_scrambling
  public :: striped_matrix_scrambling

  ! TODO:
  !   * create a subroutine generate_random_shift(shift)

contains


  ! Construct a 'Random Linear Scrambling' as described in [1].  The shift
  ! vector is optional.
  !
  subroutine random_linear_scrambling(b, x, shift)
    integer(kind=i4b), intent(in)                          :: b
    integer(kind=i4b), dimension(:,:), intent(out)         :: x
    integer(kind=i4b), dimension(:), intent(out), optional :: shift

    integer(kind=i4b) :: i, j
    integer(kind=i4b) :: temp

    x = 0

    ! Set diagonal elements randomly from [1,..,b-1].
    do i = 1, size(x, dim=1)
      call random_integer(b-1, temp)
      temp = temp + 1
      x(i, i) = temp
    end do

    ! Set non-diagonal lower-triangular elements randomly from [0,..,b-1].
    do i = 1, size(x, dim=1)
      do j = 1, i-1
        call random_integer(b, x(i,j))
      end do
    end do

    ! If a shift is asked, assign it random integers from [0,...,b-1].
    if (present(shift)) then
      shift = 0
      do i = 1, size(shift)
        call random_integer(b, shift(i))
      end do
    end if

  end subroutine random_linear_scrambling


  ! Construct a 'Random Linear Digit Scrambling' as described in [1].  The
  ! shift-vector is optional.
  !
  subroutine random_linear_digit_scrambling(b, x, shift)
    integer(kind=i4b), intent(in)                          :: b
    integer(kind=i4b), dimension(:,:), intent(out)         :: x
    integer(kind=i4b), dimension(:), intent(out), optional :: shift

    integer(kind=i4b) :: i
    integer(kind=i4b) :: temp

    x = 0
    ! Set diagonal elements randomly from [1,..,b-1].
    do i = 1, size(x, dim=1)
      call random_integer(b-1, temp)
      temp = temp + 1
      x(i, i) = temp
    end do

    ! If a shift is asked, assign it random integers from [0,...,b-1].
    if (present(shift)) then
      shift = 0
      do i = 1, size(shift)
        call random_integer(b, shift(i))
      end do
    end if

  end subroutine random_linear_digit_scrambling


  subroutine left_i_binomial_scrambling(b, x, shift)
    integer(kind=i4b), intent(in)                          :: b
    integer(kind=i4b), dimension(:,:), intent(out)         :: x
    integer(kind=i4b), dimension(:), intent(out), optional :: shift

    integer(kind=i4b) :: i, j, diagoffset, n
    integer(kind=i4b) :: temp

    n = size(x, dim=1)

    x = 0

    ! Set diagonal elements to 1 random element from [1,...,b-1]
    call random_integer(b-1, temp)
    temp = temp + 1
    do i = 1, n
      x(i, i) = temp
    end do

    ! Assign all other elements 'diagonally' random integers from [0,...,b-1]
    do diagoffset = 1, n-1
      call random_integer(b, temp)
      do j = 1, n - diagoffset
        x(j+diagoffset,j) = temp
      end do
    end do

    ! If a shift is asked, assign it random integers from [0,...,b-1].
    if (present(shift)) then
      shift = 0
      do i = 1, n
        call random_integer(b, shift(i))
      end do
    end if

  end subroutine left_i_binomial_scrambling


  subroutine striped_matrix_scrambling(b, x, shift)
    integer(kind=i4b), intent(in)                          :: b
    integer(kind=i4b), dimension(:,:), intent(out)         :: x
    integer(kind=i4b), dimension(:), intent(out), optional :: shift

    integer(kind=i4b) :: i, j
    integer(kind=i4b) :: temp

    x = 0

    ! Set the elements column by column randomly from [1,...,b-1]
    do i = 1, size(x, dim=1)
      call random_integer(b-1, temp)
      temp = temp + 1
      do j = i, size(x, dim=1)
        x(j,i) = temp
      end do
    end do

    ! If a shift is asked, assign it random integers from [0,...,b-1].
    if (present(shift)) then
      shift = 0
      do i = 1, size(shift)
        call random_integer(b, shift(i))
      end do
    end if

  end subroutine striped_matrix_scrambling


end module mod_scrambling
