module mod_debug

  use numeric_kinds
  use mod_utilities

  private

  public :: print_matrix
  public :: print_vector

  private :: print_integer_matrix
  private :: print_real_matrix
  private :: print_integer_vector

  interface print_matrix
    module procedure print_integer_matrix, print_real_matrix
  end interface print_matrix

  interface print_vector
    module procedure print_integer_vector
  end interface print_vector

contains


  subroutine print_integer_matrix(x)
    integer(kind=i4b), dimension(:,:), intent(in) :: x

    integer(kind=i4b)               :: m, w
    integer(kind=i4b)               :: row, i
    character(len=50)               :: fmtstring
    integer(kind=i4b), dimension(2) :: loc

    !print *, "--- SIMPLE MATRIX DUMP ----"
    !do row = 1, size(x, 1)
    !  print *, x(row,:)
    !end do
    !print *, "---------------------------"

    ! NEW SAFE CODE (uses one extra space if the maximum is positive)
    !m = maxval(abs(x))
    !w = get_nb_digits(m, 10) + 1

    ! OLD CODE (doesn't use extra space)
    ! Get the number of digits of the largest number in absolute value and
    ! add 1 to the width if this number is negative (for the minus-sign).
    loc = maxloc(abs(x))
    !print *, "max location is ", loc
    !print *, "loc(1) = ", loc(1)
    !print *, "loc(2) = ", loc(2)
    m = x(loc(1), loc(2))
    !print *, "max abs value is ", m
    w = get_nb_digits(m, 10)
    if (m < 0) then
      w = w + 1
    end if

    write(unit=fmtstring, fmt="(A, I0, A, I0, A)") &
      "(A, ", size(x, 2), "I", w+1, ", A)"

    do i=1,(w+1)*size(x, 2)+2
      write(unit=*, fmt="(A1)", advance="no") "-"
    end do
    write(unit=*, fmt="(A)") ""
    do row = 1, size(x, 1)
      write(unit=*, fmt=fmtstring) "[", x(row, :), "]"
    end do
    do i=1,(w+1)*size(x, 2)+2
      write(unit=*, fmt="(A1)", advance="no") "-"
    end do
    write(unit=*, fmt="(A)") ""

  end subroutine print_integer_matrix


  subroutine print_real_matrix(x, fmtstring)
    real(kind=qp), dimension(:,:), intent(in) :: x
    character(len=*), intent(in)              :: fmtstring

    integer(kind=i4b) :: row

    do row = 1, size(x, 1)
      write(unit=*, fmt=fmtstring) x(row,:)
    end do

  end subroutine print_real_matrix


  subroutine print_integer_vector(v, method)
    integer(kind=i4b), dimension(:), intent(in) :: v
    character(len=*), intent(in)                :: method

    integer(kind=i4b), dimension(1) :: m
    integer(kind=i4b)               :: w
    integer(kind=i4b), dimension(1) :: loc
    character(len=50)               :: fmtstring

    ! Get the number of digits of the largest number in absolute value and
    ! add 1 to the width if this number is negative (for the minus-sign).
    loc = maxloc(abs(v))
    m = v(loc)
    w = get_nb_digits(m(1), 10)
    if (m(1) < 0) then
      w = w + 1
    end if

    if (method == "transposed") then
      write(unit=fmtstring, fmt="(A, I0, A, I0, A)") &
        "(", size(v), "I", w+1, ")"
    else
      write(unit=fmtstring, fmt="(A, I0, A)") &
        "(I", w+1, ")"
    end if
    write(unit=*, fmt=fmtstring) v

  end subroutine print_integer_vector


end module mod_debug
