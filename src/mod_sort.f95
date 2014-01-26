! References:
!
!   [1] `Introduction to Algorithms', Thomas H. Cormen, Charles, E. Leiserson,
!       Ronald L. Rivest, The MIT Press, Eighteenth printing 1997,
!       ISBN 0-262-53091-0, 1997 printing.
!
! TODO:
!   Check why NagWare's f95 gives this warning:
!     mod_sort.f95: In function `mod_sort_MP_quicksort':
!     mod_sort.f95:36: warning: pointer type mismatch in conditional expression

module mod_sort

  use numeric_kinds

  private

  public :: quicksort
  public :: is_increasing

  private :: partition

  contains

    ! Sorts real numbers into ascending numerical order (based on the
    ! algorithm in [1]).  Worst case running time is theta(n^2), but
    ! average case running time is theta(n*log(n)) and it generally
    ! outperforms heapsort in practice (See [1], page 138).
    !
    ! If before sorting, indices = 1, 2, 3, ..., N, then after sorting
    ! we have that x_sorted(i) = x(indices(i))
    !
    recursive subroutine quicksort(x, indices)
      real(kind=qp), dimension(:), intent(inout)               :: x
      integer(kind=i4b), dimension(:), intent(inout), optional :: indices

      integer(kind=i4b) :: q

      if (size(x) > 1) then
        if (present(indices)) then
          call partition(x, q, indices)
          call quicksort(x(:q), indices(:q))
          call quicksort(x(q+1:), indices(q+1:))
        else
          call partition(x, q)
          call quicksort(x(:q))
          call quicksort(x(q+1:))
        end if
      end if

    end subroutine quicksort


    ! See [1], page 154-155 for how this private subroutine works.
    !
    subroutine partition(A, q, indices)
      real(kind=qp), dimension(:), intent(inout)               :: A
      integer(kind=i4b), intent(out)                           :: q
      integer(kind=i4b), dimension(:), intent(inout), optional :: indices

      integer(kind=i4b) :: i, j
      real(kind=qp)     :: temp
      integer(kind=i4b) :: temp_int
      real(kind=qp)     :: x


      ! Initially, set the markers just off the left and right ends of the
      ! array, so the two regions are empty.
      i = lbound(A, dim=1) - 1
      j = ubound(A, dim=1) + 1

      ! Select the first element as the `pivot' element around which to
      ! partition the array.
      x = A(i+1)

      do

        ! Make the rightmost region larger
        do
          j = j-1
          if (A(j) <= x) then
            exit
          end if
        end do

        ! Make the leftmost region larger
        do
          i = i+1
          if (A(i) >= x) then
            exit
          end if
        end do

        if (i<j) then

          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp

          ! exchange also the indices in case a user
          ! would need these...
          if (present(indices)) then
            temp_int = indices(i)
            indices(i) = indices(j)
            indices(j) = temp_int
          end if

        else

          ! At this point, the entire array has been partitioned into two
          ! subarrays such that no element of the leftmost subarray is larger
          ! than any element of the rightmost subarray, so we return j as the
          ! last bound of the leftmost subarray.

          q = j

          return
        end if
      end do

    end subroutine partition


    ! Return true if the array has increasing real numbers.
    !
    pure function is_increasing(x) result (y)

      real(kind=qp), dimension(:), intent(in) :: x
      logical                                 :: y

      integer(kind=i4b) :: i

      y = .true.

      i = 1
      do
        if (i == size(x)) then
          exit
        end if
        if (x(i+1)<x(i)) then
          y = .false.
          return
        end if
        i = i+1
      end do

    end function is_increasing

end module mod_sort
