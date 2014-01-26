! This module contains several routines that are usefull for the tests under
! the tests directory.

module mod_testing

  use numeric_kinds
  use mod_file_utils

  private

  public :: show_test_header
  public :: start_test
  public :: stop_test 
  public :: show_test_summary
  public :: compare_pointsets
  public :: write_function_data
  public :: get_format_specifier
  public :: get_nb_points
  public :: get_nb_columns

  contains

    ! Display a simple header indicating what test is going to come.
    !
    subroutine show_test_header(testname)

      character(len=*), intent(in) :: testname

      integer(kind=i4b) :: i

      write(unit=*, fmt="(A)")

      ! print left bar of #'s
      do i=1,(80-10-len_trim(testname))/2
        write(unit=*, fmt="(A1)", advance="no") "#"
      end do

      ! print middle text
      write(unit=*, fmt="(A)", advance="no") " TESTING " // testname // " "

      ! print right bar of #'s
      do i=1,ceiling((80-10-len_trim(testname))/2.0_qp)-1
        write(unit=*, fmt="(A1)", advance="no") "#"
      end do
      write(unit=*, fmt="(A1)") "#"

    end subroutine show_test_header


    subroutine start_test(message)
      character(len=*), intent(in) :: message

      write(unit=*, fmt="(A)") "=> "//message

    end subroutine start_test


    subroutine stop_test()
      
      write(unit=*, fmt="(A)") "   done."

    end subroutine stop_test


    ! Report on how many errors occured during the test.
    !
    subroutine show_test_summary(nberrors)
      integer(kind=i4b), intent(in) :: nberrors

      write(unit=*, fmt=*)

      if (nberrors > 0) then
        write(unit=*, fmt="(A, I0.0, A)") &
          "=====> ", nberrors, " errors found in this test!  Shame on you! <====="
      else
        write(unit=*, fmt="(A)") &
          "=====> No errors found during test.  Well done my friend! <====="
      end if

    end subroutine show_test_summary


    ! Return true if the pointsets are equal up to the specified
    ! absolute tolerance, else return false.
    !
    subroutine compare_pointsets(set1_file, set2_file, tol, equal)

      character(len=*), intent(in) :: set1_file, set2_file
      real(kind=qp), intent(in)    :: tol
      logical, intent(out)         :: equal

      integer(kind=i4b)                        :: unit_set1, unit_set2
      integer(kind=i4b)                        :: nb_points_set1, nb_points_set2
      integer(kind=i4b)                        :: linenb
      integer(kind=i4b)                        :: nb_cols_set1, nb_cols_set2
      real(kind=qp), dimension(:), allocatable :: nb1, nb2
      character(len=20)                        :: format_file1, format_file2
      integer(kind=i4b)                        :: ios

      equal = .false.

      ! Compare number of columns ( = dimension)
      call get_nb_columns(set1_file, nb_cols_set1)
      call get_nb_columns(set2_file, nb_cols_set2)
      if (nb_cols_set1 /= nb_cols_set2) then
        write(unit=*, fmt="(A)") &
          "ERROR: pointsets do not have the same dimension!"
        return
      end if

      ! Compare number of points ( = lines in file)
      call get_nb_points(set1_file, nb_points_set1)
      call get_nb_points(set2_file, nb_points_set2)
      if (nb_points_set1 /= nb_points_set2) then
        write(unit=*, fmt="(A)") &
          "ERROR: pointsets do not have the same number of points!"
        return
      end if

      ! Get the format specifiers for both files
      call get_format_specifier(set1_file, format_file1)
      call get_format_specifier(set2_file, format_file2)

      ! Compare the points itself (= lines in file).
      call get_unit(unit_set1)
      open(unit=unit_set1, file=set1_file, iostat=ios, status="old", &
           access="sequential", action="read", position="rewind")
      if (ios /= 0) then
        write(unit=*, fmt="(A)") "ERROR: could not open file "//set1_file// &
                                 " to compare pointsets!"
      end if
      call get_unit(unit_set2)
      open(unit=unit_set2, file=set2_file, iostat=ios, status="old", &
           access="sequential", action="read", position="rewind")
      if (ios /= 0) then
        write(unit=*, fmt="(A)") "ERROR: could not open file "//set2_file// &
                                 " to compare pointsets!"
      end if

      allocate(nb1(nb_cols_set1), nb2(nb_cols_set2))
      equal = .true.
      linenb = 1
      do
        read(unit=unit_set1, fmt=format_file1, iostat=ios) nb1
        if (ios /= 0) then
          exit
        end if
        read(unit=unit_set2, fmt=format_file2, iostat=ios) nb2
        if (ios /= 0) then
          exit
        end if
        if (any(abs(nb1-nb2) > tol)) then
          equal = .false.
          write(unit=*, fmt="(A)") ""
          write(unit=*, fmt="(A, I0, A)") "  DIFFERENCE: points on line ", linenb, " were different: "
          write(unit=*, fmt="(A)") "  Point from pointset 1:"
          write(unit=*, fmt="(A)", advance="no") "  "
          write(unit=*, fmt=format_file1) nb1
          write(unit=*, fmt="(A)") "  Point from pointset 2:"
          write(unit=*, fmt="(A)", advance="no") "  "
          write(unit=*, fmt=format_file2) nb2
          !exit
        end if
        linenb = linenb + 1
      end do

      close(unit=unit_set1)
      close(unit=unit_set2)

    end subroutine compare_pointsets


    ! Write out a table of x-values and their function values
    ! to the specified file.  The range is from mystart to myend
    ! in steps of step.
    !
    subroutine write_function_data(f, mystart, myend, step, filename)
      interface
        function f(x) result (y)
          use numeric_kinds
          real(kind=qp), intent(in) :: x
          real(kind=qp)             :: y
        end function f
      end interface
      real(kind=qp), intent(in)    :: mystart
      real(kind=qp), intent(in)    :: myend
      real(kind=qp), intent(in)    :: step
      character(len=*), intent(in) :: filename
      
      real(kind=qp) :: x
      integer       :: my_unit

      call get_unit(my_unit)
      call checked_open(my_unit, filename, "write")

      x = mystart
      do
        if (x > myend) then
          exit
        end if
        write(unit=my_unit, fmt="(2ES22.14)") x, f(x)
        x = x + step
      end do

      call checked_close(my_unit)


    end subroutine write_function_data


    ! Return the format specifier that was used for writing the points
    ! to this file.  The format specifier is obtained as follows:
    !
    !   1) Initially, the format specifier for the first coordinate of the
    !      first point is found.  It is assumed that the points use the F-edit
    !      descriptor and are of the form Fw.d.
    !
    !   2) Then the number of columns is determined.
    !
    !   3) The information from 1) and 2) is combined to the final format
    !      specifier.
    !
    subroutine get_format_specifier(my_file, my_fmt)

      character(len=*), intent(in)  :: my_file
      character(len=*), intent(out) :: my_fmt

      integer(kind=i4b) :: my_unit
      integer(kind=i4b) :: nb_cols
      integer(kind=i4b) :: startpos, dot_pos, endpos
      character(len=1)  :: my_char
      integer(kind=i4b) :: fmt_fieldwidth
      integer(kind=i4b) :: fmt_nbdigits
      integer           :: ios


      call get_nb_columns(my_file, nb_cols)

      call get_unit(my_unit)
      open(unit=my_unit, file=my_file, iostat=ios, status="old", &
           access="sequential", action="read", position="rewind")
      if (ios /= 0) then
        write(unit=*, fmt="(A)") "ERROR: could not open file "//my_file// &
                                 " to get format specifier!"
      end if

      startpos = 1
      do
        read(unit=my_unit, fmt="(A1)", advance="no") my_char
        if (my_char /= " ") then
          exit
        end if
        startpos = startpos+1
      end do
      endpos = startpos
      do
        read(unit=my_unit, fmt="(A1)", advance="no") my_char
        if (my_char == " ") then
          exit
        end if
        endpos = endpos + 1
        if (my_char == ".") then
          dot_pos = endpos
        end if
      end do

      close(unit=my_unit)

      fmt_fieldwidth = endpos
      fmt_nbdigits = endpos - dot_pos

      write(unit=my_fmt, fmt="(A, I0.0, A, I0.0, A, I0.0, A)") &
        "(", nb_cols, "F", fmt_fieldwidth, ".", fmt_nbdigits, ")"

    end subroutine get_format_specifier


    ! Return the number of points that are written into a file.
    !
    ! Conventions for the pointset-fileformat are:
    !
    !   * a row represents an s-dimensional point.
    !   * the columns represent the dimensions.
    !   * there should be no empty lines or lines containing things other than
    !     points.
    !   * there can be one or more spaces in front or after a coordinate of
    !     a point.
    !
    subroutine get_nb_points(pointset_file, nb_points)

      character(len=*), intent(in)   :: pointset_file
      integer(kind=i4b), intent(out) :: nb_points
      
      integer(kind=i4b) :: my_unit
      integer           :: ios

      call get_unit(my_unit)
      open(unit=my_unit, file=pointset_file, iostat=ios, status="old", &
           access="sequential", action="read", position="rewind")
      if (ios /= 0) then
        write(unit=*, fmt="(A)") "ERROR: could not open file "//pointset_file &
                               //" in order to determine the number of points."
      end if

      nb_points = 0
      do
        read(unit=my_unit, fmt=*, iostat=ios)
        if (ios /= 0) then
          exit
        end if
        nb_points = nb_points + 1
      end do
      close(unit=my_unit)

    end subroutine get_nb_points


    ! Return the number of columns in the first row of the file containing the
    ! pointset.  Accept this as the number of columns for the whole pointset.
    !
    ! Note:
    !   * We do not check if the remaining rows also have the same number of
    !     columns!!!
    !
    ! See also:
    !
    !   http://groups.google.be/group/comp.lang.fortran/msg/7d609924ff06d938
    !
    subroutine get_nb_columns(pointset_file, nb_columns)

      character(len=*), intent(in)   :: pointset_file
      integer(kind=i4b), intent(out) :: nb_columns

      integer(kind=i4b), parameter :: max_chars = MAX_RECORD_LENGTH
      integer(kind=i4b)            :: my_unit
      integer(kind=i4b)            :: i
      character(len=max_chars)     :: line
      integer(kind=i4b)            :: ios


      call get_unit(my_unit)
      open(unit=my_unit, file=pointset_file, iostat=ios, status="old", &
           access="sequential", action="read", position="rewind")
      if (ios /= 0) then
        write(unit=*, fmt="(A)") ""
        write(unit=*, fmt="(A)") &
          "ERROR: could not open file "//pointset_file &
          //" to determine number of columns!"
      end if

      read(unit=my_unit, fmt="(A)", iostat=ios) line
      if (ios /= 0) then
        write(unit=*, fmt="(A)") ""
        write(unit=*, fmt="(A)") &
          "ERROR: could not read file "//pointset_file &
          //" to determine number of columns!"
      end if
      nb_columns = 0
      do i = 1,max_chars-1
        if (line(i:i) /= " " .and. line(i+1:i+1) == " ") then
          nb_columns = nb_columns+1
        end if
      end do

      close(unit=my_unit, iostat=ios)
      if (ios /= 0) then
        write(unit=*, fmt="(A)") "ERROR while closing file "//pointset_file
      end if

    end subroutine get_nb_columns

end module mod_testing
