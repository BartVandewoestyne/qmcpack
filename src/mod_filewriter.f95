module mod_filewriter

  use numeric_kinds
  use mod_file_utils
  use mod_stringutil
  use mod_integration

  private

  public :: write_vector
  public :: write_sequence
  public :: write_integration_result

  character(len=50), private :: fmtstring

contains


  ! Write a single vector in a column into the specified file with the
  ! specified format descriptor.
  !
  subroutine write_vector(v, filename, fmtstring)
    real(kind=qp), dimension(:), intent(in) :: v
    character(len=*), intent(in)            :: filename
    character(len=*), intent(in)            :: fmtstring

    integer(kind=i4b) :: i
    integer           :: my_unit

    call get_unit(my_unit)
    call checked_open(my_unit, filename, "write")

    do i=1,size(v,1)
      write(unit=my_unit, fmt=fmtstring) v(i)
    end do

    call checked_close(my_unit)

  end subroutine write_vector


  ! Write a sequence of points to a file.
  !
  subroutine write_sequence(x, seqname, n, s)
    real(kind=qp), dimension(:,:), intent(in) :: x
    character(len=*), intent(in)              :: seqname
    integer(kind=i4b), intent(in)             :: n, s
      
    character(len=MAX_FILENAME_LENGTH) :: filename
    integer(kind=i4b)                  :: i
    integer                            :: my_unit

    write(unit=filename, fmt="(A, I0, A, I0, A)") &
                                        seqname//"_", n, "p_", s, "d.dat"

    call get_unit(my_unit)
    call checked_open(my_unit, trim(filename), "write")

    write(unit=my_unit, fmt="(A20, A20)") "# Type:             ", seqname 
    write(unit=my_unit, fmt="(A20, I20)") "# Number of points: ", n
    write(unit=my_unit, fmt="(A20, i20)") "# Dimension:        ", s
    !write(unit=my_unit, fmt="(A, I10)") "# Startindex: ", startindex

    write(unit=fmtstring, fmt="(A, I0, A)") "(", s, "F20.15)"
    do i=1,n
      write(unit=my_unit, fmt=fmtstring) x(i,:)
    end do

    call checked_close(my_unit)

  end subroutine write_sequence


  subroutine write_integration_result(filename, i_result, integral, dimen, &
                                      exact)
    character(len=*), intent(in)                       :: filename
    type(integration_result), dimension(:), intent(in) :: i_result
    character(len=*), intent(in), optional             :: integral
    integer(kind=i4b), intent(in), optional            :: dimen
    real(kind=qp), intent(in), optional                :: exact

    integer(kind=i4b) :: my_unit
    integer(kind=i4b) :: i

    call get_unit(my_unit)
    call checked_open(my_unit, filename, "write")

    if (present(integral)) then
      write(unit=my_unit, fmt="(A1, A25, A50)") &
        "#", "Integral:                ", integral
    end if

    if (present(dimen)) then
      write(unit=my_unit, fmt="(A1, A25, I50)") &
        "#", "Dimension:               ", dimen
    end if

    if (present(exact)) then
      write(unit=my_unit, fmt="(A1, A25, ES50.15)") &
        "#", "Exact value of integral: ", exact
    end if

    write(unit=my_unit, fmt="(A1)") "#"
    write(unit=my_unit, fmt="(A1, A14, 3A25)") "#", "n", "integral", &
                                               "absolute_error", &
                                               "relative_error"

    do i=1,size(i_result)
      write(unit=my_unit, fmt="(i15, 3(tr1, es25.15))") i_result(i)%n, &
                                                        i_result(i)%i, &
                                                        i_result(i)%abs_err, &
                                                        i_result(i)%rel_err
    end do
    call checked_close(my_unit)

  end subroutine write_integration_result

end module mod_filewriter
