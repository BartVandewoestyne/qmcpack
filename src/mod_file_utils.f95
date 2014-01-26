! Module that implements several usefull file-stuff.

module mod_file_utils

  use numeric_kinds

  private

  public :: get_unit
  public :: checked_open
  public :: checked_close

  integer(kind=i1b), public, parameter :: MAX_FILENAME_LENGTH=100
  integer(kind=i4b), public, parameter :: MAX_RECORD_LENGTH=10000
  character(len=*), public, parameter  :: FILE_SEPARATOR="/"
  integer, public                      :: unit_number

  contains


    ! Return a free unit (= a unit that has no opened file already
    ! connected with it).
    !
    ! Note that the following situation will result in the same unit:
    !
    !   call get_unit(my_first_unit)
    !   call get_unit(my_second_unit)
    !
    ! After calling get_unit, the file should immediately be opened.  Then
    ! a second unit can be requested.
    !
    subroutine get_unit(iunit)
      integer(kind=i4b), intent(out) :: iunit

      logical           :: is_open
      integer(kind=i4b) :: i
      integer           :: ios

      iunit = 0

      do i = 1, 99

        ! For NagWare's f95 we have:
        !  5 = stdin
        !  6 = stdout
        ! For g95 we have:
        !  5 = stdin
        !  6 = stdout
        ! For the F-compiler (Release 20031017) we have:
        !  0 = stderr
        !  5 = stdin
        !  6 = stdout

        if ( i /= 5 .and. i /= 6 ) then

          inquire ( unit = i, opened = is_open, iostat = ios )

          if ( ios == 0 ) then
            if ( .not. is_open ) then
              iunit = i
              return
            end if
          end if

        end if

      end do

      return

    end subroutine get_unit


    ! Open a file on the specified file unit for the specified action and
    ! check if there occured errors.
    ! 
    ! Note:
    !   If opening the file failed, then the optional argument SUCCESS can
    !   be used to let the programmer decide if he will take the
    !   responsibility to do something with the error, or if he will simply
    !   let the program stop.
    !
    subroutine checked_open(file_unit, file_name, file_action, success)
      integer(kind=i4b), intent(in)  :: file_unit
      character(len=*), intent(in)   :: file_name
      character(len=*), intent(in)   :: file_action
      logical, intent(out), optional :: success

      integer :: ios

      if (present(success)) then
        success = .true.
      end if

      open(unit=file_unit, file=file_name, iostat=ios, status="replace", &
           access="sequential", action=file_action, position="rewind")
      if (ios /= 0) then
        write(unit=*, fmt="(A)") "ERROR: could not open file "//file_name//"!"
        if (present(success)) then
          success = .false.
        else
          stop
        end if
      end if
      
    end subroutine checked_open


    ! Disconnect a file from the specified unit and perform
    ! a check if this went OK.
    !
    subroutine checked_close(file_unit, success)
      integer(kind=i4b), intent(in) :: file_unit
      logical, intent(out), optional :: success

      integer :: ios

      close(unit=file_unit, iostat=ios)

      if (present(success)) then
        success = .true.
      end if

      if (ios /= 0) then
        write(unit=*, fmt="(A, I0.0, A)") &
          "ERROR: could not close file on unit ", file_unit, "!"
        if (present(success)) then
          success = .false.
        else
          stop
        end if
      end if

    end subroutine checked_close

end module mod_file_utils
