! Module to be used for the bookkeeping of statistical data like averages
! etc...
!
! TODO:
!   * Add variance calculation.
!   * Check the performance of this module, see if it can be made more
!     efficient.

module mod_statistics

  use numeric_kinds

!  implicit none

  private

  public :: create_statistic
  public :: get_statistic
  public :: update_statistic
  public :: print_statistic
  public :: free_statistics

  type, public :: statistic
    character(len=30)                    :: sname
    integer(kind=i4b)                    :: nb_values
    real(kind=qp)                        :: current_sum
    real(kind=qp)                        :: average
    real(kind=qp)                        :: variance 
    real(kind=qp), dimension(:), pointer :: values
  end type statistic

  type(statistic), dimension(:), allocatable, private :: statistics

  contains


    ! Create a new statistic with the specified name and initialize
    ! everything to zero.
    !
    subroutine create_statistic(stats_name)

      character(len=*), intent(in) :: stats_name 

      type(statistic)                            :: newstat
      type(statistic), dimension(:), allocatable :: temp
      integer(kind=i4b)                          :: n

      newstat%sname = stats_name
      newstat%nb_values = 0
      newstat%current_sum = 0.0_qp
      newstat%average = 0.0_qp
      newstat%variance = 0.0_qp
      newstat%values => null()

      if (.not. allocated(statistics)) then
        allocate(statistics(1))
        statistics(1) = newstat
      else
        n = size(statistics)
        allocate(temp(n+1)) 
        temp(1:n) = statistics
        temp(n+1) = newstat
        deallocate(statistics)
        allocate(statistics(n+1))
        statistics = temp
        deallocate(temp)
      end if

    end subroutine create_statistic


    ! Lookup the statistic in our list of statistics and updated it via
    ! the value x.
    !
    subroutine update_statistic(stats_name, x)

      character(len=*), intent(in) :: stats_name 
      real(kind=qp), intent(in)    :: x

      integer(kind=i4b)                    :: i
      integer(kind=i4b)                    :: n
      real(kind=qp), dimension(:), pointer :: temp

      i = 1
      do
        if ((stats_name == statistics(i)%sname) .or. (i>size(statistics))) then
          exit
        end if
        i = i+1
      end do
      if (i>size(statistics)) then
        print *, "ERROR while updating statistic ", stats_name, &
                 ": statistic not found!" 
      else

        ! Add the value to the list of values
        if (associated(statistics(i)%values)) then
          n = size(statistics(i)%values)
          allocate(temp(n+1))
          temp(1:n) = statistics(i)%values
          temp(n+1) = x
          deallocate(statistics(i)%values)
          statistics(i)%values => temp
        else
          allocate(statistics(i)%values(1))
          statistics(i)%values(1) = x
        end if
        ! Update the rest of the statistics
        statistics(i)%nb_values = statistics(i)%nb_values + 1
        statistics(i)%current_sum = statistics(i)%current_sum + x
        statistics(i)%average=statistics(i)%current_sum/statistics(i)%nb_values
      end if

    end subroutine update_statistic


    subroutine print_statistic(stats_name)

      character(len=*), intent(in) :: stats_name

      type(statistic) :: stat

      call get_statistic(stats_name, stat) 
      print *, "Statistic '", stats_name, "':"
      print *, "  Average = ", stat%average
      print *, "  Values = ", stat%values
      
    end subroutine print_statistic


    ! Retrieve the statistic specified by SNAME and put
    ! it in STAT.
    !
    subroutine get_statistic(sname, stat)
      character(len=*), intent(in) :: sname
      type(statistic), intent(out) :: stat

      integer(kind=i4b) :: i

      i = 1
      do
        if ((sname == statistics(i)%sname) .or. (i>size(statistics))) then
          exit
        end if
        i = i+1
      end do
      if (i>size(statistics)) then
        print *, "ERROR: could not find statistic '", sname, "'."
      else
        stat = statistics(i)
      end if

    end subroutine get_statistic


    ! Free up memory used by all statistics.
    !
    subroutine free_statistics()
      integer(kind=i4b) :: i

      if (allocated(statistics)) then
        do i=1,size(statistics)
          deallocate(statistics(i)%values)
        end do
        deallocate(statistics)
      end if

    end subroutine free_statistics


end module mod_statistics
