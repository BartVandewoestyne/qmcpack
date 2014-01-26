! Module that implements the Haber sequence.
!
! Notes:
! 
!   * Be carefull for overflow, because the definition of the Haber sequence
!     contains 0.5*k*(k+1) which grows large quite fast if k gets large.
!
!   * It seems like we also loose precision quite fast... this is probably
!     due to the fact that we keep only the fractional part of some very
!     large numbers.  See also test_haber.f95 in the tests directory.
!
! References:
!
!   [1] 'A quasirandom approach to integration in Bayesian statistics',
!       Shaw, J., Annals of Statistics, vol. 16, pages 895-914, 1988.
!
!   [2] 'Number-theoretic Methods in Statistics', Wang, Y. and Fang, K.-T.,
!       1994, page 31.

module mod_haber

  use numeric_kinds
  use mod_primes
  use mod_utilities

  private

  public :: init_haber
  public :: haber
  public :: next_haber
  public :: free_haber
  public :: check_for_overflow

  integer(kind=i4b), private                            :: hab_s
  integer(kind=i4b), dimension(:), allocatable, private :: hab_startindex
  integer(kind=i4b), dimension(:), allocatable, private :: hab_primes
  integer(kind=i4b), private                            :: hab_point_counter

  contains


    subroutine init_haber(s, startindex, myprimes)
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), dimension(:), intent(in), optional :: startindex
      integer(kind=i4b), dimension(:), intent(in), optional :: myprimes

      hab_point_counter = 0

      ! Initialize the dimension of the points to be generated
      hab_s = s

      ! Initialize the startindices
      if (allocated(hab_startindex)) then
        if (size(hab_startindex) /= s) then
          deallocate(hab_startindex)
        end if
      end if
      if (.not. allocated(hab_startindex)) then
        allocate(hab_startindex(s))
      end if
      if (present(startindex)) then
        hab_startindex = startindex
      else
        ! Use (1,...,1) as default startindex
        hab_startindex = 1
      end if

      ! Initialize the primes to be used
      if (allocated(hab_primes)) then
        if (size(hab_primes) /= s) then
          deallocate(hab_primes)
        end if
      end if
      if (.not. allocated(hab_primes)) then
        allocate(hab_primes(s))
      end if
      if (present(myprimes)) then
        hab_primes = myprimes
      else
        ! Use the first s primes by default
        hab_primes = primes(s)
      end if

    end subroutine init_haber


    subroutine haber(n, s, startindex, myprimes, x)
      integer(kind=i4b), intent(in)                         :: n
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), dimension(:), intent(in), optional :: startindex
      integer(kind=i4b), dimension(:), intent(in), optional :: myprimes
      real(kind=qp), dimension(:, :), intent(out)           :: x

      integer(kind=i4b) :: i, j


      call init_haber(s, startindex, myprimes)

      x = spread(sqrt(real(hab_primes, kind=qp)), 1, n)

      do i=1,s
        do j=0,n-1
          x(j+1,i) = (0.5_qp*(startindex(i)+j))*(startindex(i)+j+1)*x(j+1,i)
        end do
      end do

      x = frac_part(x)

      call free_haber()

    end subroutine haber


    subroutine next_haber(x)
      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b), dimension(hab_s) :: i

      call check_for_overflow()

      i = hab_startindex + hab_point_counter

      x = frac_part((0.5_qp*i)*(i+1)*sqrt(real(hab_primes, kind=qp)))

      hab_point_counter = hab_point_counter + 1

    end subroutine next_haber


    subroutine free_haber()

      if (allocated(hab_startindex)) then
        deallocate (hab_startindex)
      end if
      if (allocated(hab_primes)) then
        deallocate (hab_primes)
      end if

    end subroutine free_haber


    ! Check if we are going to generate more points than allowed by
    ! our architecture
    subroutine check_for_overflow()

      if (huge(i4b)-(maxval(hab_startindex)+hab_point_counter+1) == 0) then
       write(unit=*, fmt="(A)") "ERROR: maximum number of Haber points reached!"
      end if

    end subroutine check_for_overflow


end module mod_haber
