! Module that implements the Hammersley point set.
!
! NOTE: there appear to be multiple definitions of the Hammersley point set
! in the literature.  This implementation allows a user to select whatever
! definition he wants to use.
!
! References:
!
!   [1] Niederreiter, Harald, `Random Number Generation and Quasi-Monte
!       Carlo Methods'
!   [2] `On the efficiency of certain quasi-random sequences of points in
!       evaluating multi-dimensional integrals', Halton, J. H., Numerische
!       Mathematik 2, pages 84-90, 1960.
!   [3] `Number-theoretic Methods in Statistics', Wang, Y. and Fang, K.-T.,
!       1994, page 28.

module mod_hammersley

  use numeric_kinds
  use mod_primes
  use mod_radical_inverse

  private

  public :: hammersley
  public :: free_hammersley

  private :: init_hammersley

  integer(kind=i4b), private                            :: ham_s
  integer(kind=i4b), dimension(:), allocatable, private :: ham_primes
  integer(kind=i4b), dimension(:), allocatable, private :: ham_startindex
  integer(kind=i4b), dimension(:), allocatable, private :: ham_leap

contains

  subroutine init_hammersley(s, myprimes, startindex, leap)
    integer(kind=i4b), intent(in) :: s

    ! These arrays should have dimensions s-1
    integer(kind=i4b), dimension(:), intent(in), optional :: myprimes
    integer(kind=i4b), dimension(:), intent(in), optional :: startindex
    integer(kind=i4b), dimension(:), intent(in), optional :: leap

    ham_s = s

    ! Initialize the startindices
    if (allocated(ham_startindex)) then
        if(size(ham_startindex) /= s-1) then
          deallocate(ham_startindex)
        end if
    end if
    if (.not. allocated(ham_startindex)) then
      allocate(ham_startindex(s-1))
    end if
    if (present(startindex)) then
      ham_startindex = startindex
    else
      ! Use (1,...,1) as default startindex
      ham_startindex = 1
    end if
    
    ! Initialize the primes
    if (allocated(ham_primes)) then
        if(size(ham_primes) /= s-1) then
          deallocate(ham_primes)
        end if
    end if
    if (.not. allocated(ham_primes)) then
      allocate(ham_primes(s-1))
    end if
    if (present(myprimes)) then
      ham_primes = myprimes
    else
      ham_primes = primes(s-1)
    end if

    ! Initialize the leaps
    if (allocated(ham_leap)) then
        if(size(ham_leap) /= s-1) then
          deallocate(ham_leap)
        end if
    end if
    if (.not. allocated(ham_leap)) then
      allocate(ham_leap(s-1))
    end if
    if (present(leap)) then
      ham_leap = leap
    else
      ham_leap = 1
    end if
    
  end subroutine init_hammersley


  ! The Hammersley point set.
  !
  ! Notes:

  !   * There appear to be different definitions of this point set
  !     in the literature.  See [1], page 31, [2], page 85 and [3] page 29.
  !     Currently, by setting specific values for 'def', a user can
  !     choose between the definition in [1], [2] or [3].  Note that the
  !     startindices should be chosen appropriately to match the definition.
  !
  !   * MYPRIMES, STARTINDEX and LEAP should have size s-1.
  !
  subroutine hammersley(n, s, myprimes, startindex, leap, def, x)
    integer(kind=i4b), intent(in)                         :: n
    integer(kind=i4b), intent(in)                         :: s
    integer(kind=i4b), dimension(:), intent(in), optional :: myprimes
    integer(kind=i4b), dimension(:), intent(in), optional :: startindex
    integer(kind=i4b), dimension(:), intent(in), optional :: leap 
    character(len=*), intent(in)                          :: def
    real(kind=qp), dimension(:, :), intent(out)           :: x

    integer(kind=i4b) :: i, j, k

    call init_hammersley(s, myprimes, startindex, leap)

    ! The first dimension depends on the used definition of the Hammersley set.
    select case (def)
      case ("Halton")
        x(:, 1) = (/ (real(k, kind=qp)/n, k=1,n) /)
      case ("Niederreiter", "Shaw")
        x(:, 1) = (/ (real(k, kind=qp)/n, k=0,n-1) /)
      case ("Fang_Wang")
        x(:, 1) = (/ ( (2*k-1.0_qp)/(2*n), k=1,n ) /)
    end select

    ! The last dimensions are Halton numbers
    do i=2,ham_s
      do j=0,n-1
        x(j+1,i) = radical_inverse(ham_startindex(i-1) + j*ham_leap(i-1), &
                                   ham_primes(i-1))
      end do
    end do  

    call free_hammersley()

  end subroutine hammersley


  subroutine free_hammersley()
    
    if (allocated(ham_startindex)) then
      deallocate(ham_startindex)
    end if
    if (allocated(ham_primes)) then
      deallocate(ham_primes)
    end if
    if (allocated(ham_leap)) then
      deallocate(ham_leap)
    end if

  end subroutine free_hammersley
  
end module mod_hammersley
