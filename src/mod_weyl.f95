! Module that implements different kinds of Weyl-sequences, also called
! 'gp-sets'.
!
! References:
!
!   [1] "Uniform Randon Numbers: Theory And Practice", Shu Tezuka, p. 43.
!
!   [2] 'Number-theoretic Methods in Statistics', Wang, Y. and Fang, K.-T.,
!       1994, page 27.
!
! TODO:
!   * Add Nested Weyl sequence and Shuffled Nested Weyl sequence.

module mod_weyl

  use numeric_kinds
  use mod_primes
  use mod_utilities
  use mod_constants

  private

  public :: init_weyl
  public :: next_weyl
  public :: init_square_root_seq
  public :: next_square_root
  public :: square_root_seq
  public :: init_gp_set_b
  public :: next_gp_set_b
  public :: gp_set_b
  public :: init_gp_set_c
  public :: next_gp_set_c
  public :: gp_set_c
  public :: free_weyl


  real(kind=qp), dimension(:), allocatable, private     :: wey_irrationals
  integer(kind=i4b), dimension(:), allocatable, private :: wey_n
  integer(kind=i4b), dimension(:), allocatable, private :: wey_step


  contains

    subroutine init_weyl(s, irrationals, n_start, step)
      integer(kind=i4b), intent(in)                         :: s
      real(kind=qp), dimension(:), intent(in)               :: irrationals
      integer(kind=i4b), dimension(:), intent(in), optional :: n_start
      integer(kind=i4b), dimension(:), intent(in), optional :: step

      if (allocated(wey_irrationals)) then
        deallocate(wey_irrationals)
      end if
      allocate(wey_irrationals(s))

      if (allocated(wey_n)) then
        deallocate(wey_n)
      end if
      allocate(wey_n(s))

      if (allocated(wey_step)) then
        deallocate(wey_step)
      end if
      allocate(wey_step(s))

      wey_irrationals = irrationals

      if (present(n_start)) then
        wey_n = n_start
      else
        wey_n = 1
      end if

      if (present(step)) then
        wey_step = step
      else
        wey_step = 1
      end if

    end subroutine init_weyl


    ! Get next element from the Weyl sequence.
    !
    subroutine next_weyl(x)

      real(kind=qp), dimension(:), intent(out) :: x

      ! Check for overflow
      !if (maxval(wey_n) > huge(qp)/maxval(wey_irrationals)) then
      !  stop "N is too big to be usable for this architecture..."
      !end if

      x = frac_part(wey_n*wey_irrationals)

      wey_n = wey_n + wey_step

    end subroutine next_weyl


    ! Initialize square-root sequence
    !
    subroutine init_square_root_seq(s, myprimes, n_start, step)
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), dimension(:), intent(in), optional :: myprimes
      integer(kind=i4b), dimension(:), intent(in), optional :: n_start
      integer(kind=i4b), dimension(:), intent(in), optional :: step

      real(kind=qp), dimension(s) :: local_sqrt_primes

      if (present(myprimes)) then
        local_sqrt_primes = sqrt(real(myprimes, kind=qp))
      else
        local_sqrt_primes = sqrt(real(primes(s), kind=qp))
      end if

      call init_weyl(s, local_sqrt_primes, n_start, step)

    end subroutine init_square_root_seq


    ! Get next square-root sequence element
    !
    subroutine next_square_root(x)

      real(kind=qp), dimension(:), intent(out) :: x

      call next_weyl(x)

    end subroutine next_square_root


    ! The Square Root Sequence
    ! (A Weyl-sequence with sqrt(prime(s)) as the irrational numbers.)
    !
    subroutine square_root_seq(n, s, n_start, p, x)
      integer(kind=i4b), intent(in)                         :: n
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), dimension(:), intent(in), optional :: n_start
      integer(kind=i4b), dimension(:), intent(in), optional :: p
      real(kind=qp), dimension(:,:), intent(out)            :: x

      integer(kind=i4b)               :: i, j
      integer(kind=i4b), dimension(s) :: loc_n_start
      integer(kind=i4b), dimension(s) :: loc_p

      if (.not. present(n_start)) then
        loc_n_start = 1
      else
        loc_n_start = n_start
      end if
    
      if (.not. present(p)) then
        loc_p = primes(s)
      else
        loc_p = p
      end if

      do i=1,s
       x(:,i) = frac_part((loc_n_start(i)+(/(j,j=0,n-1)/))*sqrt(real(loc_p(i), kind=qp)))
      end do

      call free_weyl()

    end subroutine square_root_seq


    ! Initialize gp set (b) from [2]
    !
    subroutine init_gp_set_b(s, p, n_start, step)
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), intent(in)                         :: p
      integer(kind=i4b), dimension(:), intent(in), optional :: n_start
      integer(kind=i4b), dimension(:), intent(in), optional :: step

      integer(kind=i4b)           :: i
      real(kind=qp), dimension(s) :: gmma

      gmma = (/ (p**(i/(s+1.0_qp)), i=1,s) /)
      call init_weyl(s, gmma, n_start, step)

    end subroutine init_gp_set_b


    subroutine gp_set_b(n, s, p, n_start, step, x)
      integer(kind=i4b), intent(in)                         :: n
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), intent(in)                         :: p
      integer(kind=i4b), dimension(:), intent(in), optional :: n_start
      integer(kind=i4b), dimension(:), intent(in), optional :: step
      real(kind=qp), dimension(:,:), intent(out)            :: x

      integer(kind=i4b)               :: i, j
      real(kind=qp), dimension(s)     :: gmma
      integer(kind=i4b), dimension(s) :: loc_n_start
      integer(kind=i4b), dimension(s) :: loc_step

      if (.not. present(n_start)) then
        loc_n_start = 1
      else
        loc_n_start = n_start
      end if

      if (.not. present(step)) then
        loc_step = 1
      else
        loc_step = step
      end if

      gmma = (/ (p**(i/(s+1.0_qp)), i=1,s) /)

      do i=1,s
        x(:,i) = frac_part((loc_n_start(i)+(/(j*loc_step(i),j=0,n-1)/))*gmma(i))
      end do
      
    end subroutine gp_set_b


    ! Initialize gp set (c) from [2]
    !
    subroutine init_gp_set_c(s, p, n_start, step)
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), intent(in)                         :: p
      integer(kind=i4b), dimension(:), intent(in), optional :: n_start
      integer(kind=i4b), dimension(:), intent(in), optional :: step

      integer(kind=i4b)           :: i
      real(kind=qp), dimension(s) :: gmma

      if (p < 2*s+3) then
        print *, "ERROR: in init_gp_set_c(), p should be a prime larger"
        print *, "       or equal to 2*s+3 !!!"
        print *, "       Could not initialize gp_set_c!!!"
      else
        gmma = frac_part( (/ (2*cos(2*pi*i/p), i=1,s) /) )
        call init_weyl(s, gmma, n_start, step)
      end if

    end subroutine init_gp_set_c


    subroutine gp_set_c(n, s, p, n_start, step, x)
      integer(kind=i4b), intent(in)                         :: n
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), intent(in)                         :: p
      integer(kind=i4b), dimension(:), intent(in), optional :: n_start
      integer(kind=i4b), dimension(:), intent(in), optional :: step
      real(kind=qp), dimension(:,:), intent(out)            :: x

      integer(kind=i4b)               :: i, j
      real(kind=qp), dimension(s)     :: gmma
      integer(kind=i4b), dimension(s) :: loc_n_start
      integer(kind=i4b), dimension(s) :: loc_step

      if (p < 2*s+3) then
        print *, "ERROR: in gp_set_c(), p should be a prime larger"
        print *, "       or equal to 2*s+3 !!!"
        print *, "       Could not generate points for gp_set_c!!!"
      else
        if (.not. present(n_start)) then
          loc_n_start = 1
        else
          loc_n_start = n_start
        end if

        if (.not. present(step)) then
          loc_step = 1
        else
          loc_step = step
        end if

        gmma = frac_part( (/ (2*cos(2*pi*i/p), i=1,s) /) )

        do i=1,s
          x(:,i) = frac_part((loc_n_start(i)+loc_step(i)*(/(j,j=0,n-1)/))*gmma(i))
        end do
      end if

      call free_weyl()

    end subroutine gp_set_c



    ! Get next element from the gp set (b) from [2]
    !
    subroutine next_gp_set_b(x)

      real(kind=qp), dimension(:), intent(out) :: x

      call next_weyl(x)

    end subroutine next_gp_set_b


    ! Get next element from the gp set (c) from [2]
    !
    subroutine next_gp_set_c(x)

      real(kind=qp), dimension(:), intent(out) :: x

      call next_weyl(x)

    end subroutine next_gp_set_c


    subroutine free_weyl()

      if (allocated(wey_irrationals)) then
        deallocate(wey_irrationals)
      end if
      if (allocated(wey_n)) then
        deallocate(wey_n)
      end if
      if (allocated(wey_step)) then
        deallocate(wey_step)
      end if

    end subroutine free_weyl


end module mod_weyl
