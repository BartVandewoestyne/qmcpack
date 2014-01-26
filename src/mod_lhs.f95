! Module that implements several versions of Latin Hypercube Sampling.
! 
! References:
!
!   [1] `Randomly Permuted (t,m,s)-Nets and (t,s)-Sequences', Owen, Art B.,
!       in 'Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing',
!       Niederreiter, Harald and P. J.-S. Shiue, Springer-Verlag, New York,
!       1995, pages 299-317.
!
!   [2] `RandQMC User's guide, A package for randomized quasi-Monte Carlo
!       methods in C', Christiane Lemieux, Mikolaj Cieslak, Kristopher Luttmer,
!       University of Calgary, January 13, 2004.

module mod_lhs

  ! Note: There is no next_lhs() function here, since the points depend on n.

  use numeric_kinds
  use mod_utilities

  private

  public :: lhs

  contains

    ! lhs_type must be one of the following strings:
    !
    !   'Patterson' - Patterson (1954)
    !   'MCB'       - McKay, Conover and Beckman (1979)
    !   'Modified'  - Modified Latin Hypercube Sampling (See [2], page 9)
    !
    subroutine lhs(n, s, lhs_type, x)
      integer(kind=i4b), intent(in)              :: n
      integer(kind=i4b), intent(in)              :: s
      character(len=*), intent(in)               :: lhs_type
      real(kind=qp), dimension(:,:), intent(out) :: x

      real(kind=qp), dimension(n, s)  :: myrand_n
      real(kind=qp), dimension(s)     :: myrand_1
      integer(kind=i4b), dimension(n) :: dummy
      integer(kind=i4b)               :: i

      do i=1,s
        call randperm(n, dummy)
        x(:,i) = dummy-1.0_qp
      end do

      select case (lhs_type)

        case ("MCB")

          call random_number(myrand_n)
          x = (x + myrand_n)/n

        case ("Patterson")

          x = (x + 0.5_qp)/n

        case ("Modified")

          call random_number(myrand_1)
          x = frac_part( x/n + spread(myrand_1, 1, n) )

        case default

          stop "Please specify a correct Latin Hypercube Sample type!"

      end select

    end subroutine lhs

end module mod_lhs
