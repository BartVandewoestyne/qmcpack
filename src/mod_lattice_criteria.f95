! Module that implements different methods to calculate several criteria for
! measuring the `goodness' of a lattice rule.
!
! References:
!
!   [1] `Lattice Methods for Multiple Integration', Sloan, I. H. and Joe, S.,
!       Oxford Science Publications, 1994, ISBN 0 19 853472 8.
!
!   [2] `On computing the lattice rule criterion R', Joe, S. and Sloan, I. H.,
!       Mathematics of Computation, vol. 59, pages 557-568, 1992.

module mod_lattice_criteria

  use numeric_kinds
  use f_worstcase
  use mod_function
  use mod_utilities
  use mod_constants

  private

  public :: P_2
  public :: P_4
  public :: P_6
  public :: P_infinity
  public :: Ptilde_2
  public :: R

  private :: F

contains


  ! Calculate P_2(z, N).  See [1], page 69-73.
  !
  function P_2(z, N) result(p)

    integer(kind=i4b), dimension(:), intent(in) :: z
    integer(kind=i4b), intent(in)               :: N
    real(kind=qp)                               :: p

    type(functionparams) :: params
    real(kind=qp)        :: my_sum
    integer(kind=i4b)    :: i

    my_sum = 0.0_qp
    do i=0,N-1
      my_sum = my_sum + f_worstcase2_1p(frac_part( (real(i, kind=qp)/N)*z ), params)
    end do
    p = my_sum/N - 1

  end function P_2


  ! Calculate P_4(z, N).  See [1], page 69-73.
  !
  function P_4(z, N) result(p)

    integer(kind=i4b), dimension(:), intent(in) :: z
    integer(kind=i4b), intent(in)               :: N
    real(kind=qp)                               :: p

    type(functionparams) :: params
    real(kind=qp)        :: my_sum
    integer(kind=i4b)    :: i

    my_sum = 0.0_qp
    do i=0,N-1
      my_sum = my_sum + f_worstcase4_1p(frac_part( (real(i, kind=qp)/N)*z ), params)
    end do
    p = my_sum/N - 1

  end function P_4


  ! Calculate P_6(z, N).  See [1], page 69-73.
  !
  function P_6(z, N) result(p)

    integer(kind=i4b), dimension(:), intent(in) :: z
    integer(kind=i4b), intent(in)               :: N
    real(kind=qp)                               :: p

    type(functionparams) :: params
    real(kind=qp)        :: my_sum
    integer(kind=i4b)    :: i

    my_sum = 0.0_qp
    do i=0,N-1
      my_sum = my_sum + f_worstcase6_1p(frac_part( (real(i, kind=qp)/N)*z ), params)
    end do
    p = my_sum/N - 1

  end function P_6


  ! Calculate P_infinity(z, N).  See [1], page 73.
  !
  function P_infinity(z, N) result(p)

    integer(kind=i4b), dimension(:), intent(in) :: z
    integer(kind=i4b), intent(in)               :: N
    real(kind=qp)                               :: p

    type(functionparams) :: params
    real(kind=qp)        :: my_sum
    integer(kind=i4b)    :: i

    my_sum = 0.0_qp
    do i=0,N-1
      my_sum = my_sum + f_worstcase_infinity_1p(frac_part( (real(i, kind=qp)/N)*z ), params)
    end do
    p = my_sum/N - 1

  end function P_infinity


  ! Calculate Ptilde_2(z, N).  See [1], page 73, formula (4.17).
  !
  function Ptilde_2(z, N) result(p)

    integer(kind=i4b), dimension(:), intent(in) :: z
    integer(kind=i4b), intent(in)               :: N
    real(kind=qp)                               :: p

    type(functionparams) :: params
    real(kind=qp)        :: my_sum
    integer(kind=i4b)    :: i

    my_sum = 0.0_qp
    do i=0,N-1
      my_sum = my_sum + f_saltykov_1p(frac_part( (real(i, kind=qp)/N)*z ), params)
    end do
    p = my_sum/N - 1

  end function Ptilde_2


  ! Calculate R(z, N).  This subroutine currently requires O(N^2) operations
  ! because we are not using the more efficient method presented in [2].
  !
  ! See [1], formula (4.23) page 76 and [2].
  !
  function R(z, N) result (res)

    integer(kind=i4b), dimension(:), intent(in) :: z
    integer(kind=i4b), intent(in)               :: N
    real(kind=qp)                               :: res

    real(kind=qp)     :: my_sum
    integer(kind=i4b) :: i

    if (modulo(N,2) == 1) then

      my_sum = 0.0_qp
      do i=0,N-1
        my_sum = my_sum + F( N, frac_part( (real(i, kind=qp)/N)*z ) )
      end do
      res = my_sum/N - 1

    else

      print *, "ERROR: R(z, N) is currently only working for odd N!"

    end if

  end function R


  ! This function is necessary for calculating the criterion R(z, N).
  !
  function F(N, x) result (res)

    integer(kind=i4b), intent(in)           :: N
    real(kind=qp), dimension(:), intent(in) :: x
    real(kind=qp)                           :: res

    real(kind=qp), dimension(size(x)) :: mysum
    integer(kind=i4b)                 :: h
    
    if (modulo(N,2) == 1) then ! Odd N

      mysum = 0.0_qp
      do h=1,(N-1)/2
        mysum = mysum + (cos(2*pi*h*x))/h
      end do
      res = product(2*mysum + 1)

    else ! Even N

      print *, "ERROR: F(N, x) for the calculation of R(z, N) is currently only working for odd N."

    end if

  end function F


end module mod_lattice_criteria
