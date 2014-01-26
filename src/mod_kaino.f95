! This module implements Kaino's method for multidimentional numerical
! integration of certain periodic integrands.  It has a theoretical
! convergence order of O(1/(N^{q+1})).
!
! Note that this method is a closed integration method, so an increase from
! I_N to I_{N+1} is not possible.
!
! Notes:
!   * In his paper, Kaino makes an assumption on the irrationals, namely
!     that
!               1/N << irrationals << 1
!
! References:
!
!   [1] `Another Note on Haselgrove's Method for Numerical Integration',
!       Kaino K., Journal of the Korean Physical Society, Vol. 40, No. 6,
!       June 2002, pp. 1010-1014.

module mod_kaino

  use numeric_kinds
  use mod_function
  use mod_constants
  use mod_number_theory
  use mod_utilities
  use mod_periodize

  private

  public :: kaino

contains


  ! Kaino's integration method.  Note that we
  ! only return one value for the integral at N and not the
  ! incremental values for each N because this is a closed method.
  !
  subroutine kaino(q, N, irrationals, func, params, periodizer, integral_value)
    integer(kind=i4b), intent(in)           :: q
    integer(kind=i4b), intent(in)           :: N
    real(kind=qp), intent(in), dimension(:) :: irrationals
    interface
      function func(x, fparams) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function func
    end interface
    type(functionparams), intent(in)        :: params
    character(len=*), intent(in)            :: periodizer
    real(kind=qp), intent(out)              :: integral_value

    integer(kind=i4b) :: k
    real(kind=qp)     :: f
    real(kind=qp)     :: my_sum
    real(kind=qp)     :: a_q

    ! First calculate the normalizing constant for the weight function
    if (modulo(q,2) == 0) then
      a_q = real(double_factorial(q), kind=qp)/double_factorial(q-1)
    else
      a_q = ((pi/2)*double_factorial(q))/double_factorial(q-1)
    end if

    ! Now do the summation
    my_sum = 0.0_qp
    do k=1,N
      call periform(func, params, frac_part(k*irrationals), periodizer, f)
      my_sum = my_sum + a_q*( (sin((pi*k)/N))**q )*f
    end do
    integral_value = (my_sum/N)

  end subroutine kaino


end module mod_kaino
