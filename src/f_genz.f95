! Module that implements Alan Genz's testfunctions.
!
! References:
!
!   [1] 'High-Dimensional Numerical Integration on Parallel Computers',
!       Master's thesis Rudolf Schuerer, University of Salzburg, may 2001.
!
!   [2] 'Testing multidimensional integration routines', Alan Genz, Tools,
!       methods and languages for scientific and engineering computation,
!       Elsevier North-Holland, 1984, ISBN 0-444-87570-0.

module f_genz

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_special_functions
  use mod_constants

  private

  public :: f_oscillatory_1p      ! Genz's function 1
  public :: f_oscillatory_exact
  public :: f_productpeak_1p      ! Genz's function 2
  public :: f_productpeak_exact
!  public :: f_cornerpeak_1p       ! Genz's function 3
!  public :: f_cornerpeak_exact
  public :: f_gaussian_1p         ! Genz's function 4
  public :: f_gaussian_exact
  public :: f_c0_1p               ! Genz's function 5
  public :: f_c0_np
  public :: f_c0_exact
  public :: f_discontinuous_1p    ! Genz's function 6
  public :: f_discontinuous_exact

  contains


    ! Genz function 1: OSCILLATORY
    !
    ! Note:
    !
    !  * We follow the definition from [1] and [2].  The function implemented
    !    in Genz's code is slightly different (all parameters a are taken to
    !    be 1).
    !
    function f_oscillatory_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      real(kind=qp), dimension(size(x)) :: a, u

      a = params%a
      u = params%u

      y = cos(2*pi*u(1)+dot_product(a, x))

    end function f_oscillatory_1p

    ! See [1] for a derivation of this exact result.
    ! Note:
    !
    !   * This is the exact value for the integral of the function as defined
    !     in [1] and [2].  Genz's implementation contains code that gives the
    !     exact value for a slighty different function (probably where all
    !     parameters a are taken to be 1).
    !
    function f_oscillatory_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: a, u
      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      a = params%a
      u = params%u
      lb = bounds%lb
      ub = bounds%ub

      ! TODO: this is the expression from Rudolf Schuerer's thesis.  Check
      !       if it results in the same value as Genz's code for the exact
      !       value.  It probably won't, since Genz implemented a slighty
      !       other function...
      i = cos( 2*pi*u(1)+sum(0.5*a*(ub-lb)) )*product( 2/a*sin(0.5*a*(ub-lb)) )

    end function f_oscillatory_exact


    ! Genz function 2: PRODUCT PEAK
    !
    ! Note:
    !
    !   * For dimensions > 50, the numbers arising during the approximation of
    !     the integral by most algorithms become too small to be represented
    !     by floating point numbers.  Thus, this integrand family is not used
    !     for dimensions > 40.
    !
    function f_productpeak_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      real(kind=qp), dimension(size(x)) :: a, u

      a = params%a
      u = params%u

      y = product( 1/( a**(-2)+(x-u)**2 ) )

    end function f_productpeak_1p

    function f_productpeak_exact(params, bounds) result(i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: a, u
      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      a = params%a
      u = params%u
      lb = bounds%lb
      ub = bounds%ub

      i = product( a*( atan(a*(ub-u))-atan(a*(lb-u)) ) )

    end function f_productpeak_exact


    ! Genz function 4: GAUSSIAN
    !
    function f_gaussian_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = exp(-sum(params%a*params%a*(x-params%u)*(x-params%u)))

    end function f_gaussian_1p

    ! See code from Alan Genz on how we obtain this exact value,
    ! or the expression in Rudolf Schuerer's master's thesis.
    ! Note: Rudolf's thesis contains the expression with erf(x)
    ! while the code from Alan Genz uses phi(x) (the Gaussian probability
    ! distribution!)
    !
    function f_gaussian_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: a, u
      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub
      integer(kind=i4b)                         :: j

      a = params%a
      u = params%u
      lb = bounds%lb
      ub = bounds%ub

      i = 1.0_qp

      do j=1, size(a)
        i = i*( sqrt(pi)/a(j) )                          &
            * (  phi( a(j)*(ub(j)-u(j))*sqrt(2.0_qp) )   &
               - phi( a(j)*(lb(j)-u(j))*sqrt(2.0_qp) ) )
      end do

    end function f_gaussian_exact


    ! Genz function 5: C^0 FUNCTION
    !
    function f_c0_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = exp(-sum(params%a*abs(x-params%u)))

    end function f_c0_1p

    ! IMPORTANT: the rows for x should represent the different dimensions!!!
    function f_c0_np(x, params) result (y)
      real(kind=qp), dimension(:,:), intent(in) :: x
      type(functionparams), intent(in)          :: params
      real(kind=qp), dimension(size(x, dim=2))  :: y

      real(kind=qp), dimension(size(x, dim=1))  :: a, u
      integer(kind=i4b)                         :: N

      N = size(x, dim=2)
      a = params%a
      u = params%u

      y = exp(-sum( spread(a,2,N)*abs(x-spread(u,2,N)), dim=1 ))
    
    end function f_c0_np

    function f_c0_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds

      real(kind=qp), dimension(size(bounds%lb)) :: a, u
      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub
      real(kind=qp)                             :: i

      a = params%a
      u = params%u
      lb = bounds%lb
      ub = bounds%ub

      i = product((2.0_qp-exp(a*(lb-u))-exp(a*(u-ub)))/a)

    end function f_c0_exact


    ! Genz's function 6: DISCONTINUOUS
    !
    ! Note:
    !
    !   * Genz's definition sets this function to zero if
    !
    !			x_1 > u_1 OR x_2 > u_2
    !
    !     Looking at his implementation however, it is doubtfull that this
    !     condition is checked *only* for the first two dimensions.  It appears
    !     as it is being checked for all dimensions...
    !     This version only checks the first 2 dimensions and follows Genz's
    !     mathematical definition, not his implementation.
    !
    function f_discontinuous_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = exp(dot_product(params%a, x))

      ! Check only the first 2 dimensions, not all!
      if ( x(1)>params%u(1) .or. x(2)>params%u(2)) then
        y = 0
      end if

    end function f_discontinuous_1p

    ! TODO:
    !
    !   * Check if the value this function returns is indeed the exact value.
    !
    function f_discontinuous_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: a, u
      real(kind=qp), dimension(size(bounds%lb)) :: ub_prime
      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      a = params%a
      u = params%u
      lb = bounds%lb
      ub = bounds%ub
      ub_prime = ub

      ! First set the new values for the upper bounds for the first two
      ! dimensions (See [1], page 30-31).
      ! Note that we do *not* set it for *all* dimensions because we follow the
      ! exact mathematical definition of Genz's papers.
      ub_prime(1) = max( lb(1), min(ub(1),u(1)) )
      ub_prime(2) = max( lb(2), min(ub(2),u(2)) )

      i = product( ( exp(a*ub_prime) - exp(a*lb) )/a )

    end function f_discontinuous_exact

end module f_genz
