! Module that implements some of Ronald's simple and quick testfunctions
! (obtained via private communication).

module f_ronald

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_constants
  use mod_utilities
  use mod_bernoulli
  use mod_number_theory

  private

  public :: f_ronald01_1p
  public :: f_ronald01_exact
  public :: f_ronald02_1p
  public :: f_ronald02_exact
  public :: f_ronald03_1p
  public :: f_ronald03_exact
  public :: f_ronald04_1p
  public :: f_ronald04_exact
  public :: f_ronald05_1p
  public :: f_ronald05_exact
  public :: f_ronald_polynom6_1p
  public :: f_ronald_polynom6_bernoulli5_1p
  public :: f_ronald_polynom6_exact
  public :: f_sine_power_1p
  public :: f_sine_power_deriv
  public :: f_sine_power_bernoulli5_1p
  public :: f_sine_power_exact
  public :: f_sine_power_1d_exact
  public :: f_dtai_1p
  public :: generate_dtai_parameters

  contains

    !                   2
    !       f(x) := x[1]
    !
    function f_ronald01_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = x(1)**2

    end function f_ronald01_1p

    function f_ronald01_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = (ub(1)**3-lb(1)**3)*product(ub(2:)-lb(2:))/3

    end function f_ronald01_exact


    !                  s
    !               --------'
    !              '  |  |        2
    !      f(x) :=    |  |    x[i]
    !                 |  |
    !                 |  |
    !                i = 1
    !
    function f_ronald02_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(x*x)

    end function f_ronald02_1p

    function f_ronald02_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product((ub**3-lb**3)/3)

    end function f_ronald02_exact


    !                  s
    !               --------'
    !              '  |  |        6
    !      f(x) :=    |  |    x[i]
    !                 |  |
    !                 |  |
    !                i = 1
    !
    function f_ronald_polynom6_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(x**6)

    end function f_ronald_polynom6_1p

    !                  s
    !               --------'
    !              '  |  |   
    !      f(x) :=    |  |  ... some lengthy polynomial ...  
    !                 |  |
    !                 |  |
    !                i = 1
    !
    function f_ronald_polynom6_bernoulli5_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product((x**2*(x**2*(x*(6*x-18)+15)-3)+1)/6.0)

    end function f_ronald_polynom6_bernoulli5_1p

    function f_ronald_polynom6_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product((ub**7-lb**7)/7)

    end function f_ronald_polynom6_exact


    !                  s
    !               --------'
    !              '  |  |
    !      f(x) :=    |  |    sin(2 pi x)
    !                 |  |
    !                 |  |
    !                i = 1
    !
    function f_ronald03_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(sin(2*pi*x))

    end function f_ronald03_1p

    function f_ronald03_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( (cos(2*pi*lb)-cos(2*pi*ub))/(2*pi) )

    end function f_ronald03_exact


    !                  s
    !               --------'
    !              '  |  |
    !      f(x) :=    |  |    cos(2 pi x)
    !                 |  |
    !                 |  |
    !                i = 1
    !
    function f_ronald04_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(cos(2*pi*x))

    end function f_ronald04_1p

    function f_ronald04_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( (sin(2*pi*ub)-sin(2*pi*lb))/(2*pi) )

    end function f_ronald04_exact


    !              s
    !           --------'
    !          '  |  |
    !  f(x) :=    |  |    (1 + sin(2 Pi x))
    !             |  |
    !             |  |
    !            i = 1
    !
    function f_ronald05_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(1 + sin(2*pi*x))

    end function f_ronald05_1p

    function f_ronald05_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( ub - lb + (cos(2*pi*lb)-cos(2*pi*ub))/(2*pi) )

    end function f_ronald05_exact


    !                      s
    !                    /===\            k
    !                     ! !      i pi x  i
    !                     ! !  sin(------)
    !                     ! !        4
    !                    i = 1
    !
    function f_sine_power_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      integer(kind=i4b) :: i

      y = product( (/ ((sin((pi*i*x(i))/4))**params%k(i), i=1,size(x)) /) )

    end function f_sine_power_1p

    function f_sine_power_exact(params, bounds) result (res)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: res

      integer(kind=i4b) :: i

      res = 1.0_qp
      do i=1,size(bounds%lb)
        res = res*f_sine_power_1d_exact(i, params%k(i), &
                                        bounds%lb(i), bounds%ub(i))
      end do
      
    end function f_sine_power_exact


    ! Return the exact value of the one-dimensional integral
    !                     b
    !                    /
    !                    [     k %pi i x
    !                    I  sin (-------) dx
    !                    ]          4
    !                    /
    !                     a
    !
    ! with k smaller than 6.
    !
    function f_sine_power_1d_exact(i, k, a, b) result (res)
      integer(kind=i4b), intent(in) :: i
      integer(kind=i4b), intent(in) :: k
      real(kind=qp), intent(in)     :: a
      real(kind=qp), intent(in)     :: b

      real(kind=qp) :: res

      select case (k)

        case (1)

          res = 4*cos(pi*a*i/4.0)/(pi*i)-4*cos(pi*b*i/4.0)/(pi*i)

        case (2)

          res = &
                (2*sin(pi*a*i/2.0)-pi*a*i)/(pi*i)/2.0-(2*sin(pi*b*i/ &
                 2.0)-pi*b*i)/(pi*i)/2.0

        case (3)

          res = &
                cos(3.0*pi*b*i/4.0)/(pi*i)/3.0-3*cos(pi*b*i/4.0)/ &
        (pi*i)-cos(3.0*pi*a*i/4.0)/(pi*i)/3.0+3*cos(pi*a*i/4.0)/ &
        (pi*i)

        case (4)
          
          res = &
(sin(pi*b*i)-8*sin(pi*b*i/2.0)+3*pi*b*i)/(pi*i)/ &
     8.0-(sin(pi*a*i)-8*sin(pi*a*i/2.0)+3*pi*a*i)/(pi*i)/8.0

        case (5)

          res = &
-cos(5.0*pi*b*i/4.0)/(pi*i)/20.0+5.0*cos(3.0*pi*b*i/4.0)/ &
     (12.0*pi*i)+(-5.0)*cos(pi*b*i/4.0)/ &
     (2.0*pi*i)+cos(5.0*pi*a*i/4.0)/(pi*i)/ &
     20.0+(-5.0)*cos(3.0*pi*a*i/4.0)/ &
     (12.0*pi*i)+5.0*cos(pi*a*i/4.0)/(2.0*pi*i)

        case (6)

          res = &
(sin(3.0*pi*a*i/2.0)-9*sin(pi*a*i)+45*sin(pi*a*i/ &
     2.0)-15*pi*a*i)/(pi*i)/48.0-(sin(3.0*pi*b*i/ &
     2.0)-9*sin(pi*b*i)+45*sin(pi*b*i/2.0)-15*pi*b*i)/(pi*i)/ &
     48.0

        case default

          print *, "ERROR: cannot return exact value for integral for k>6"

      end select

    end function f_sine_power_1d_exact

    ! Return the mth derivative of f_sine_power in one dimension
    ! with parameter j in the argument of the sine and power k:
    !
    !                    m
    !                   d       k %pi j x
    !                   --- (sin (-------))
    !                     m          4
    !                   dx
    !
    !
    function f_sine_power_deriv(x, m, j, k) result (y)
      real(kind=qp), intent(in)     :: x
      integer(kind=i4b), intent(in) :: m ! the order of the derivative
      integer(kind=i4b), intent(in) :: j ! the parameter in the argument of
                                         ! the sine
      integer(kind=i4b), intent(in) :: k ! the power of the sine
      real(kind=qp)                 :: y

      integer(kind=i4b) :: i

      select case (modulo(k,2))

        ! even
        case (0)
           y = (pi*j/4)**m * (-1)**(k/2)/(2**(k-1)) &
                 * sum( (/ ( (-1)**i*nchoosek(k,i)*(k-2*i)**m &
                             *cos((k-2*i)*pi*j*x/4+m*pi/2) , i=0,(k-2)/2 ) /) )

        ! odd
        case (1)
           y = (pi*j/4)**m * (-1)**((k-1)/2)/(2**(k-1)) &
                 * sum( (/ ( (-1)**i*nchoosek(k,i)*(k-2*i)**m &
                              *sin((k-2*i)*pi*j*x/4+m*pi/2), i=0,(k-1)/2 ) /) )

      end select

    end function f_sine_power_deriv

    ! The Bernoulli-periodized version of f_sine_power
    ! for alpha=5.
    !
    function f_sine_power_bernoulli5_1p(x, params) result(y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y
      
      integer(kind=i4b) :: i, j
      integer(kind=i4b) :: alpha

      alpha = 5

      y = product((/((sin(pi*i*x(i)/4))**params%k(i) &
                     + sum((/ (bernoulli_poly(j+1,x(i))/integer_factorial(j+1) &
                              *(f_sine_power_deriv(0.0_qp, j, i, params%k(i)) &
                            - f_sine_power_deriv(1.0_qp, j, i, params%k(i))), &
                           j=0,alpha-1 ) /)), &
                  i=1,size(x))/))

    end function f_sine_power_bernoulli5_1p


    ! Ronald's DTAI-function.
    !
    !    f(w) := max(0, w x  + b , ..., w x  + b )
    !                      1    1          n    n
    !
    ! The vectors \vec{x}_1,...,\vec{x}_n must be stored as the columns
    ! of params%x.
    !
    function f_dtai_1p(x, params) result(y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      ! Method 1
      !y = maxval( (/ 0.0_qp, matmul(x,params%x) + params%b /) )

      ! Method 2
      !y = maxval (matmul(x,params%x) + params%b)
      !y = max (0.0_qp, y)

      ! Method 3 (probably fastest)
      y = max(0.0_qp, maxval(matmul(x,params%x) + params%b))

    end function f_dtai_1p

    ! Generate parameters x and b for the DTAI function.
    !
    ! \vec{x}_1,...,\vec{x}_n are in [-1,1] and stored as the columns
    ! of x.
    !
    ! The scalars b_1,...,b_n are in the interval [-sqrt(n),sqrt(n)].
    !
    subroutine generate_dtai_parameters(n, x, b)
      integer(kind=i4b), intent(in)              :: n
      real(kind=qp), dimension(:,:), intent(out) :: x
      real(kind=qp), dimension(:), intent(out)   :: b

      call random_number(x)
      x = 2*x-1

      call random_number(b)
      b = (2*b-1)*sqrt(real(n, kind=qp))

    end subroutine generate_dtai_parameters


end module f_ronald
