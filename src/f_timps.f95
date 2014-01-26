! Module that implements Tim Pillard's testfunctions.
!
module f_timps

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_constants
  use mod_number_theory

  private
 
  public :: f_bivariate_normal_1p
  public :: f_bivariate_normal_np

  public :: f_timps01_1p
!  public :: f_timps01_exact
  public :: f_timps02_1p
  public :: f_timps03_np
  public :: f_timps04_np
  public :: f_timps05_np
  
  ! Functions based on the cubebased Gaussian exp(-max(x)^2)
  public :: f_timps_cube01_1p
  public :: f_timps_cube01_np
  public :: f_timps_cube01_exact_inf
  public :: f_timps_cube02_1p
  public :: f_timps_cube02_np
  public :: f_timps_cube02_exact_inf

  public :: f_timps_cap_1p
  public :: rho_one_np

contains


  ! The bivariate normal distribution.
  !
  ! Note: use this only when the dimension is 2!!!
  !
  function f_bivariate_normal_1p(x, params) result(rho)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: rho

    rho = 0.5_qp/pi*exp(-0.5_qp*(x(1)*x(1)+x(2)*x(2)))

  end function f_bivariate_normal_1p


  ! The bivariate normal distribution.
  !
  ! Note: use this only when the dimension is 2!!!
  !
  function f_bivariate_normal_np(x, params) result(rho)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: rho

    integer(kind=i4b) :: i

    do i=1,size(x,2)
      rho(i) = 0.5_qp/pi*exp(-0.5_qp*(x(1,i)*x(1,i)+x(2,i)*x(2,i)))
    end do

  end function f_bivariate_normal_np


  !                          s
  !                       --------'
  !                      '  |  |        1
  ! f(x[1], ..., x[s]) =    |  |    ---------
  !                         |  |        2
  !                         |  |    x[i]  + 1
  !                        i = 1
  !
  function f_timps01_1p(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    y = 1/(product(x*x+1.0_qp))

  end function f_timps01_1p


!  function f_timps01_exact(params, bounds) result (i)
!    type(functionparams), intent(in)    :: params
!    type(integrationbounds), intent(in) :: bounds
!    real(kind=qp)                       :: i
!
!    real(kind=qp), dimension(size(bounds%lb)) :: lb, ub
!
!    lb = bounds%lb
!    ub = bounds%ub
!
!    i = product( (exp(2*ub) - exp(2*lb) + lb - ub)/(exp(2.0_qp)-1) )
!
!  end function f_timps01_exact


  function f_timps02_1p(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    if (all(x<0.5_qp**(1.0_qp/size(x))) .and. all(x>-0.5_qp**(1.0_qp/size(x)))) then
       y = 1.0_qp/2.0_qp**(size(x)-1.0_qp)
    else
       y = 0.0_qp
    end if

  end function f_timps02_1p


  function f_timps03_np(x, params) result(f)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: f

    integer(kind=i4b) :: i

    do i=1,size(x,2)
      if (x(1,i) > x(2,i)) then
        f(i) = 1
      else
        f(i) = 0
      end if
    end do

  end function f_timps03_np


  ! A test function from Tim Pillards.
  !
  function f_timps04_np(x, params) result(f)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: f

    f = 1/(product(abs(x)+1, 1)**3)

  end function f_timps04_np


  ! A test function from Tim Pillards.
  !
  !     exp(-abs(prod(XY,2)))
  !
  function f_timps05_np(x, params) result (f)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: f

    f = exp(-product(abs(x)+1, 1))

  end function f_timps05_np


  ! Cube-based, more or less same value at border of a o-centered cube.
  !
  ! f = exp(-max(|x|)^2)
  !
  ! Convergence to one dimension as:
  !
  ! int_{R^d} f = int_{R+} exp(-|x1|^2)*|x1|^(d-1)/(d-1)!
  !
  ! exact value is 2^d GAMMA(d/2+1)
  !
  ! This is equal to 2^d (d/2)! for d even.
  !
  ! TODO: check of dit hetzelfde is als wat in artikel staat!!!
  !  
  function f_timps_cube01_1p(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    y = exp(-maxval(abs(x), dim=1)**2)
    
  end function f_timps_cube01_1p

  function f_timps_cube01_np(x, params) result (y)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: y

    y = exp(-maxval(abs(x), dim=1)**2)
    
  end function f_timps_cube01_np

  ! Return the exact integral value for integration over R^d of
  ! this testfunction.  The exact value is
  !
  !     2^d * GAMMA(d/2 + 1)
  !
  ! Note: exact values for uneven dimensions are not implemented yet!
  !
  function f_timps_cube01_exact_inf(d) result (res)
    integer(kind=i4b), intent(in) :: d
    real(kind=qp)                 :: res

    if (modulo(d,2) == 0) then
      res = 2**d*integer_factorial(d/2)
    else if (d == 15) then
      res = 459879458.192864464029960696693_qp ! computed with Maple
    else
      print *, "TODO: exact integral value for this dimension is not implemented yet!"
    end if
    
  end function f_timps_cube01_exact_inf


  ! Cube-based, more or less same value at border of a o-centered cube.
  !
  ! f = rho*g
  !
  ! with
  !
  !   rho = exp(-max(|x|)^2) and
  !   g = 1
  !
  ! Convergence to one dimension as:
  !
  ! int_{R^d} f = int_{R+} exp(-|x1|^2)*|x1|^(d-1)/(d-1)!
  !  
  ! TODO: check of dit hetzelfde is als wat in artikel staat!!!
  !
  function f_timps_cube02_1p(x, params) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    type(functionparams), intent(in)        :: params
    real(kind=qp)                           :: y

    y = exp(-maxval(x**2, dim=1))/product(1+abs(x), dim=1)
    
  end function f_timps_cube02_1p

  function f_timps_cube02_np(x, params) result (y)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: y

    y = exp(-maxval(x**2, dim=1))/product(1+abs(x), dim=1)
    
  end function f_timps_cube02_np


  ! Return the exact value for integration over R^d of
  ! this testfunction.  The exact value is
  !
  !   2^d * d * integral_{0}^{inf} ((ln(1+x))^(d-1)/(1+x))*exp(-x^2) dx
  !
  ! Currently, this routine returns some hardcoded values that
  ! were calculated with Maple using the command:
  !
  !   evalf[50](2^d*d*Int((ln(1+x))^(d-1)/(1+x)*exp(-x^2), x=0..infinity, method=_CCquad));
  !
  ! One could also have used '_Dexp', '_Gquad' or '_Sinc' as methods to compute
  ! these values using Maple.
  !
  function f_timps_cube02_exact_inf(d) result (res)
    integer(kind=i4b), intent(in) :: d
    real(kind=qp)                 :: res

    real(kind=qp), dimension(15), parameter :: exact = &
      (/ 1.2102673050066891634883422486284870827346823365527_qp, & ! d=1
         1.6999206629245814551143090404636506978456267246735_qp, & ! d=2
         2.6321607347218313738178489061830450227210009405594_qp, & ! d=3
         4.3787650664691875013523810506161681229118178403861_qp, & ! d=4
         7.7070230769572591335288147052477046999331383412394_qp, & ! d=5
         14.207563034653382941512489240242934279372405340754_qp, & ! d=6
         27.235786862297209788181140709621318235233029147997_qp, & ! d=7
         54.005393370409109547262944331167046373374623990061_qp, & ! d=8
         110.31399015434862789320142961361153956191761083819_qp, & ! d=9
         231.37060703330671975643736895050414093696293644387_qp, & ! d=10
         496.96297512018285292717485463765552855263934343782_qp, & ! d=11
         1090.7636692880249045304792353095301121186627486226_qp, & ! d=12
         2441.9247241993393396439814166906143772606735452871_qp, & ! d=13
         5567.3882104123291907797604438751188255375986300328_qp, & ! d=14
         12909.375249813430032126088212747385073893323787884_qp /) ! d=15

    res = exact(d)

  end function f_timps_cube02_exact_inf


  function f_timps_cap_1p(x) result (y)
    real(kind=qp), dimension(:), intent(in) :: x
    real(kind=qp)                           :: y

    if (all(x<1.0_qp) .and. all(x>-1.0_qp)) then
       y = 1.0_qp
    else
       y = 0.0_qp
    end if

  end function f_timps_cap_1p


  function rho_one_np(x, params) result(rho)
    real(kind=qp), dimension(:,:), intent(in) :: x
    type(functionparams), intent(in)          :: params
    real(kind=qp), dimension(size(x,2))       :: rho

    rho = 1.0_qp

  end function rho_one_np

end module f_timps
