! Implements testfunctions found in different papers.
!
! References:
!
!   [1] ALGORITHM 647: `Implementation and Relative Efficiency of Quasirandom
!       Sequence Generators', Fox, L. Bennett, ACM Transactions on
!       Mathematical Software, Vol. 12, nr 4, 1986, pages 362--376.
!
!   [2] `Generating and Testing the Modified Halton Sequences', Emanouil I.
!       Atanassov and Mariya K. Durchova in `Revised Papers from the 5th
!       International Conference on Numerical Methods and Applications',
!       Lecture Notes in Computer Science 2542, pp. 91-98, 2003.
!
!   [3] `One more experiment on estimating high-dimensional integrals by
!       quasi-Monte Carlo methods', Sobol, I. M. and Asotsky, D. I.,
!       Mathematics and Computers in Simulation, Vol. 62, nr 3-6, 2003,
!       pages 255-263.
!
!   [4] `A Method for Numerical Integration', C. B. Haselgrove, Mathematics
!       of Computation, Vol. 15, 1961, pages 323-337.
!
! TODO:
!
!   * Implement f_haselgrove61_exact for general hyperrectangles instead of
!     just for the unit hypercube.

module f_testfunctions

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_special_functions

  private

  public :: f_simple_product_1p
  public :: f_simple_product_exact
  public :: f_haselgrove61_1p
  public :: f_haselgrove61_exact
  public :: f_atanassov03_f2_1p
  public :: f_atanassov03_f2_exact
  public :: f_atanassov03_f3_1p
  public :: f_atanassov03_f3_exact
  public :: f_bart_george_f01_1p
  public :: f_bart_george_f01_exact
  public :: f_fox86_1p
  public :: f_fox86_np
  public :: f_fox86_exact
  public :: f_sobol03asotsky_1p
  public :: f_sobol03asotsky_exact
  public :: f_wang_hickernell

  contains


    ! Mainly used for debugging purposes.
    !
    function f_simple_product_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(x)

    end function f_simple_product_1p

    function f_simple_product_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb
      real(kind=qp), dimension(size(bounds%ub)) :: ub 

      lb = bounds%lb
      ub = bounds%ub

      i = product(0.5_qp*(ub**2-lb**2))

    end function f_simple_product_exact


    ! See [4]
    function f_haselgrove61_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = exp(-product(x))
      
    end function f_haselgrove61_1p

    ! The exact value of this integral is:
    !
    !                            /    s                               \
    !                            | --------'     (k + 1)       (k + 1)|
    !                          k |'  |  |    b[r]        - a[r]       |
    !                      (-1)  |   |  |    -------------------------|
    !             infinity       |   |  |              k + 1          |
    !              -----         |   |  |                             |
    !               \            \  r = 1                             /
    !                )     --------------------------------------------
    !               /                           k!
    !              -----
    !              k = 0
    !
    function f_haselgrove61_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      integer(kind=i4b)                         :: k
      real(kind=qp)                             :: up, term

      print *, "WARNING: beware of rounding errors when using the exact"
      print *, "value of f_haselgrove61 for an integral over a general"
      print *, "rectangle."

      term = product(bounds%ub-bounds%lb)

      ! up is equal to (-1)**k / k!
      up = 1.0_qp
      i = 0.0_qp
      k = 0
      do

        if (abs(term) < epsilon(1.0_qp)) then
          exit
        else
          k = k+1
          i = i + term

          ! Calculate next term to be added
          up = -up / k
          term = up * product( (bounds%ub**(k+1)-bounds%lb**(k+1))/(k+1) )

       end if

     end do

    end function f_haselgrove61_exact


    function f_atanassov03_f2_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(x**3 + 0.75)
      
    end function f_atanassov03_f2_1p

    function f_atanassov03_f2_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      ! This function doesn't use parameters
      
      real(kind=qp), dimension(size(bounds%lb)) :: lb
      real(kind=qp), dimension(size(bounds%ub)) :: ub

      lb = bounds%lb
      ub = bounds%ub

      i = product((ub**4+3*ub-(lb**4+3*lb))*0.25_qp)
      
    end function f_atanassov03_f2_exact


    function f_atanassov03_f3_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      real(kind=qp)                           :: myprod
      integer(kind=i4b)                       :: j

      y = 0.0_qp

      ! NOTE: according to the paper, this should be -1, but the implementation
      ! of Atanassov and the exact value for the integral indicate that we
      ! should use -1 here in order to get a positive result.
      myprod = -1.0_qp

      do j=1, size(x)
        myprod = -myprod*x(j)
        y = y + myprod
      end do

    end function f_atanassov03_f3_1p

    function f_atanassov03_f3_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      integer(kind=i4b) :: dimen

      ! No parameters appear in the integration result here.

      ! TODO: check the value of the integral when integrating
      ! over other ranges than the unit cube.

      dimen = size(bounds%ub)

      i = 1.0_qp/3*(1.0_qp-(-1.0_qp)**dimen/2.0_qp**dimen) 

    end function f_atanassov03_f3_exact


    function f_bart_george_f01_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y
      
      real(kind=qp), dimension(size(x)) :: a, u

      a = params%a
      u = params%u

      y = product(1-a*abs(x-u))
      
    end function f_bart_george_f01_1p

    function f_bart_george_f01_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(params%a))  :: a, u
      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      a = params%a
      u = params%u
      lb = bounds%lb
      ub = bounds%ub

      i = product(-a*u*u-lb-a*u*lb+0.5_qp*a*lb*lb+ub+a*u*ub-0.5_qp*a*ub*ub)

    end function f_bart_george_f01_exact


    function f_fox86_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      real(kind=qp), dimension(size(x)) :: a, u

      a = params%a
      u = params%u

      y = product(abs(a*x-u))
      
    end function f_fox86_1p

    function f_fox86_np(x, params) result (y)
      real(kind=qp), dimension(:,:), intent(in) :: x
      type(functionparams), intent(in)          :: params
      real(kind=qp), dimension(size(x, dim=2))  :: y

      real(kind=qp), dimension(size(x, dim=1)) :: a, u
      integer(kind=i4b)                        :: N

      N = size(x, dim=2)
      a = params%a
      u = params%u

      y = product(abs(spread(a,2,N)*x-spread(u,2,N)), dim=1)
    
    end function f_fox86_np

    ! The exact value of the integral in [lb,ub].
    ! Assumptions:
    !   lb <= u/a <= ub
    function f_fox86_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp) :: i

      real(kind=qp), dimension(size(bounds%lb)) :: a, u
      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      a = params%a
      u = params%u
      lb = bounds%lb
      ub = bounds%ub

      ! TODO: check if this integral value is correct.
      ! TODO: check why changing the value of i doesn't have
      !       any influence on wether pub_mascagni_chi reports
      !       successfull or not.
      i = product(0.5_dp*(2*u*u+lb*lb*a*a-2*u*lb*a+ub*ub*a*a-2*u*ub*a)/a)
      !i = 1

    end function f_fox86_exact


    ! See [3]
    ! If params%c are all taken to be the same value c, then
    ! if the product dimension*c is large, the function has a
    ! sharp peak at the point (1,...,1).
    !
    function f_sobol03asotsky_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x 
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product(1+params%c*(x-0.5_qp))

    end function f_sobol03asotsky_1p

    !
    function f_sobol03asotsky_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      ! The exact value in [0,1]^s is 1.
      i = 1.0_qp

      ! TODO: check what the exact value is for other bounds than the
      !       unit cube!!!

    end function f_sobol03asotsky_exact


    function f_wang_hickernell(x, params) result (y)
      real(kind=qp), dimension(:,:), intent(in) :: x
      type(functionparams), intent(in)          :: params
      real(kind=qp), dimension(size(x, dim=2))  :: y

      real(kind=qp), dimension(size(x, dim=1))  :: a
      integer(kind=i4b)                         :: N

      N = size(x, dim=2)
      a = params%a

      y = product( (abs(4*x-2) + spread(a,2,N))/(1+spread(a,2,N)), dim=1 )

    end function f_wang_hickernell


end module f_testfunctions
