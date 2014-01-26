! Module that implements Mathe and Wei's algorithm for Quasi-Monte Carlo
! integration over R^d.
!
! Usage:
!
!   * First, initialize the method with
!
!       init_mathe_wei(init_V, init_w)
!
!     where V is the initial width of the hypercubes and w is the
!     growth rate of the widths of the hypercubes.
!       
!   * Secondly, you have to calculate the maximum size of the low-discrepancy
!     point-set that the algorithm will need, using the subroutine
!
!       calculate_max_pointset_size(N, d, weight_type, s, res)
!
!     where
!
!       N               required number of function evaluations
!       d               dimension
!       weight_type     the type of weight function, which can be one of
!                         'MVN' : Multivariate Normal Distributions
!                         'ECD' : Elliptically Contoured Distributions
!                         'M_SR': the most general class for the weight
!                                 function rho
!       s               decay parameter (provide a dummy value if this is
!                       not applicable)
!       res             the result variable that will contain the maximum size
!
!   * After that, you can calculate an approximation for an integral using
!
!       mathe_wei(N, x, f, params_f, rho, params_rho, weight_type, s, res)
!
!     where
!
!       x               a low-discrepancy pointset that is large enough and
!                       which has its different points specified as its columns.
!       f               the function
!       params_f        some parameters for the function (if they exist)
!       rho             the weight function
!       params_rho      some parameters for the weight function (if they exist)
!       res             the resulting integral value
!       
! Notes:
!
!       * The general formula for the ECD and M_SR case for m are confirmed
!         by Peter Mathe.
!
!       * Mathe also confirms that in the ECD and M_SR case, we may replace
!         the 2 by a certain w in the formula for the n_j:
!
!           "The initial volume V appears, but due to the construction of the
!           $m_j$ this cancels, and one simply has to replace 2 by w.  BTW,
!           this is dicussed in §~5, around (23), when using Fibonacci numbers."
!
! References:
!
!   [1] `Quasi-Monte Carlo integration over R^d', Peter Mathe and Gang Wei,
!        Mathematics of Computation, Volume 73, Number 246, Pages 827-841, 2004.
!
! TODO:
!
!   * Check how the formulas for MVN look like with general V and w.
!
!   * Instead of passing the points x of the algorithm in one block (see
!     the mathe_wei(...) subroutine, check if we can do it in a
!     point-by-point way and use the f_1p() function evaluation methods.
!
module mod_mathewei

  use numeric_kinds
  use mod_function

  private

  public :: init_mathe_wei
  public :: calculate_max_pointset_size
  public :: mathe_wei

  private :: calculate_m
  private :: calculate_nj
  private :: calculate_qj
  private :: evaluate_s_mn
  private :: indicator_I

  real(kind=qp), private :: V = 2.0_qp ! initial width of a hypercube
                                       ! (default is 2)
  real(kind=qp), private :: w = 2.0_qp ! growth of the widths of the cubes
                                       ! (default is 2)


contains

  ! Initialize Mathe and Wei's algorithm with a certain initial width V of a
  ! hypercube and a certain initial growth w of the widths of the cubes.  If
  ! this initialization is forgotten or not done, then by default V=2 and w=2.
  !
  subroutine init_mathe_wei(init_V, init_w)
    real(kind=qp), intent(in) :: init_V
    real(kind=qp), intent(in) :: init_w

    V = init_V
    w = init_w

  end subroutine init_mathe_wei


  ! Determine the maximum number of low-discrepancy points over [0,1] that we
  ! will need for this algorithm and return it in the argument res.
  !
  subroutine calculate_max_pointset_size(N, d, weight_type, s, res)
    integer(kind=i4b), intent(in)  :: N
    integer(kind=i4b), intent(in)  :: d
    character(len=*), intent(in)   :: weight_type
    integer(kind=i4b), intent(in)  :: s
    integer(kind=i4b), intent(out) :: res

    integer(kind=i4b)                            :: m
    integer(kind=i4b), dimension(:), allocatable :: nj

    call calculate_m(N, d, s, weight_type, m)

    allocate(nj(0:m))
    call calculate_nj(N, d, s, m, weight_type, nj)
    res = maxval(nj)
    deallocate(nj)

  end subroutine calculate_max_pointset_size


  ! Run Mathe and Wei's algorithm and use x as low-discrepancy points.  The
  ! columns of x represent the different points, the rows of x represent the
  ! different dimensions.
  ! TODO:
  !   * Instead of passing the points of the algorithm in one block x,
  !     check if we can do it in a point-by-point way and use the f_1p()
  !     function evaluation methods.
  !
  subroutine mathe_wei(N, x, f, params_f, rho, params_rho, weight_type, s, res)
    integer(kind=i4b), intent(in)             :: N
    real(kind=qp), dimension(:,:), intent(in) :: x
    interface
      function f(x, params_f) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:,:), intent(in) :: x
        type(functionparams), intent(in)          :: params_f
        real(kind=qp), dimension(size(x,2))       :: y             ! TODO: CHECK IF THIS SHOULD BE size(x,1) or size(x,2)
      end function f
    end interface
    type(functionparams), intent(in)          :: params_f
    interface
      function rho(x, params_rho) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:,:), intent(in) :: x
        type(functionparams), intent(in)          :: params_rho
        real(kind=qp), dimension(size(x,2))       :: y             ! TODO: CHECK IF THIS SHOULD BE size(x,1) or size(x,2)
      end function rho
    end interface
    type(functionparams), intent(in)          :: params_rho
    character(len=*), intent(in)              :: weight_type
    integer(kind=i4b), intent(in)             :: s
    real(kind=qp), intent(out)                :: res


    integer(kind=i4b)                            :: m, d
    integer(kind=i4b), dimension(:), allocatable :: nj
    real(kind=qp), dimension(:), allocatable     :: qj

    print *, "N = requested number of function evaluations = ", N
    d = size(x, dim=1)

    call calculate_m(N, d, s, weight_type, m)
    print *, "m = ", m
    allocate(nj(0:m), qj(0:m))

    call calculate_nj(N, d, s, m, weight_type, nj)
    call calculate_qj(m, d, qj)
    print *, "nj = ", nj
    print *, "sum(nj) = total amount of function evaluations = ", sum(nj)
    print *, "qj = Volume(Q_j) = ", qj

    call evaluate_s_mn(m, nj, f, params_f, rho, params_rho, qj, x, res)
    
    deallocate(nj, qj)

  end subroutine mathe_wei


  ! Calculate m, which determines the number of cubes to use in the algorithm.
  ! This number depends on the total number N of points to be used, the
  ! dimension d and the speed s at which the weight function decays.  See also
  ! Theorem 2 on page 834 and Corollary 2 on page 836 in [1].
  !
  subroutine calculate_m(N, d, s, weight_type, m)
    integer(kind=i4b), intent(in)  :: N
    integer(kind=i4b), intent(in)  :: d ! dimension
    integer(kind=i4b), intent(in)  :: s ! decay-factor
    character(len=*), intent(in)   :: weight_type
    integer(kind=i4b), intent(out) :: m


    select case (weight_type)

      case ("MVN") ! Multivariate Normal

        ! Due to the rapid decay of this weight, fewer cubes are necessary in
        ! the hierarchy.  We may even choose to use only one hypercube to apply
        ! the QMC algorithm.  See Section 4.1 page 836 in [1].
        m = floor(0.5_qp*log(log(real(N, kind=qp))))

      case ("ECD") ! Elliptically Contoured Distributions

        ! The class of Elliptically Contoured Distributions is a subset of
        ! M_s(R) for R large enough.  For this class, we have that s>d.  (See
        ! also page 835 in [1]).
        if (s>d) then
          m = ceiling(1.0_qp/(s-d)*log(real(N, kind=qp))/log(w))
        else
          stop "ERROR: s must be larger than d for Elliptically Contoured Distributions!"
        end if

      case ("M_SR")
      
        ! This is the most general case of weights with prescribed decay at
        ! infinity (See Section 3.2 in [1]).
        m = ceiling(1.0_qp/(s-d)*log(real(N, kind=qp))/log(w))

      case default

        stop "ERROR: unknown weight function type!"

      end select

  end subroutine calculate_m


  ! Calculate all the volumes for Q_0 up to Q_m with dimension d.
  ! See page 828 in [1].
  !
  subroutine calculate_qj(m, d, q)
    integer(kind=i4b), intent(in)             :: m
    integer(kind=i4b), intent(in)             :: d
    real(kind=qp), dimension(0:), intent(out) :: q

    integer(kind=i4b) :: j

    q = (/ ( (V*w**j)**d, j=0,m ) /)

  end subroutine calculate_qj


  ! Calculate the number of LDS-points n_j for each cube.  This number depends
  ! on the total preferred number N of points to be used, the dimension d and
  ! the speed s at which the weight function decays.  See also Theorem 2 in [1].
  !
  subroutine calculate_nj(N, d, s, m, weight_type, nj)
    integer(kind=i4b), intent(in)                :: N
    integer(kind=i4b), intent(in)                :: d ! dimension
    integer(kind=i4b), intent(in)                :: s ! decay factor
    integer(kind=i4b), intent(in)                :: m
    character(len=*), intent(in)                 :: weight_type
    integer(kind=i4b), dimension(:), intent(out) :: nj

    integer(kind=i4b)             :: j
    real(kind=qp), dimension(0:m) :: dummy

    select case (weight_type)

      case ("MVN") ! Multivariate Normal

        ! See Section 4.1 page 836 in [1].
        dummy = (/ (exp(-2**(2.0_qp*j-1)), j=0,m) /)
        nj = ceiling(N*dummy/sum(dummy))

      case ("ECD") ! Elliptically Contoured Distributions
        
        ! The class of Elliptically Contoured Distributions is a subset
        ! of M_s(R) for R large enough.  For this class, we have that s>d.
        ! (See also page 835 in [1]).
        if (s>d) then
          dummy = (/ (w**(-j*(s-d)), j=0,m) /)
          nj = ceiling(N*dummy/sum(dummy))
        else
          stop "ERROR: s must be larger than d for Elliptically Contoured Distributions!"
        end if

      case ("M_SR")

        ! This is the most general case of weights with prescribed decay
        ! at infinity (See Section 3.2 in [1]).

        dummy = (/ (w**(-j*(s-d)), j=0,m) /)
        nj = ceiling(N*dummy/sum(dummy))

      case default

        stop "ERROR: unknown weight function type."

    end select

  end subroutine calculate_nj


  ! The indicator function for I_j (see also page 828 of [1]).
  ! The columns of y represent the different points to be checked.
  !
  function indicator_I(j, y) result (res)
    integer(kind=i4b), intent(in)               :: j
    real(kind=qp), dimension(:,:), intent(in)   :: y
    integer(kind=i4b), dimension(size(y,dim=2)) :: res

    res = 0 

    if (j==0) then
      where (maxval(abs(y), dim=1) < 0.5_qp*V)
        res = 1
      end where
    else
      where (maxval(abs(y), dim=1) > 0.5_qp*V*w**(j-1) &
        .and. maxval(abs(y), dim=1) < 0.5_qp*V*w**j)
        res = 1
      end where
    end if

  end function indicator_I


  ! Evaluate the sum formula (1) from page 2 in [1].
  !
  subroutine evaluate_s_mn(m, n, f, params_f, rho, params_rho, q, x, res)
    integer(kind=i4b), intent(in)                :: m
    integer(kind=i4b), dimension(0:), intent(in) :: n
    interface
      function f(x, params_f) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:,:), intent(in) :: x
        type(functionparams), intent(in)          :: params_f
        real(kind=qp), dimension(size(x,2))       :: y           !TODO: CHECK IF THIS SHOULD BE size(x,1) or size(x,2)
      end function f
    end interface
    type(functionparams), intent(in)             :: params_f
    interface
      function rho(x, params_rho) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:,:), intent(in) :: x
        type(functionparams), intent(in)          :: params_rho
        real(kind=qp), dimension(size(x,2))       :: y           !TODO: CHECK IF THIS SHOULD BE size(x,1) or size(x,2)
      end function rho
    end interface
    type(functionparams), intent(in)             :: params_rho
    real(kind=qp), dimension(0:), intent(in)     :: q
    real(kind=qp), dimension(:,:), intent(in)    :: x
    real(kind=qp), intent(out)                   :: res


    real(kind=qp), dimension(:,:), allocatable :: y
    real(kind=qp)                              :: temp
    integer(kind=i4b)                          :: j

    res = 0.0_qp

    do j=0,m ! loop over all m+1 areas/cubes

      ! Transform the first n_j points of the low-discrepancy pointset x to
      ! the current area/cube.
      allocate(y(size(x, dim=1), 1:n(j)))
      y = -0.5_qp*V*w**j + 0.5_qp*V*w**(j+1.0_qp)*x(:,1:n(j))

      ! Evaluate the part for I_j and add it to the result.
      temp = q(j)/n(j)*sum(f(y, params_f)*rho(y, params_rho)*indicator_I(j,y))
      res = res + temp
      deallocate(y)

      write(unit=*, fmt="(A, I0.2, A, ES22.15, A, I0.2, A, F22.15)") &
        "I_", j, ": ", temp, "       Q_", j, ": ", res

    end do

  end subroutine evaluate_s_mn


end module mod_mathewei
