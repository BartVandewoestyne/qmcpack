! Module that implements Sugihara and Murota's modification of Niederreiter's
! generalization of Haselgrove's method for numerical integration.  Sugihara
! and Murota show in [1] that O(N^(-k)) convergence can be reached, but in
! [2], Kaino improves this and shows that O(N^(-(k+1))) convergence can be
! reached.
!
! References:
!
!   [1] Sugihara, Masaaki and Murota, Kazuo, `A Note on Haselgrove's Method
!       for Numerical Integration', Mathematics of Computation, volume 39,
!       number 160, October 1982, pages 549-554.
!
!   [2] `Another Note on Haselgrove's Method for Numerical Integration',
!       Kaino K., Journal of the Korean Physical Society, Vol. 40, No. 6,
!       June 2002, pp. 1010-1014.
!
module mod_sugihara_murota

  use numeric_kinds
  use mod_utilities
  use mod_function
  use mod_periodize
  use mod_special_functions

  private

  public :: init_sugihara_murota
  public :: next_sugihara_murota
  public :: init_sugihara_murota_k2_weyl
  public :: next_sugihara_murota_k2_weyl
  public :: init_sugihara_murota_k2_general
  public :: next_sugihara_murota_k2_general
  private :: calculate_Aq

  ! The irrationals used in the method.
  real(kind=qp), dimension(:), allocatable, private     :: irrationals

  ! The parameter q which represents the asymptotic order of convergence.
  integer(kind=i4b), private                            :: q

  ! The name of the periodizer that will be used.
  character(len=MAX_LENGTH_PERIODIZER_NAME), private    :: periodizer

  ! The N from the method (representing the number of points used).
  integer(kind=i4b), private                            :: counter

  ! Some other internal variables used in this algorithm.
  real(kind=qp), dimension(:), allocatable, private     :: T
  real(kind=qp), private                                :: a, b, c, Aq
  integer(kind=i4b), dimension(:), allocatable, private :: sm_B

contains


  ! Initialize Sugihara and Murota's q'th order method for the given function.
  !
  subroutine init_sugihara_murota(q_init, irrationals_init, func, params, periodizer_init)
    integer(kind=i4b), intent(in)           :: q_init
    real(kind=qp), dimension(:), intent(in) :: irrationals_init
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
    character(len=*), intent(in)            :: periodizer_init

    real(kind=qp) :: f


    counter = 1
    q = q_init

    if (allocated(irrationals)) then
      deallocate(irrationals)
    end if
    allocate(irrationals(size(irrationals_init)))
    irrationals = irrationals_init

    periodizer = periodizer_init


    ! The binomial coefficients with the appropriate power of (-1).
    if (allocated(sm_B)) then
      deallocate(sm_B)
    end if
    allocate(sm_B(0:q))
    call binomial_coefficients(q, sm_B)
    sm_B(1::2) = -sm_B(1::2)

    ! The sums
    if (allocated(T)) then
      deallocate(T)
    end if
    allocate(T(0:q))
    call periform(func, params, frac_part(irrationals), periodizer, f)
    T = sm_B*f

    ! Calculate the constant that is in front of the sum
    Aq = calculate_Aq(q)

  end subroutine init_sugihara_murota


  ! Initialize Sugihara and Murota's second order method for the given function.
  !
  subroutine init_sugihara_murota_k2_weyl(irrationals_init, func, params, periodizer)
    real(kind=qp), dimension(:), intent(in) :: irrationals_init
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

    if (allocated(irrationals)) then
      deallocate(irrationals)
    end if
    allocate(irrationals(size(irrationals_init)))

    irrationals = irrationals_init

    counter = 1

    call periform(func, params, frac_part(irrationals), periodizer, a)
    b = -2*a
    c = a

  end subroutine init_sugihara_murota_k2_weyl


  ! Initialize Sugihara and Murota's second order method for the given function
  ! and for general pointsets.
  !
  subroutine init_sugihara_murota_k2_general(x, func, params, periodizer)
    real(kind=qp), dimension(:,:), intent(in) :: x
    interface
      function func(x, fparams) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function func
    end interface
    type(functionparams), intent(in)          :: params
    character(len=*), intent(in)              :: periodizer

    counter = 1

    call periform(func, params, x(1,:), periodizer, a)
    b = -2*a
    c = a

  end subroutine init_sugihara_murota_k2_general


  ! Return the next estimate when applying the loop with
  ! Sugihara and Murota's weight-function with general q.
  !
  subroutine next_sugihara_murota(func, params, integral_value)
    interface
      function func(x, fparams) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function func
    end interface
    type(functionparams), intent(in) :: params
    real(kind=qp), intent(out)       :: integral_value

    real(kind=qp)     :: f_new
    integer(kind=i4b) :: i


    counter = counter + 1

    call periform(func, params, &
                  frac_part(counter*irrationals), periodizer, f_new)

    ! update the T-vector
    T = T*( ((counter-1.0_qp)/counter)**(q+(/ (i, i=0,q) /)) ) + sm_B*f_new

    integral_value = (sum(T)*Aq/counter)

  end subroutine next_sugihara_murota


  ! Return the next estimate when applying the loop with
  ! Sugihara and Murota's weight-function with k=2.
  !
  subroutine next_sugihara_murota_k2_weyl(func, params, periodizer, integral_value)
    interface
      function func(x, fparams) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function func
    end interface
    type(functionparams), intent(in) :: params
    character(len=*), intent(in)     :: periodizer
    real(kind=qp), intent(out)       :: integral_value

    real(kind=qp) :: f_new

    counter = counter + 1

    call periform(func, params, &
                  frac_part(counter*irrationals), periodizer, f_new)

    a = ((counter-1.0_qp)/counter)**2*a + f_new
    b = ((counter-1.0_qp)/counter)**3*b - 2*f_new
    c = ((counter-1.0_qp)/counter)**4*c + f_new

    integral_value = ((a+b+c)/counter)*30

  end subroutine next_sugihara_murota_k2_weyl


  ! Return the next estimate when applying the loop with
  ! Sugihara and Murota's weight-function with k=2 and
  ! for general sequences.
  !
  subroutine next_sugihara_murota_k2_general(x, func, params, periodizer, integral_value)
    real(kind=qp), dimension(:,:), intent(in) :: x
    interface
      function func(x, fparams) result (y)
        use numeric_kinds
        use mod_function
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function func
    end interface
    type(functionparams), intent(in)          :: params
    character(len=*), intent(in)              :: periodizer
    real(kind=qp), intent(out)                :: integral_value

    real(kind=qp) :: f_new

    counter = counter + 1

    call periform(func, params, x(counter,:), periodizer, f_new)

    a = ((counter-1.0_qp)/counter)**2*a + f_new
    b = ((counter-1.0_qp)/counter)**3*b - 2*f_new
    c = ((counter-1.0_qp)/counter)**4*c + f_new

    integral_value = ((a+b+c)/counter)*30

  end subroutine next_sugihara_murota_k2_general


  ! Calculate the factor A_q that is in front of the
  ! weight function.
  !
  function calculate_Aq(q) result (Aq)
    integer(kind=i4b), intent(in) :: q
    real(kind=qp)                 :: Aq

    integer(kind=i4b) :: i

    Aq = product(  (/ ((2*q+1.0_qp-i)/(q-i), i=0,q-1) /) ) * (q+1)

  end function calculate_Aq

end module mod_sugihara_murota
