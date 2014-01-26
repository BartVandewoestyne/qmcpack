! Module that implements Niederreiter's improvement on Haselgrove's method.
!
! References: 
!
!   [1] Niederreiter, H., `Application of Diophantine approximations to
!       numerical integration' in `Diophantine Approximation and its
!       Applications', Osgood, Charles F., Academic Press 1973, pages 120-199,
!       ISBN 0-12-528650-3.

module mod_haselgrove_niederreiter

  use numeric_kinds
  use mod_function
  use mod_utilities
  use mod_periodize

  private

  public :: init_nied_q2_weyl
  public :: next_nied_q2_weyl
  public :: init_nied_q2_general
  public :: next_nied_q2_general
  public :: init_nied_q3_weyl
  public :: next_nied_q3_weyl
  public :: init_nied_q3_general
  public :: next_nied_q3_general
  public :: init_nied_weyl
  public :: next_nied_weyl

  ! The irrationals used in the method.
  real(kind=qp), dimension(:), allocatable, private       :: irrationals

  ! The parameter q which is also the asymptotic order of convergence.
  integer(kind=i4b), private                              :: q

  ! The name of the periodizer that will be used.
  character(len=MAX_LENGTH_PERIODIZER_NAME), private      :: periodizer

  ! The N from the method.
  integer(kind=i8b), private                              :: counter

  ! Sums and variables used in the 2nd order method.
  real(kind=qp), private                                  :: S1, S2, S
  real(kind=qp), private                                  :: f_new1, &
                                                             f_new2, &
                                                             f_new3, &
                                                             f_drop

  ! Differences and variables used for the 3rd order method.
  real(kind=qp), private                                  :: D1, D2, D3
  real(kind=qp), private                                  :: f_1n3, &
                                                             f_2n3, &
                                                             f_2n4, &
                                                             f_2n5, &
                                                             f_2n6, &
                                                             f_3n6, &
                                                             f_3n7, &
                                                             f_3n8, &
                                                             f_3n9

  ! Differences and coefficients a_{N,n}^q for the q'th order method.
  real(kind=qp), dimension(:), allocatable, private       :: d
  integer(kind=i4b), dimension(:,:), allocatable, private :: a, a_prev


  contains


    ! Initialize Niederreiter's improved method for q=2 which uses the
    ! Weyl sequence.
    !
    subroutine init_nied_q2_weyl(irrationals_init, periodizer_init, func, params)
      real(kind=qp), dimension(:), intent(in) :: irrationals_init
      character(len=*), intent(in)            :: periodizer_init
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

      if (allocated(irrationals)) then
        deallocate(irrationals)
      end if
      allocate(irrationals(size(irrationals_init)))
      irrationals = irrationals_init

      periodizer = periodizer_init

      counter = 1

      call periform(func, params, frac_part(0*irrationals), periodizer, S2)
      S1 = 0
      f_drop = S2
      f_new1 = 0
      f_new2 = 0

    end subroutine init_nied_q2_weyl


    ! Initialize Niederreiter's improved algorithm for q=2 and
    ! for general sequences with points that are independent of N.
    !
    subroutine init_nied_q2_general(periodizer_init, x, func, params)
      character(len=*), intent(in)              :: periodizer_init
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

      periodizer = periodizer_init
      counter = 1

      call periform(func, params, x(1,:), periodizer, S2)
      S1 = 0
      f_drop = S2
      f_new1 = 0
      f_new2 = 0

    end subroutine init_nied_q2_general


    ! Initialize Niederreiter's improved algorithm for q=3 and
    ! with the Weyl sequence.  The first two steps are unrolled.
    !
    subroutine init_nied_q3_weyl(irrationals_init, periodizer_init, func, params)
      real(kind=qp), dimension(:), intent(in) :: irrationals_init
      character(len=*), intent(in)            :: periodizer_init
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

      real(kind=qp) :: f0, f1, f2, f3


      if (allocated(irrationals)) then
        deallocate(irrationals)
      end if
      allocate(irrationals(size(irrationals_init)))
      irrationals = irrationals_init

      periodizer = periodizer_init

      ! The first 2 steps are unrolled.
      counter = 2

      call periform(func, params, frac_part(0*irrationals), periodizer, f0)
      call periform(func, params, frac_part(1*irrationals), periodizer, f1)
      call periform(func, params, frac_part(2*irrationals), periodizer, f2)
      call periform(func, params, frac_part(3*irrationals), periodizer, f3)

      D1 = 3*f1 + 3*f2 + f3
      D2 = D1 - f0
      S = D1 + f0

      print *, "D2 = ", D2
      print *, "D1 = ", D1
      print *, "D0 = ", S

      f_2n4 = f0
      f_2n3 = f1
      f_3n6 = f0
      f_new1 = f1
      f_new2 = f2
      f_new3 = f3

    end subroutine init_nied_q3_weyl


    ! Initialize Niederreiter's improved algorithm for q=3 and
    ! for general sequences with points that are independent of N.
    ! The first two steps are unrolled.
    !
    subroutine init_nied_q3_general(periodizer_init, x, func, params)
      real(kind=qp), dimension(:,:), intent(in) :: x
      character(len=*), intent(in)              :: periodizer_init
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

      real(kind=qp) :: f0, f1, f2, f3

      periodizer = periodizer_init

      ! The first 2 steps are unrolled.
      counter = 2

      call periform(func, params, x(1,:), periodizer, f0)
      call periform(func, params, x(2,:), periodizer, f1)
      call periform(func, params, x(3,:), periodizer, f2)
      call periform(func, params, x(4,:), periodizer, f3)

      D1 = 3*f1 + 3*f2 + f3
      D2 = D1 - f0
      S = D1 + f0

      f_2n4 = f0
      f_2n3 = f1
      f_3n6 = f0
      f_new1 = f1
      f_new2 = f2
      f_new3 = f3

    end subroutine init_nied_q3_general


    ! Perform an extra step for Niederreiter's improved method with q=2 and
    ! using the Weyl sequence.
    !
    subroutine next_nied_q2_weyl(func, params, integral_value)
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

      counter = counter + 1

      call periform(func, params, &
             frac_part( 2*(counter-1)   *irrationals), periodizer, f_new1)
      call periform(func, params, &
             frac_part((2*(counter-1)-1)*irrationals), periodizer, f_new2)
      call periform(func, params, &
             frac_part(   (counter-2)   *irrationals), periodizer, f_drop)
      S1 = S2 - f_drop + f_new2
      S2 = S1 - f_new1
      S = S + S1 + S2

      integral_value = (S/counter)/counter

    end subroutine next_nied_q2_weyl


    ! Do an extra step in Niederreiter's improved algorithm for q=2 and
    ! with points that are not dependent on N.
    !
    subroutine next_nied_q2_general(x, func, params, integral_value)
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
      real(kind=qp), intent(out)                :: integral_value

      counter = counter + 1
      call periform(func, params, x(2*(counter-1)  +1,:), periodizer, f_new1)
      call periform(func, params, x(2*(counter-1)-1+1,:), periodizer, f_new2)
      call periform(func, params, x(   counter-2   +1,:), periodizer, f_drop)
      S1 = S2 - f_drop + f_new2
      S2 = S1 + f_new1
      S = S + S1 + S2

      integral_value = (S/counter)/counter

    end subroutine next_nied_q2_general


    ! Perform an extra step for Niederreiter's improved method with q=3 and
    ! using the Weyl sequence.
    !
    subroutine next_nied_q3_weyl(func, params, integral_value)
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

      counter = counter + 1

      ! calculate new function values
      call periform(func, params, &
                frac_part(  (counter-3)*irrationals), periodizer, f_1n3)
      f_2n6 = f_2n4
      f_2n5 = f_2n3
      call periform(func, params, &
                frac_part((2*counter-4)*irrationals), periodizer, f_2n4)
      call periform(func, params, &
                frac_part((2*counter-3)*irrationals), periodizer, f_2n3)
      f_3n9 = f_3n6
      f_3n8 = f_new1
      f_3n7 = f_new2
      f_3n6 = f_new3
      call periform(func, params, &
                frac_part((3*counter-5)*irrationals), periodizer, f_new1)
      call periform(func, params, &
                frac_part((3*counter-4)*irrationals), periodizer, f_new2)
      call periform(func, params, &
                frac_part((3*counter-3)*irrationals), periodizer, f_new3)

      ! Update the differences.
      D3 = 3*f_1n3 - 3*(f_2n6 + 3*(f_2n5 + f_2n4) + f_2n3) + f_3n9 + f_new3 &
           + 3*(f_3n8 + f_new2) + 6*(f_3n7 + f_new1) + 7*f_3n6
      D2 = D2 + D3
      D1 = D1 + D2
      S = S + D1

      !print *, "D3 = ", D3
      !print *, "D2 = ", D2
      !print *, "D1 = ", D1
      !print *, "S = ", S

      ! returned integral value
      integral_value = ((S/counter)/counter)/counter

      !print *, "I = ", integral_value

    end subroutine next_nied_q3_weyl


    ! Do an extra step in Niederreiter's improved algorithm for q=2 and
    ! with points that are not dependent of N.
    !
    subroutine next_nied_q3_general(x, func, params, integral_value)
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
      real(kind=qp), intent(out)                :: integral_value


      counter = counter + 1

      ! calculate new function values
      call periform(func, params, x(counter-3+1,:), periodizer, f_1n3)
      f_2n6 = f_2n4
      f_2n5 = f_2n3
      call periform(func, params, x(2*counter-4+1,:), periodizer, f_2n4)
      call periform(func, params, x(2*counter-3+1,:), periodizer, f_2n3)
      f_3n9 = f_3n6
      f_3n8 = f_new1
      f_3n7 = f_new2
      f_3n6 = f_new3
      call periform(func, params, x(3*counter-5+1,:), periodizer, f_new1)
      call periform(func, params, x(3*counter-4+1,:), periodizer, f_new2)
      call periform(func, params, x(3*counter-3+1,:), periodizer, f_new3)

      ! update the differences
      D3 = 3*f_1n3 - 3*(f_2n6 + 3*(f_2n5 + f_2n4) + f_2n3) + f_3n9 + f_new3 &
           + 3*(f_3n8 + f_new2) + 6*(f_3n7 + f_new1) + 7*f_3n6
      D2 = D2 + D3
      D1 = D1 + D2
      S = S + D1

      ! returned integral value
      integral_value = ((S/counter)/counter)/counter


    end subroutine next_nied_q3_general


    ! Initialize Niederreiter's improved method with the specified
    ! order of convergence q and the specified irrationals.
    !
    subroutine init_nied_weyl(q_init, irrationals_init, func, params, periodizer_init)
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


      integer(kind=i4b) :: i, j, k, curr_q
      real(kind=qp)     :: f

      q = q_init

      if (allocated(irrationals)) then
        deallocate(irrationals)
      end if
      allocate(irrationals(size(irrationals_init)))
      irrationals = irrationals_init

      periodizer = periodizer_init


      ! Allocate storage for the calculation of the coefficients a_{N,n}^(q).
      ! We will need these coefficients for the initialization and also in each
      ! step to calculate the q't order backward difference D_{N,q}^{(q)}.
      allocate(a(0:q,0:q*(q-1)), a_prev(0:q,0:q*(q-1)))

      ! Allocate space for the zero'th up to the q'th order differences.
      allocate(d(0:q))

      ! The first q-1 steps are unrolled.
      counter = q-1

      ! Initialize for N=0..q-1 (or N=0..q ??? CHECK) and for q=0.  Note that for q=0 all the
      ! coefficients are 1, except when N=0, then the coefficient is 0.
      a_prev = 0
      a_prev(1:,0) = 1

      ! Now step through the q-values and use the recurrence in q to
      ! update the table.
      do curr_q=1,q

        a = 0

        ! Use the recurrence in q to update a from a_prev (q-1) to a (q).
        ! Calculate one extra row because we will need it in the calculation
        ! of D_q later on.  Note also that we do not have to update the row a(0,:)
        ! because we know that for N=0 we always have a value of 0.
        do i = 1,q
          do j = 0,curr_q*(i-1)
            do k = max(0, j-(curr_q-1)*(i-1)), min(i-1, j)
              a(i,j) = a(i,j) + a_prev(i, j-k)
            end do
          end do
        end do

        a_prev = a

      end do

      ! Find the zero'th order differences D_{i,0}^{(q)} for i=0,...,q-1.  Note
      ! that we do not have to calculate D_{0,0}^{(q)} because we know it is 0
      ! because N=0 here.
      d = 0
      do i=1,q-1
        do j=0,q*(i-1)
          call periform(func, params, frac_part(j*irrationals), periodizer, f)
          d(q-1-i) = d(q-1-i) + a(i,j)*f
        end do
      end do

      ! Now update the table of the differences so that d(i) = D_{q-1,i}^{(q)} = i'th order
      ! backward difference.  Note that the inner loop here runs from q-1 to i in order
      ! to avoid overwriting certain values in the wrong order.
      do i=1,q-1
        do j=q-1,i,-1
          d(j) = d(j-1)-d(j)
        end do
      end do

    end subroutine init_nied_weyl


    ! Perform an extra step for Niederreiter's improved method with
    ! general q and using the Weyl sequence.
    !
    subroutine next_nied_weyl(func, params, integral_value)
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

      real(kind=qp)                     :: temp, f
      integer(kind=i4b), dimension(0:q) :: binom
      integer(kind=i4b)                 :: j, k


      counter = counter + 1

      call binomial_coefficients(q, binom)

      ! Calculate q'th order difference D_{counter,q}^{(q)}.
      d(q) = 0
      do j=1,q
        temp = 0
        do k = 0,q*(j-1)
          call periform(func, params,                             &
                        frac_part((k+j*(counter-q))*irrationals), &
                        periodizer, f)
          temp = temp + a(j,k)*f
        end do
        d(q) = d(q) + binom(j)*(-1.0)**(q-j)*temp
      end do

      ! Update the other differences.
      do j=q-1,0,-1
        d(j) = d(j) + d(j+1)
      end do

      ! Calculate the final integral value.
      integral_value = d(0)
      do j=1,q
        integral_value = integral_value/counter
      end do

    end subroutine next_nied_weyl

end module mod_haselgrove_niederreiter
