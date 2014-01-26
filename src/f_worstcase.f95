! Module containing the worst possible functions in the class E_alhpa(1), when
! integrating with lattices (See [1], page 69-73.)
!
! References:
!
!  [1] Sloan, I. H. and Joe, S., `Lattice Methods for Multiple Integration',
!      Oxford Science Publications, 1994, ISBN 0-19-853472-8.

module f_worstcase

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_constants

  private

  public :: f_worstcase2_1p
  public :: f_worstcase2_exact
  public :: f_worstcase4_1p
  public :: f_worstcase4_exact
  public :: f_worstcase6_1p
  public :: f_worstcase6_exact
  public :: f_worstcase_infinity_1p
  public :: f_worstcase_infinity_exact
  public :: f_saltykov_1p
  public :: f_saltykov_exact

contains


    !                   2   2
    !   f(x) := 1 + 2 Pi  (x  - x + 1/6)
    !
    ! See [1], page 73.
    !
    function f_worstcase2_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( 1 + 2*pi*pi*(x*x - x + 1.0_qp/6.0_qp) )

    end function f_worstcase2_1p

    function f_worstcase2_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( 2.0_qp/3.0_qp*pi*pi*(ub*ub*ub-lb*lb*lb)   &
                    - pi*pi*(ub*ub-lb*lb)                    &
                    + ub-lb + 1.0_qp/3.0_qp*pi*pi*(ub-lb) )

    end function f_worstcase2_exact


    !                 4          2        2
    !               Pi  (1 - 30 x  (1 - x) )
    !   f(x) := 1 + ------------------------
    !                          45
    !
    ! See [1], page 73.
    !
    function f_worstcase4_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( 1 + pi**4*(1 - 30*x*x*(1-x)*(1-x))/45 )

    end function f_worstcase4_1p

    function f_worstcase4_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i =  product( - 2*pi**4*(ub**5 - lb**5)/15    &
                    +   pi**4*(ub**4 - lb**4)/3     &
                    - 2*pi**4*(ub**3 - lb**3)/9     &
                    +          ub    - lb           &
                    +   pi**4*(ub - lb)/45        )

    end function f_worstcase4_exact


    !                   6          2        4        5       6
    !               2 pi  (1 - 21 x  + 105 x  - 126 x  + 42 x )
    !   f(x) := 1 + -------------------------------------------
    !                                   945
    !
    ! See [1], page 73.
    !
    function f_worstcase6_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( 1 + (2*pi**6/945)*(1 - 21*x**2 + 105*x**4    &
                                        - 126*x**5 + 42*x**6) )

    end function f_worstcase6_1p

    function f_worstcase6_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product(   (4*pi**6/315)*(ub**7 - lb**7)   &
                   - (2*pi**6/45 )*(ub**6 - lb**6)   &
                   + (2*pi**6/45 )*(ub**5 - lb**5)   &
                   - (2*pi**6/135)*(ub**3 - lb**3)   &
                   +                ub    - lb       &
                   + (2*pi**6/945)*(ub - lb)       )

    end function f_worstcase6_exact


    !                 s
    !              --------'
    !             '  |  |
    !     f(x) :=    |  |    (1 + 2 cos(2 Pi x[i]))
    !                |  |
    !                |  |
    !               i = 1
    !
    ! See [1], page 73.
    !
    function f_worstcase_infinity_1p(x, params) result (y)
      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( 1 + 2*cos(2*pi*x) )

    end function f_worstcase_infinity_1p

    function f_worstcase_infinity_exact(params, bounds) result (i)
      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( (ub-lb) + (sin(2*pi*ub)-sin(2*pi*lb))/pi )

    end function f_worstcase_infinity_exact


    !                  s
    !               --------'
    !              '  |  |                   2
    !      f(x) :=    |  |    (3 (2 x[i] - 1) )
    !                 |  |
    !                 |  |
    !                i = 1
    !
    ! See [1], page 75.
    !
    function f_saltykov_1p(x, params) result (y)

      real(kind=qp), dimension(:), intent(in) :: x
      type(functionparams), intent(in)        :: params
      real(kind=qp)                           :: y

      y = product( 3*(2*x-1)*(2*x-1) )

    end function f_saltykov_1p

    function f_saltykov_exact(params, bounds) result (i)

      type(functionparams), intent(in)    :: params
      type(integrationbounds), intent(in) :: bounds
      real(kind=qp)                       :: i

      real(kind=qp), dimension(size(bounds%lb)) :: lb, ub

      lb = bounds%lb
      ub = bounds%ub

      i = product( 4*(ub*ub*ub - lb*lb*lb) - 6*(ub*ub - lb*lb) + 3*(ub - lb) )

    end function f_saltykov_exact

end module f_worstcase
