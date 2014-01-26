! Module that implements testfunctions from Laurent Vanbeylen.

module f_vanbeylen

  use numeric_kinds
  use mod_function
  use mod_integration
  use mod_constants

  private

  public  :: f_vanbeylen01_1p
  public  :: init_y_t
  private :: h01
  private :: f01
 
  real(kind=qp), parameter, private :: a = 1,      &
                                       b = 1,      &
                                       alpha = 1,  &
                                       beta = 1,   &
                                       lambda = 1, &
                                       mu = 1

  real(kind=qp), dimension(:), allocatable, private :: y

  contains

    subroutine init_y_t(N)
      integer(kind=i4b), intent(in) :: N

      allocate(y(0:N-1))
      call random_number(y)

    end subroutine init_y_t


    function h01(x_t, x_tm1) result(h)
      real(kind=qp), intent(in) :: x_t, x_tm1
      real(kind=qp)             :: h

      h = a*x_t + b*x_tm1

    end function h01


    function f01(x_t, y_t) result(f)
      real(kind=qp), intent(in) :: x_t, y_t
      real(kind=qp)             :: f

      f = y_t + alpha*x_t + beta*(x_t**3)

    end function f01


    ! This is a function for integration over the real line!
    !
    function f_vanbeylen01_1p(x) result (f)
      real(kind=qp), dimension(0:), intent(in) :: x
      real(kind=qp)                            :: f

      integer(kind=i4b) :: t

      f = exp( - 0.5*((h01(x(0), 0.0_qp ))**2)/lambda &
               - 0.5*((f01(x(0), y(0)))**2)/mu )

      do t=1,size(x)-1
        f = f*exp( - 0.5*((h01(x(t), x(t-1)))**2)/lambda &
                   - 0.5*((f01(x(t), y(t)  ))**2)/mu )
      end do

    end function f_vanbeylen01_1p


end module f_vanbeylen
