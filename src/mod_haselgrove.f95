! Module that implements a method for numerical integration on periodic
! functions having certain requirements on the Fourier coefficients.
!
! References:
!
!   [1] "A Method for Numerical Integration", C. B. Haselgrove, Mathematics
!       of Computation, Vol. 15, 1961, pages 323-337.
!
! TODO:
!   * Check why we can't exactly reproduce Haselgrove's results.
!   * Implement the version that reduces rounding errors and avoids
!     the fact that S2 becomes too large.

module mod_haselgrove

  use numeric_kinds
  use mod_function
  use mod_constants
  use mod_utilities
  use mod_integration
  use f_testfunctions

  private

  public :: haselgrove
  public :: init_haselgrove

  real(kind=qp), dimension(:), allocatable, private :: has_irrationals


  contains

    subroutine init_haselgrove(s, alpha)

      integer(kind=i4b), intent(in)                     :: s
      real(kind=qp), dimension(:), intent(in), optional :: alpha

      if (allocated(has_irrationals)) then
        deallocate(has_irrationals)
      end if
      allocate(has_irrationals(s))

      select case (s)

        case (1)

          !has_irrationals = (/ 0.73258893_qp /)
          has_irrationals = (/ 0.83969144_qp /)

        case (2)

          !has_irrationals = (/ 0.62055505_qp, 0.22610245_qp /)
          has_irrationals = (/ 0.59734470_qp, 0.92828094_qp /)

        case (3)

          !has_irrationals = (/ 0.96498949_qp, 0.81091316_qp, 0.46960090_qp /)
          has_irrationals = (/ 0.74235492_qp, 0.57387033_qp, 0.32279917_qp /)

        case (4)

          !has_irrationals = (/ 0.62366851_qp, 0.04150108_qp, 0.48574769_qp, &
          !                      0.27210703_qp /)
          has_irrationals = (/ 0.17665781_qp, 0.71327190_qp, 0.98875216_qp, &
                               0.60299793_qp /)

        case (5)

          !has_irrationals = (/ 0.95734608_qp, 0.86730270_qp, 0.09724025_qp, &
          !                      0.31301950_qp, 0.48476582_qp /)
          has_irrationals = (/ 0.44810200_qp, 0.53589831_qp, 0.56039410_qp, &
                               0.83630131_qp, 0.22148205_qp /)

        case (6)

          !has_irrationals = (/ 0.43657951_qp, 0.59185199_qp, 0.05024400_qp, &
          !                      0.84373919_qp, 0.3810400_qp, 0.75808683_qp /)
          has_irrationals = (/ 0.10613747_qp, 0.40278232_qp, 0.88772556_qp, &
                               0.43554826_qp, 0.17219381_qp, 0.63794472_qp /)

        case (7)

          !has_irrationals = (/ 0.80638723_qp, 0.22584927_qp, 0.72510075_qp, &
          !                      0.51310685_qp, 0.11080509_qp, 0.60161858_qp, &
          !                      0.92715171_qp /)
          has_irrationals = (/ 0.58505729_qp, 0.50196855_qp, 0.77797734_qp, &
                               0.60504620_qp, 0.62193588_qp, 0.84244165_qp, &
                               0.64543976_qp /)

        case (8)

          !has_irrationals = (/ 0.73750248_qp, 0.08314415_qp, 0.84753682_qp, &
          !                      0.88989711_qp, 0.80254484_qp, 0.27951501_qp, &
          !                      0.67340402_qp, 0.53040927_qp /)
          has_irrationals = (/ 0.23975940_qp, 0.01544979_qp, 0.57794809_qp, &
                               0.81182909_qp, 0.78068912_qp, 0.62319488_qp, &
                               0.70710061_qp, 0.60389317_qp /)

        case default

          has_irrationals = alpha

      end select

    end subroutine init_haselgrove


    subroutine haselgrove(f, params, N, r, res)

      interface
        function f(x, fparams) result (y)
          use numeric_kinds
          use mod_function
          real(kind=qp), dimension(:), intent(in) :: x
          type(functionparams), intent(in)        :: fparams
          real(kind=qp)                           :: y
        end function f
      end interface
      type(functionparams), intent(in) :: params
      integer(kind=i4b), intent(in)    :: N
      integer(kind=i4b), intent(in)    :: r
      real(kind=qp), intent(out)       :: res

      integer(kind=i4b)                :: i, result_index
      real(kind=qp)                    :: S1, S2, S2_old, lastrow
      real(kind=qp), dimension(N/1000) :: abs_err
      real(kind=qp), dimension(N/1000) :: cumres
      real(kind=qp)                    :: exact
      type(integrationbounds)          :: bounds


      call allocate_integration_bounds(bounds, size(has_irrationals))
      call set_unit_hypercube(bounds)
      exact =  f_haselgrove61_exact(params, bounds)
      !write(unit=*, fmt="(A, F10.8)") "Exact value for integral is: ", exact
      print *, "Exact value for integral is: ", exact

      select case (r)

        case (1)

          S1 = f(0.0_qp*2.0_qp*pi*has_irrationals, params)
          result_index = 0
          do i=1,N

            ! This is as Haselgrove described it, but it does not work?
            ! TODO:

            ! This works, but is it as Haselgrove meant it?
            S1 = S1 + 2.0_qp*f(frac_part(i*has_irrationals), params)

            if (modulo(i, 1000) == 0) then
              result_index = result_index + 1
              cumres(result_index) = S1/(2*i+1)
              abs_err(result_index) = cumres(result_index)-exact
            end if

          end do
          res = S1/(2*N+1)


        ! TODO: when N gets large (say 100000) then probably S2 gets too large and
        !       rounding errors break the computations...
        case (2)

          S2 = f(0.0_qp*2.0_qp*pi*has_irrationals, params)
          lastrow = 0.0_qp
          result_index=0
          do i=1,N

            S2_old = S2

            ! This is as Haselgrove described it, but it does not work?
            !lastrow = lastrow + 2*f(2*abs(frac_part(0.5_qp*i*has_irrationals)), params)

            ! This works, but is it as Haselgrove meant it?
            lastrow = lastrow + 2.0_qp*f(frac_part(i*has_irrationals), params)
            S2 = S2_old + lastrow
            print *, "S2 = ", S2

            if (modulo(i, 1000) == 0) then
              result_index = result_index + 1
              cumres(result_index) = S2/(i+1)**2
              abs_err(result_index) = cumres(result_index)-exact
            end if

          end do
          res = S2/(N+1)**2

        case (3)

          print *, "WARNING: Haselgrove's method for r=3 not implemented yet!"

        case (4)

          print *, "WARNING: Haselgrove's method for r=4 not implemented yet!"

      end select

      open(unit=20, file="haselgrove.dat", &
           status="replace", access="sequential", action="write")
      do i=1,N/1000
        write(unit=20, fmt="(I5, F12.8, ES25.15e2)") i*1000, cumres(i), abs_err(i)
      end do
      close(unit=20)

    end subroutine haselgrove

end module mod_haselgrove
