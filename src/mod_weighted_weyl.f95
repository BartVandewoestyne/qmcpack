! This is an experimental module, used for experimenting with
! different weight functions in combination with the Weyl sequence.

module mod_weighted_weyl

  use numeric_kinds
  use mod_function
  use mod_constants
  use mod_special_functions
  use mod_utilities
  use mod_periodize
  use mod_weightfunction

  private

  public :: weighted_weyl

contains


  subroutine weighted_weyl(weightfunction_name, N, irrationals, func, params, periodizer, integral_value)
    character(len=*), intent(in)            :: weightfunction_name
    integer(kind=i4b), intent(in)           :: N
    real(kind=qp), intent(in), dimension(:) :: irrationals
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
    real(kind=qp), intent(out)              :: integral_value

    integer(kind=i4b) :: k
    real(kind=qp)     :: f
    real(kind=qp)     :: my_sum


    my_sum = 0.0_qp

    select case (weightfunction_name)

      case("SagSzekeres")

        do k=1,N-1
          call periform(func, params, frac_part(k*irrationals), periodizer, f)
          my_sum = my_sum + sag_szekeres(real(k, kind=qp)/N)*f
        end do
        integral_value = my_sum/(N-1)

      case("IMT")

        do k=1,N-1
          call periform(func, params, frac_part(k*irrationals), periodizer, f)
          my_sum = my_sum + imt(real(k, kind=qp)/N)*f
        end do
        integral_value = my_sum/(N-1)

      case("SugiharaMurotaGaussian")

        do k=1,N
          call periform(func, params, frac_part(k*irrationals), periodizer, f)
          my_sum = my_sum + sugihara_murota_gaussian(real(k, kind=qp)/N, 1.0_qp, k, size(irrationals))*f
        end do
        integral_value = my_sum/N

      case default
        stop "ERROR: unknown weightfunction name!"

    end select

  end subroutine weighted_weyl


end module mod_weighted_weyl
