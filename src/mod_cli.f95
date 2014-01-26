! Module implementing some general purpose Command Line Interface (CLI)
! routines.
!
module mod_cli

  use numeric_kinds
  use mod_sobol
  use mod_halton
  use mod_weyl
  use mod_extensible_lattice
  use mod_primes

  private

  public :: ask_dimension 
  public :: ask_nb_points 
  public :: ask_testfunction_name
  public :: ask_periodizer_name
  public :: ask_weightfunction_name
  public :: ask_pointset_name
  public :: create_pointset

contains

  ! Ask for the dimension.
  !
  subroutine ask_dimension(s)

    integer(kind=i4b), intent(out) :: s

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the dimension of the problem: "
    read *, s

  end subroutine ask_dimension


  ! Ask for the number of points.
  !
  subroutine ask_nb_points(n)

    integer(kind=i4b), intent(out) :: n

    write(unit=*, fmt="(A)", advance="no") &
      "Enter the number of points: "
    read *, n

  end subroutine ask_nb_points


  ! Ask for the name of a testfunction and print the selected function to
  ! screen.
  !
  subroutine ask_testfunction_name(func_name)

    character(len=*), intent(out) :: func_name

    write(unit=*, fmt="(A)") "Available testfunctions:"
    write(unit=*, fmt="(A)") "  f_simple_product (product(x))"
    write(unit=*, fmt="(A)") "  f_warnock01 (product((2*exp(2*x)-1)/(exp(2)-1)))"
    write(unit=*, fmt="(A)") "  f_warnock02 (product( 3*x**5 + (sin(6*pi*x))**2 ))"
    write(unit=*, fmt="(A)") "  f_warnock02_bernoulli3"
    write(unit=*, fmt="(A)") "  f_warnock02_bernoulli5"
    write(unit=*, fmt="(A)") "  f_warnock03 (product( 30*x**29 ))"
    write(unit=*, fmt="(A)") "  f_fox86"
    write(unit=*, fmt="(A)") "  f_nuyens01 (product(1+sin(x)))"
    write(unit=*, fmt="(A)") "  f_ronald01 (x(1)**2)"
    write(unit=*, fmt="(A)") "  f_ronald02 (product(x**2))"
    write(unit=*, fmt="(A)") "  f_ronald_polynom6 (product(x**6))"
    write(unit=*, fmt="(A)") "  f_ronald_polynom6_bernoulli5 (product(x**6))"
    write(unit=*, fmt="(A)") "  f_ronald03 (product(sin(2*pi*x)))"
    write(unit=*, fmt="(A)") "  f_ronald04 (product(cos(2*pi*x)))"
    write(unit=*, fmt="(A)") "  f_ronald05 (product(1+sin(2*pi*x)))"
    write(unit=*, fmt="(A)") "  f_worstcase2"
    write(unit=*, fmt="(A)") "  f_worstcase4"
    write(unit=*, fmt="(A)") "  f_worstcase6"
    write(unit=*, fmt="(A)") "  f_sine_power"
    write(unit=*, fmt="(A)") "  f_sine_power_bernoulli5"
    write(unit=*, fmt="(A)", advance="no") "Enter function name: "
    read *, func_name
    write(unit=*, fmt="(A)") "You are running your tests on the following testfunction:"
    select case (func_name)

      case ("f_simple_product")

         write(unit=*, fmt="(A)") "f(x) = product(x)"

      case ("f_ronald05")

         write(unit=*, fmt="(A)") "f(x) = product(1+sin(2*pi*x))"

      case default

         write(unit=*, fmt="(A)") "f(x) = TODO"

    end select

  end subroutine ask_testfunction_name


  ! Ask for the name of a periodizer.
  !
  subroutine ask_periodizer_name(periodizer_name)
  
    character(len=*), intent(out) :: periodizer_name

    write(unit=*, fmt="(A)") "Available periodizing methods:"
    write(unit=*, fmt="(A)") "  none"
    write(unit=*, fmt="(A)") "  Laurie_m2 (1996)"
    write(unit=*, fmt="(A)") "  Sidi_m2 (1993)"
    write(unit=*, fmt="(A)") "  Mori (1978)"
    write(unit=*, fmt="(A)") "  DE (1974)"
    write(unit=*, fmt="(A)") "  IMT (1970)"
    write(unit=*, fmt="(A)") "  TANH (1964)"
    write(unit=*, fmt="(A)") "  Korobov_m1 (1963)"
    write(unit=*, fmt="(A)") "  Korobov_m2 (1963)"
    write(unit=*, fmt="(A)") "  Korobov_m3 (1963)"
    write(unit=*, fmt="(A)") "  Korobov_m4 (1963)"
    write(unit=*, fmt="(A)") "  Korobov_m5 (1963)"
    write(unit=*, fmt="(A)") "  Korobov_m6 (1963)"
    write(unit=*, fmt="(A)") "  Korobov_m7 (1963)"
    write(unit=*, fmt="(A)") "  Korobov_m8 (1963)"
    write(unit=*, fmt="(A)", advance="no") &
      "Enter the name of the periodizer you want to use: "
    read *, periodizer_name

  end subroutine ask_periodizer_name


  ! Ask for the name of a weightfunction
  !
  subroutine ask_weightfunction_name(weightfunction_name)
  
    character(len=*), intent(out) :: weightfunction_name

    write(unit=*, fmt="(A)") "Available weightfunctions:"
    write(unit=*, fmt="(A)") "  IMT"
    write(unit=*, fmt="(A)") "  SagSzekeres"
    write(unit=*, fmt="(A)") "  SugiharaMurotaGaussian"
    write(unit=*, fmt="(A)", advance="no") &
      "Enter the name of the weightfunction you want to use: "
    read *, weightfunction_name

  end subroutine ask_weightfunction_name


  ! Ask for the name of a pointset
  !
  subroutine ask_pointset_name(pointset_name)

    character(len=*), intent(out) :: pointset_name

    write(unit=*, fmt="(A)") "Available pointsets:"
    write(unit=*, fmt="(A)") "  halton"
    write(unit=*, fmt="(A)") "  sobol"
    write(unit=*, fmt="(A)") "  weyl"
    write(unit=*, fmt="(A)") "  extensible_lattice"
    write(unit=*, fmt="(A)", advance="no") "Enter the pointset to use: "
    read *, pointset_name

  end subroutine ask_pointset_name


  ! Create a pointset of the specified
  ! type and size and put it into x.
  !
  subroutine create_pointset(pointset_name, x, startindex)
    character(len=*), intent(in)                :: pointset_name
    real(kind=qp), dimension(:,:), intent(out)  :: x
    integer(kind=i4b), dimension(:), intent(in) :: startindex

    type(soboltype)                              :: mysoboltype
    real(kind=qp), dimension(:), allocatable     :: irrationals
    integer(kind=i4b)                            :: i
    integer(kind=i4b), dimension(:), allocatable :: z

    integer(kind=i4b) :: N, s

    N = size(x, dim=1)
    s = size(x, dim=2)

    select case (pointset_name)

      case ("halton")

        call init_halton(s, startindex)
        do i = 1,N
          call next_halton(x(i,:))
        end do

      case ("sobol")

        mysoboltype%poly_order = "JoeKuo"
        mysoboltype%dirnumbers = "JoeKuo"
        mysoboltype%use_ones = .true.
        mysoboltype%use_antonov_saleev = .true.
        call init_sobol(N, s, startindex, mysoboltype)
        do i=1,N
          call next_sobol(x(i,:))
        end do

      case ("weyl")

        allocate(irrationals(s))
        irrationals = sqrt(real(primes(s), kind=qp))
        call init_weyl(s, irrationals, startindex)
        do i=1,N
          call next_weyl(x(i,:))
        end do
        deallocate(irrationals)

      case ("extensible_lattice")

        allocate(z(s))
        call get_vector_nuyens_kuo(s, z)
        call init_extensible_lattice(s, z)
        do i=1,N
          call next_extensible_lattice_gray(x(i,:))
        end do
        deallocate(z)

      case default

        write(unit=*, fmt="(A)") "ERROR: uknown pointset!"

    end select
  
  end subroutine create_pointset

end module mod_cli
