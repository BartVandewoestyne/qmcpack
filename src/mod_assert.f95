! Module for debugging, testing and error checking that uses an assert-like
! mechanism.

module mod_assert

  use numeric_kinds

  private

  public :: init_assert
  public :: assert
  public :: get_nb_assert_errors
  public :: increase_nb_assert_errors

  private :: assert_qp
  private :: assert_array_qp
  private :: assert_int
  private :: assert_array_int
  private :: assert_logical

  integer(kind=i4b), private           :: nb_assert_errors
  real(kind=qp), parameter, private    :: default_tol = epsilon(1.0_qp)
  character(len=*), parameter, private :: INDENT = "     "

  interface assert
    module procedure assert_qp, assert_array_qp, assert_int, assert_array_int, assert_logical
  end interface assert

contains

    subroutine init_assert()

      nb_assert_errors = 0

    end subroutine init_assert


    subroutine increase_nb_assert_errors(i)
      integer(kind=i4b), intent(in) :: i

      nb_assert_errors = nb_assert_errors + i

    end subroutine increase_nb_assert_errors


    subroutine assert_int(a, b)
      integer, intent(in) :: a, b

      if (a /= b) then
        print *, ""
        write(unit=*, fmt="(A, I0)") INDENT//"ERROR: found    : ", a
        write(unit=*, fmt="(A, I0)") INDENT//"       expecting: ", b
        print *, ""
        nb_assert_errors = nb_assert_errors + 1
      end if

    end subroutine assert_int
    

    subroutine assert_array_int(a, b)
      integer(kind=i4b), dimension(:), intent(in) :: a, b

      if ( any(a /= b) ) then

        print *, ""
        write(unit=*, fmt="(A)") INDENT//"ERROR:"
        write(unit=*, fmt="(A)") INDENT//"   found    :"
        write(unit=*, fmt="(A, I0)") INDENT, a
        write(unit=*, fmt="(A)") INDENT//"    expecting:"
        write(unit=*, fmt="(A, I0)") INDENT, b
        print *, ""

        nb_assert_errors = nb_assert_errors + 1

      end if

    end subroutine assert_array_int
    

    subroutine assert_qp(a, b, max_rel_err)
      real(kind=qp), intent(in)           :: a, b
      real(kind=qp), intent(in), optional :: max_rel_err

      real(kind=qp) :: tol

      if (present(max_rel_err)) then
        tol = max_rel_err
      else
        tol = epsilon(b)
      end if

      if ( abs(a-b) > tol*max(abs(a),abs(b)) ) then

        print *, ""
        write(unit=*, fmt="(A)") INDENT//"ERROR:"
        print *, INDENT//"    found    : ", a
        print *, INDENT//"    expecting: ", b
        print *, ""

        nb_assert_errors = nb_assert_errors + 1

      end if

    end subroutine assert_qp


    subroutine assert_array_qp(a, b, max_rel_err)
      real(kind=qp), dimension(:), intent(in) :: a, b
      real(kind=qp), intent(in), optional     :: max_rel_err

      real(kind=qp) :: tol

      if (present(max_rel_err)) then
        tol = max_rel_err
      else
        tol = epsilon(b)
      end if

      if ( any(abs(a-b) > tol*(max(abs(a),abs(b)))) ) then

        print *, ""
        write(unit=*, fmt="(A)") INDENT//"ERROR: "
        write(unit=*, fmt="(A)") INDENT//"  found:"
        print *, INDENT, a
        write(unit=*, fmt="(A)") INDENT//"  expecting:"
        print *, INDENT, b
        print *, ""

        nb_assert_errors = nb_assert_errors + 1

      end if

    end subroutine assert_array_qp
    

    subroutine assert_logical(a, b)
      logical, intent(in) :: a, b

      if (.not. (a .eqv. b)) then
        print *, INDENT//"ERROR: found ", a, ", expecting ", b
        nb_assert_errors = nb_assert_errors + 1
      end if

    end subroutine assert_logical


    function get_nb_assert_errors() result (res)
      integer(kind=i4b) :: res

      res = nb_assert_errors

    end function get_nb_assert_errors

end module mod_assert
