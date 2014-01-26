! Module containing some numeric and character constants.

module mod_constants

  use numeric_kinds

  private

  ! Numeric constants
  real(kind=qp), parameter, public :: pi = &
    3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865133_qp

  ! Character constants
  !character(len=*), parameter, public :: backslash = "\\"(1:1)
  character(len=*), parameter, public :: backslash = char(92)

end module mod_constants
