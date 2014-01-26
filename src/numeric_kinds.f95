module numeric_kinds

  ! Named constants for 1, 2, 4 and 8 byte signed integers.  On most
  ! architectures, integers are represented in two's complement notation,
  ! so they come from the closed interval
  !
  !                     [-2^{n-1}...2^{n-1}-1]
  !
  ! where n is the number of bits used.  Therefore, the maximum absolute
  ! value that can occur is 2^{n-1}, which leads to the following table:
  !
  !   1 byte  => max 2^7  = 1.28e+02                  so 10^2 is assured
  !   2 bytes => max 2^15 = 3.2768e+04,               so 10^4 is assured
  !   4 bytes => max 2^31 = 2.147483648e+09,          so 10^9 is assured
  !   8 bytes => max 2^63 = 9.223372036854775808e+18, so 10^18 is assured
  !
  integer, parameter, public :: i1b = selected_int_kind(2)
  integer, parameter, public :: i2b = selected_int_kind(4)
  integer, parameter, public :: i4b = selected_int_kind(9)
  integer, parameter, public :: i8b = selected_int_kind(18)

  ! Named constants for single, double and quadruple precision reals.
  ! If quadruple precision is not available, make it double precision.
  integer, parameter, public :: sp = kind(1.0)
  integer, parameter, public :: dp = selected_real_kind(2*precision(1.0_sp))
  integer, parameter, public :: qp_preferred = &
                                        selected_real_kind(2*precision(1.0_dp))

  ! Uncomment this for the normal situation...
  integer, parameter, public :: qp = (1+sign(1,qp_preferred))/2*qp_preferred+ &
                                     (1-sign(1,qp_preferred))/2*dp

  ! Use these for the precision test...
  !integer, parameter, public :: qp = dp
  !integer, parameter, public :: qp2 = dp
  integer, parameter, public :: qp2 = (1+sign(1,qp_preferred))/2*qp_preferred+ &
                                      (1-sign(1,qp_preferred))/2*dp


  ! Some experimenting resulted in the following table:
  !
  ! For Intel(R) Pentium(R) 4 CPU 2.40GHz:
  !
  ! Compiler
  !   => <best integer kind>
  !   => <best real kind>
  !-------------------------------------------------------------------------
  ! Fortran Company/NAG F compiler Release 20031017
  !  => i8b
  !  => dp
  !
  ! G95 (GCC 4.0.1 (g95!) Aug 17 2005)
  !  => i8b
  !  => dp
  !
  ! NAGWare Fortran 95 compiler Release 5.0(361)
  !  => i8b
  !  => dp
  !
  ! GNU Fortran 95 (GCC 4.1.0 20050818 (experimental))
  !  =>
  !  =>
  !
  ! Intel(R) Fortran Compiler for 32-bit applications, Version 7.0   Build
  ! 20021028Z
  !  => 
  !  =>
  !
  ! Intel(R) Fortran Compiler for 32-bit applications, Version 9.0    Build
  ! 20050430Z
  !  => i8b
  !  => qp

end module numeric_kinds
