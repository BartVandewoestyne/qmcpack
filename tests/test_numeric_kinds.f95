program test_numeric_kinds

  use qmcpack

  real              :: x_real
!  double precision  :: x_double_precision ! OBSOLETE!
  real(kind=sp)     :: x_sp
  real(kind=dp)     :: x_dp
  real(kind=qp)     :: x_qp
  real(kind=qp2)    :: x_qp2
  integer(kind=i1b) :: x_i1b
  integer(kind=i2b) :: x_i2b
  integer(kind=i4b) :: x_i4b
  integer(kind=i8b) :: x_i8b
  integer(kind=i4b) :: nb_errors

  call show_test_header("NUMERIC_KINDS")

  nb_errors = 0

  ! This is simply to not have compiler warnings about used but not set
  ! variables.
  x_real = 1.0
!  x_double_precision = 1.d0 ! OBSOLETE
  x_sp = 1.0_sp
  x_dp = 1.0_dp
  x_qp = 1.0_qp
  x_qp2 = 1.0_qp2
  x_i1b = 1
  x_i2b = 1
  x_i4b = 1
  x_i8b = 1

  write(unit=*, fmt=*)
  write(unit=*, fmt="(A)") "  ### INTEGER INFORMATION ###"
  if (i1b==-1) then
    write(unit=*, fmt="(A)") "  i1b is not available on this computer"
  end if
  if (i2b==-1) then
    write(unit=*, fmt="(A)") "  i2b is not available on this computer"
  end if
  if (i4b==-1) then
    write(unit=*, fmt="(A)") "  i4b is not available on this computer"
  end if
  if (i8b==-1) then
    write(unit=*, fmt="(A)") "  i8b is not available on this computer"
  end if
  write(unit=*, fmt="(A, I0.0)") "  Radix for i1b is: ", radix(x_i1b) 
  write(unit=*, fmt="(A, I0.0)") "  Radix for i2b is: ", radix(x_i2b) 
  write(unit=*, fmt="(A, I0.0)") "  Radix for i4b is: ", radix(x_i4b) 
  write(unit=*, fmt="(A, I0.0)") "  Radix for i8b is: ", radix(x_i8b) 
  write(unit=*, fmt="(A, I0.0, A, I0.0)") "  Digits in radix ", radix(x_i1b), " for i1b is: ", digits(x_i1b) 
  write(unit=*, fmt="(A, I0.0, A, I0.0)") "  Digits in radix ", radix(x_i2b), " for i2b is: ", digits(x_i2b) 
  write(unit=*, fmt="(A, I0.0, A, I0.0)") "  Digits in radix ", radix(x_i4b), " for i4b is: ", digits(x_i4b) 
  write(unit=*, fmt="(A, I0.0, A, I0.0)") "  Digits in radix ", radix(x_i8b), " for i8b is: ", digits(x_i8b) 
  write(unit=*, fmt="(A, I0.0)") "  Largest integer for i1b is: ", huge(x_i1b)
  write(unit=*, fmt="(A, I0.0)") "  Largest integer for i2b is: ", huge(x_i2b)
  write(unit=*, fmt="(A, I0.0)") "  Largest integer for i4b is: ", huge(x_i4b)
  write(unit=*, fmt="(A, I0.0)") "  Largest integer for i8b is: ", huge(x_i8b)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal exponent range for i1b is: ", range(x_i1b)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal exponent range for i2b is: ", range(x_i2b)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal exponent range for i4b is: ", range(x_i4b)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal exponent range for i8b is: ", range(x_i8b)

  write(unit=*, fmt=*)
  write(unit=*, fmt="(A)") "  ### REAL INFORMATION ###"
  if (sp==-1) then
    write(unit=*, fmt="(A)") "  sp is not available on this computer"
  end if
  if (dp==-1) then
    write(unit=*, fmt="(A)") "  dp is not available on this computer"
  end if
  if (qp==-1) then
    write(unit=*, fmt="(A)") "  qp is not available on this computer, so"
    write(unit=*, fmt="(A)") "  making it dp."
  end if
  if (qp2==-1) then
    write(unit=*, fmt="(A)") "  qp2 is not available on this computer, so"
    write(unit=*, fmt="(A)") "  making it dp."
  end if
  write(unit=*, fmt="(A, I0.0)") "  Radix for sp is: ", radix(x_sp) 
  write(unit=*, fmt="(A, I0.0)") "  Radix for dp is: ", radix(x_dp) 
  write(unit=*, fmt="(A, I0.0)") "  Radix for qp is: ", radix(x_qp) 
  write(unit=*, fmt="(A, I0.0)") "  Radix for qp2 is: ", radix(x_qp2) 
  write(unit=*, fmt="(A, I0.0)") "  Digits for sp is: ", digits(x_sp) 
  write(unit=*, fmt="(A, I0.0)") "  Digits for dp is: ", digits(x_dp) 
  write(unit=*, fmt="(A, I0.0)") "  Digits for qp is: ", digits(x_qp) 
  write(unit=*, fmt="(A, I0.0)") "  Digits for qp2 is: ", digits(x_qp2) 
  write(unit=*, fmt="(A, I0.0)") "  Maximum exponent e_max for sp is: ", maxexponent(x_sp) 
  write(unit=*, fmt="(A, I0.0)") "  Maximum exponent e_max for dp is: ", maxexponent(x_dp) 
  write(unit=*, fmt="(A, I0.0)") "  Maximum exponent e_max for qp is: ", maxexponent(x_qp) 
  write(unit=*, fmt="(A, I0.0)") "  Maximum exponent e_max for qp2 is: ", maxexponent(x_qp2) 
  write(unit=*, fmt="(A, I0.0)") "  Minimum exponent e_min for sp is: ", minexponent(x_sp) 
  write(unit=*, fmt="(A, I0.0)") "  Minimum exponent e_min for dp is: ", minexponent(x_dp) 
  write(unit=*, fmt="(A, I0.0)") "  Minimum exponent e_min for qp is: ", minexponent(x_qp) 
  write(unit=*, fmt="(A, I0.0)") "  Minimum exponent e_min for qp2 is: ", minexponent(x_qp2) 
  write(unit=*, fmt="(A, ES60.50e5)") "  Epsilon for sp is: ", epsilon(x_sp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Epsilon for dp is: ", epsilon(x_dp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Epsilon for qp is: ", epsilon(x_qp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Epsilon for qp2 is: ", epsilon(x_qp2)
  write(unit=*, fmt="(A, ES60.50e5)") "  Smallest normalized positive real number for sp is: ", tiny(x_sp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Smallest normalized positive real number for dp is: ", tiny(x_dp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Smallest normalized positive real number for qp is: ", tiny(x_qp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Smallest normalized positive real number for qp2 is: ", tiny(x_qp2)
  write(unit=*, fmt="(A, ES60.50e5)") "  Largest real for sp is: ", huge(x_sp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Largest real for dp is: ", huge(x_dp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Largest real for qp is: ", huge(x_qp)
  write(unit=*, fmt="(A, ES60.50e5)") "  Largest real for qp2 is: ", huge(x_qp2)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal precision for real is: ", precision(x_real)
!  print *, "Equivalent decimal precision for double_precision is: ", precision(x_double_precision) ! OBSOLETE
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal precision for sp is: ", precision(x_sp)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal precision for dp is: ", precision(x_dp)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal precision for qp is: ", precision(x_qp)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal precision for qp2 is: ", precision(x_qp2)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal exponent range for sp is: ", range(x_sp)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal exponent range for dp is: ", range(x_dp)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal exponent range for qp is: ", range(x_qp)
  write(unit=*, fmt="(A, I0.0)") "  Equivalent decimal exponent range for qp2 is: ", range(x_qp2)

  call show_test_summary(nb_errors)

end program test_numeric_kinds
