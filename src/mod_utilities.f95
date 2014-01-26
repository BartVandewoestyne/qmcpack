! Module containing some general-purpose functions and subroutines.
!
! References:
!
!   [1] `ALGORITHM 659: Implementing Sobol's Quasirandom Sequence
!       Generator', Bratley, Paul and Fox, L. Bennett, 1988, page 93.
!
!   [2] `Computational Investigations of Low-Discrepancy Sequences in
!       Simulation Algorithms for Bayesian Networks', Cheng, Jian and
!       Druzdel, Marek J., UAI-2000.
!
!   [3] Halton's 1960 paper also has a formula.
!
!   [4] Donald E. Knuth, `The Art of Computer Programming, Volume 2:
!       Seminumerical Algorithms, Third Edition', ISBN 0-201-89684-2.
!
!   [5] Press, William H., Teukolsky, Saul A., Vetterling, William T.,
!       Flannery, Brian P., `Numerical Recipes in Fortran 77', second
!       edition, Cambridge University Press, 1992, ISBN 0-521-43064-X.
!
!   [6] `Random Number Generation and Monte Carlo Methods', second
!       edition, James E. Gentle, Springer 2003, ISBN 0-387-00178-6.
!
module mod_utilities

  use numeric_kinds
  use mod_number_theory

  private

  public :: get_nb_digits
  public :: find_rightmost_zerobitpos
  public :: btest_real
  public :: graycode
  public :: b_ary_graycode
  public :: randperm
  public :: random_integer
  public :: frac_part
  public :: get_size_n_vector
  public :: create_n_vector
  public :: binomial_coefficients
  public :: nchoosek
  public :: search_permutation_polynomial
  public :: resample_indices

contains


  ! Return the number of digits of n in base b.  If n is zero, then return 1
  ! for the number of digits.
  !
  function get_nb_digits(n, b) result (m)

    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b), intent(in) :: b
    integer(kind=i4b)             :: m

    integer(kind=i4b) :: temp

    ! For the practical ranges of N we consider, this is probably the preferred
    ! method when working in bases > 2.  The amount of work for this routine
    ! increases with the value of N, so asymptotically it is slower than O(N).
    ! However, for the practical ranges of N we consider, the advantage of this
    ! method is definitely showing up if N is quite low and b is high.
    m = 1
    temp = n
    do
      temp = temp/b  ! we rely on integer division here!
      if (temp == 0) then
        exit
      end if
      m = m+1
    end do

    ! For the practical ranges of N we consider, this is probably the preferred
    ! method when working in base 2.  For each value of N, the amount of work
    ! will be the same (calculating two LOG's and a CEILING) so this is O(N).
    !if (n==0) then
    !  m = 1
    !else
    !  m = ceiling( log(real(n+1, kind=dp))/log(real(b, kind=dp)) , kind=i4b)
    !end if
    
  end function get_nb_digits


  ! Find the rightmost zero bit-position in the base 2 representation of n.
  ! The least significant bit has position 1 (not zero!).
  !
  ! Note:
  !   * This routine is extensively being used when generating Sobol points
  !     using the Antonov-Saleev technique.  It should be optimized as much
  !     as possible!
  !
  function find_rightmost_zerobitpos(n) result (pos)

    ! Note:
    !   I have two ways of implementing this.  One method uses the BTEST
    !   intrinsic and the other one uses the IAND intrinsic.  Which one is
    !   faster seems to depend on the compiler...

    ! The variable named `temp' is only necessary when using method 2.
    ! integer(kind=i4b) :: temp

    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b)             :: pos



                        !!!!!!!!!!!!
                        ! METHOD 1 !
                        !!!!!!!!!!!!

    ! This appeared to be the fasted method for F and f95 (with -O4)
    pos = 0
    do
      if (.not. btest(n, pos)) then
        exit
      end if
      pos = pos + 1
    end do
    ! the BTEST intrinsic uses an index of 0 for the least significant digit
    ! while we start from 1 for our indexing, so add an extra 1 here.
    pos = pos + 1

                        !!!!!!!!!!!!
                        ! METHOD 2 !
                        !!!!!!!!!!!!

    ! This method appeared to be a little faster when compiling with the
    ! Intel(R) Fortran Compiler for 32-bit applications, Version 7.0.
    !pos = 1
    !temp = n
    !do
    !  if (iand(temp, 1) == 0) exit
    !  temp = temp/2  ! we rely on integer division here!
    !  pos = pos + 1
    !end do

  end function find_rightmost_zerobitpos


  ! For a number x in [0, 1], return true if the bit in the base 2
  ! representation of x at position k after the binary point is one, return
  ! false otherwise.  Positions after the binary point start counting at k=-1
  ! and run downwards as -2, -3, ...
  !
  ! Example:
  !
  !   Assuming x = 0.6875 (base 10) = 0.1011 (base 2)
  !
  !       btest_real(0.6875, -1) will return TRUE
  !       btest_real(0.6875, -2) will return FALSE
  !       btest_real(0.6875, -3) will return TRUE
  !       btest_real(0.6875, -4) will return TRUE
  !
  function btest_real(x, k) result (res)

    real(kind=qp), intent(in)     :: x
    integer(kind=i4b), intent(in) :: k
    logical                       :: res

    if (x<0 .or. x>1) then
      print *, "ERROR: btest_real(x, k) only operates on x-values in [0,1]"
    end if

    if (modulo(floor(x*(2**(-k))), 2) == 1) then
      res = .true.
    else
      res = .false.
    end if
    
  end function btest_real


  ! Return the fractional part of x.
  !
  elemental function frac_part(x) result (y)

    real(kind=qp), intent(in) :: x
    real(kind=qp)             :: y
  
    y = x - floor(x)

  end function frac_part


  ! Return the Gray code of n.
  !
  elemental function graycode(n) result (y)

    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b)             :: y

    y = ieor(n, ishft(n, -1))

  end function graycode


  ! Return the b-ary graycode of n.  If n is represented in base b as
  !
  !             n = n_m*b^{m-1} + ... + n_2*b + n_1
  !
  ! then the b-ary graycode of n is
  !
  !             n = g_m*b^{m-1} + ... + g_2*b + g_1
  !
  ! and the digits g_1,...,g_m are obtained via (for example with m=4)
  !
  !     [ g_1 ]   [1    -1     0     0] [ n_1 ]
  !     [     ]   [                   ] [     ]
  !     [ g_2 ]   [0     1    -1     0] [ n_2 ]
  !     [     ] = [                   ] [     ]  MODULO B
  !     [ g_3 ]   [0     0     1    -1] [ n_3 ]
  !     [     ]   [                   ] [     ]
  !     [ g_4 ]   [0     0     0     1] [ n_4 ]
  !
  !
  elemental function b_ary_graycode(n, b) result (y)

    integer(kind=i4b), intent(in) :: n
    integer(kind=i4b), intent(in) :: b
    integer(kind=i4b)             :: y

    integer(kind=i4b) :: temp
    integer(kind=i4b) :: n_prev, n_curr
    integer(kind=i4b) :: b_power

    temp = n
    n_prev = modulo(temp, b)
    b_power = 1
    y = 0
    do
      if (n_prev == temp) then
        exit
      end if
      temp = temp/b
      n_curr = modulo(temp, b)
      y = y + modulo(n_prev-n_curr, b)*b_power
      n_prev = n_curr
      b_power = b_power*b
    end do
    y = y + n_prev*b_power

  end function b_ary_graycode


  ! Return a random permutation of {1,...,N}.
  !
  ! See Algorithm 3.4.2P on page 145 of [4].  This subroutine uses
  ! the algorithm without the exchange operation as described in the text
  ! at the bottom of page 145.
  !
  ! Note:
  !   * Other algorithms can be found at:
  !
  !       1) http://ftp.cac.psu.edu/pub/ger/fortran/hdk/byterand.for
  !       2) http://www.icsi.berkeley.edu/~storn/DE_FORTRAN90.f90
  !       3) In the book [6] on page 220.
  !
  subroutine randperm(N, p)

    integer(kind=i4b), intent(in)                :: N
    integer(kind=i4b), dimension(:), intent(out) :: p

    integer(kind=i4b) :: j, k
    real(kind=qp)     :: u

    p = 0

    do j=1,N

      call random_number(u)
      k = floor(j*u) + 1

      p(j) = p(k)
      p(k) = j

    end do

  end subroutine randperm


  ! Returns a random integer from the closed interval 0..max_int-1
  !
  subroutine random_integer(max_int, res)

    integer(kind=i4b), intent(in)  :: max_int
    integer(kind=i4b), intent(out) :: res

    real(kind=qp) :: x

    call random_number(x)

    res = floor(max_int*x)

  end subroutine random_integer


  ! Return the size of a rank 1 array containing some values that are acceptable
  ! for a log-log plot.
  !
  function get_size_n_vector(N, b) result(vec_size)
    integer(kind=i4b), intent(in) :: N
    integer(kind=i4b), intent(in) :: b
    integer(kind=i4b)             :: vec_size

    vec_size = floor(log(real(N, kind=qp))/log(real(b, kind=qp)))

  end function get_size_n_vector


  ! Return a rank 1 array containing some values that are acceptable
  ! for a log-log plot.
  !
  subroutine create_n_vector(b, vec)
    integer(kind=i4b), intent(in)                :: b
    integer(kind=i4b), dimension(:), intent(out) :: vec
    
    integer(kind=i4b) :: i

    vec = (/ (b**i, i=1,size(vec)) /)

  end subroutine create_n_vector


  ! Return a rank-1 array containing at position k the binomial
  ! coefficient binomial(n, k-1).  The coefficients are calculated
  ! using the recurrence relation:
  !
  !             binomial(n, k+1) = (n-k)/(k+1)*binomial(n, k)
  !
  ! See also [1], page 209.
  !
  subroutine binomial_coefficients(n, coef)
    integer(kind=i4b), intent(in)                 :: n
    integer(kind=i4b), dimension(0:), intent(out) :: coef

    integer(kind=i4b) :: k

    coef(0) = 1
    do k=0,n-1
      coef(k+1) = (n-k)*coef(k)/(k+1)
    end do

  end subroutine binomial_coefficients


  ! Return the binomial coefficient
  !
  !       (n)       n!
  !       ( ) = -----------
  !       (k)   k! (n - k)!
  !
  ! Note that the current (naive) implementation is prone to overflow errors,
  ! severely restricting the range of input values.
  !
  function nchoosek(n, k) result (c)
    integer(kind=i4b), intent(in) :: n, k
    integer(kind=i4b)             :: c

    c = integer_factorial(n)/(integer_factorial(k)*(integer_factorial(n-k)))
  end function nchoosek


  subroutine search_permutation_polynomial(p, u, v, w)
    integer(kind=i4b), intent(in)  :: p
    integer(kind=i4b), intent(out) :: u, v, w

    integer(kind=i4b)                   :: i, j, k, x
    integer(kind=i4b), dimension(0:p-1) :: res
    logical                             :: found

    
    do i=0,p-1
      do j=0,p-1
        do k=0,p-1

          found = .true.

          ! Check if polynomial returns different results for all x
          do x=0,p-1
            res(x) = modulo(modulo(u*x*x, p) + modulo(v*x, p) + w, p)
            if ( any(res(0:x-1) == res(x)) ) then
              found = .false.
              exit
            end if
          end do

          if (found) then
            u = i
            v = j
            w = k
            print *, "u = ", i, ", v = ", j, ", w = ", k
            return
          end if

        end do
      end do
    end do

    if (.not. found) then
      print *, "found no good polynomial..."
    end if

  end subroutine search_permutation_polynomial


  ! Resample the indices from i that have the weights from w and start
  ! from the uniform numbers in u to do the inverse CDF.
  ! The new samples are stored in j.
  !
  ! TODO: CHECK THIS (bartv, 19 april 2007)
  !
  subroutine resample_indices(w, u, ind_out)
    real(kind=qp), dimension(:), intent(in)      :: w
    real(kind=qp), dimension(:), intent(in)      :: u
    integer(kind=i4b), dimension(:), intent(out) :: ind_out

    integer(kind=i4b) :: i, j, nb_samples_out
    real(kind=qp)     :: cumsum

    nb_samples_out = size(u)

    cumsum = 0.0_qp
    j = 0
    do i = 1,nb_samples_out
      do 
        if (cumsum > u(i)) then
          exit
        end if
        j = j + 1
        cumsum = cumsum + w(j)
      end do
      ind_out(i) = j
    end do

  end subroutine resample_indices


end module mod_utilities
