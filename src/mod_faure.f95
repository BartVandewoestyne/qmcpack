! Module containing different implementations of the Faure and scrambled Faure
! sequence.  Currently available Faure scramblings are
!
!  * "Faure"                            Original Faure (See [1]) 
!  * "GFaure"                           GFaure
!  * "RandomLinearScrambling"           Random Linear scrambling
!  * "RandomLinearDigitScrambling"      Random Linear Digit scrambling
!  * "LeftIbinomial"                    Left I-binomial scrambling
!  * "StripedMatrix"                    Striped Matrix scrambling
!  * "InversivePre"                     Inversive + Random Linear Digit
!                                       Scrambling
!  * "InversivePost"                    Random Linear Digit Scrambling
!                                       + Inversive
!
! To be implemented in the future:
!
!  * RightIbinomial                     Right I-binomial scrambling
!  * FaureTezuka                        ??? Multiplication on the right ???
!
! Note that there are actually three versions of the Faure sequence:
!
!  * Original Faure (See [1])
!  * Generalized Faure by Faure (See [?])
!  * Generalized Faure by Tezuka and Tokuyama (GFaure in Tezuka's book?  See
!    also slides of Hongmei Chi)
!
! Notes:
!   * currently, only fau_n(1) is used instead of the whole
!     array.  Probably we will be able to replace the array by a
!     scalar because a scalar suffices.
!
! References:
!
!   [1] Faure, H.  `Discrepance de suites associees a un systeme de
!       numeration (en dimension s)',  Acta Arithmetica XLI (1982), 337-351.
!
!   [2] B. L. Fox.  `Algorithm 647: Implementation and relative efficiency
!       of quasirandom sequence generator.'  ACM Trans. Math.
!       Software, 12:362-376, 1986.
!       => The generator of this paper only allows generation of points
!          in 40 dimensions.
!       => The QS-values are hardcoded up to the first 40 dimensions.
!
!   [3] Thiemard, Eric.  `Economic Generation of Low-Discrepancy Sequences
!       with a b-ary Gray Code', Technical Report RO981201, Departement de
!       Mathematiques, Ecole Polytechnique Federale de Lausanne, CH-1015
!       Lausanne, Switzerland. (See http://roso.epfl.ch/papers/grayfaure/)
!
!   [4] Hee Sun Hong, Fred J. Hickernell, `Algorithm 823: Implementing
!       Scrambled Digital Sequences', ACM Transactions on Mathematical
!       Software, Vol. 29, No. 2, June 2003, Pages 95--109.
!
!   [8] `Quasi-Monte Carlo Methods in Numerical Finance', Corwin Joy, Phelim P.
!       Boyle, Ken Seng Tan, Management Science, Vol. 42, No. 6, June 1996.
!
!
! Interesting links:
!
!   http://www.csit.fsu.edu/~burkardt/f_src/faure/faure.html
!
!   http://mint.sbg.ac.at/desc_SFaure-ori.html
!
!   http://www.vni.com/products/imsl/jmsl/v30/api/com/imsl/stat/FaureSequence.html
!   http://www.iro.umontreal.ca/~simardr/ssj/doc/html/umontreal/iro/lecuyer/hups/FaureSequence.html (based on Fox's implementation?)
!
! TODO:
!   * Check the points returned by these routines for correctness.

module mod_faure

  use numeric_kinds
  use mod_primes
  use mod_radical_inverse
  use mod_utilities
  use mod_scrambling
  use mod_debug
  use mod_number_theory

  private

  public :: init_faure_bratley_fox
  public :: next_faure_bratley_fox
  public :: init_faure
  public :: next_faure
  public :: init_faure_gray
  public :: next_faure_gray
  public :: free_faure

  private :: construct_pascal_power_matrices
  private :: print_matrices
  private :: print_vectors
  private :: generate_scrambling
  private :: apply_scrambling
  private :: next_faure_left_matrix_scramble
  private :: next_faure_inversive_pre
  private :: next_faure_inversive_post

  integer(kind=i4b), private                                :: s
  integer(kind=i4b), private                                :: qs
  real(kind=qp), private                                    :: inverse_qs
  character(len=50), private                                :: scrambletype
  integer(kind=i4b), dimension(:,:), allocatable, private   :: coef
  integer(kind=i4b), dimension(:,:,:), allocatable, private :: pascal_powers
  integer(kind=i4b), dimension(:,:,:), allocatable, private :: c_eff
  integer(kind=i4b), dimension(:,:), allocatable, private   :: shift_eff
  integer(kind=i4b), dimension(:,:,:), allocatable, private :: gfaure_a
  integer(kind=i4b), dimension(:,:,:), allocatable, private :: owen_matrix
  integer(kind=i4b), dimension(:,:), allocatable, private   :: owen_shift
  integer(kind=i4b), dimension(:,:), allocatable, private   :: graycode_bits
  integer(kind=i4b), dimension(:), allocatable, private     :: n_digits
  integer(kind=i4b), dimension(:), allocatable, private     :: temp_digits
  integer(kind=i4b), dimension(:), allocatable, private     :: fau_n
  integer(kind=i4b), dimension(:), allocatable, private     :: step
  integer(kind=i4b), private                                :: max_digits

contains


  ! Initialize the data that is necessary to generate a scrambled Faure
  ! sequence in non-Graycode ordering.
  !
  subroutine init_faure(maxpoints, init_s, init_scrambletype, init_base, init_startindex, init_step)
    integer(kind=i4b), intent(in)                         :: maxpoints
    integer(kind=i4b), intent(in)                         :: init_s
    character(len=*), intent(in), optional                :: init_scrambletype
    integer(kind=i4b), intent(in), optional               :: init_base
    integer(kind=i4b), dimension(:), intent(in), optional :: init_startindex
    integer(kind=i4b), dimension(:), intent(in), optional :: init_step

    s = init_s

    if (present(init_base)) then
      if (init_base >= s) then
        qs = init_base
      else
        write(unit=*, fmt="(A)") &
          "ERROR: base for Faure sequence must be greater than or equal to the dimension!"
        stop
      end if
    else
      ! Determine the first prime larger or equal to s
      qs = prime_ge(s)
    end if

    ! Set the n-value for which we start generating points (this allows
    ! to skip an initial segment of the (scrambled) Faure sequence).
    if (allocated(fau_n)) then
      deallocate(fau_n)
    end if
    allocate(fau_n(s))
    if (present(init_startindex)) then
      fau_n = init_startindex
    else
      fau_n = qs**4-1
    end if

    ! Set a possible step value
    if (allocated(step)) then
      deallocate(step)
    end if
    allocate(step(s))
    if (present(init_step)) then
      step = init_step
    else
      step = 1
    end if

    ! Set the type of scrambling to use
    if (present(init_scrambletype)) then
      scrambletype = init_scrambletype
    else
      scrambletype = "None"
    end if

    ! Calculate an upper bound on the number of digits we are ever gonna need.
    ! Add 1 for safety???
    max_digits = get_nb_digits(maxval(fau_n+(maxpoints-1)*step), qs)+1

    ! Allocate and construct the `Pascal Power matrices'
    if (allocated(pascal_powers)) then
      deallocate(pascal_powers)
    end if
    allocate(pascal_powers(0:s-1, 0:max_digits-1, 0:max_digits-1))
    pascal_powers = 0
    call construct_pascal_power_matrices(pascal_powers)

    ! Construct effective generator matrix and shift vector
    ! depending on the type of scrambling used

    if (allocated(c_eff)) then
      deallocate(c_eff)
    end if
    allocate(c_eff(s, 0:max_digits-1, 0:max_digits-1))

    if (allocated(shift_eff)) then
      deallocate(shift_eff)
    end if
    allocate(shift_eff(s, 0:max_digits-1))

    if (allocated(n_digits)) then
      deallocate(n_digits)
    end if
    allocate(n_digits(max_digits))

    if (allocated(temp_digits)) then
      deallocate(temp_digits)
    end if
    allocate(temp_digits(max_digits))

    call generate_scrambling(scrambletype)
    call apply_scrambling(scrambletype)

    call print_matrices(c_eff)

  end subroutine init_faure


  ! Return the next point of the Faure sequence, scrambled with
  ! the correct scrambling method.
  !
  subroutine next_faure(x)
    real(kind=qp), dimension(:), intent(out) :: x

    select case (scrambletype)

      case ("None",                        &
            "RandomLinearScrambling",      &
            "RandomLinearDigitScrambling", &
            "StripedMatrix",               &
            "LeftIbinomial",               &
            "GFaure")

        call next_faure_left_matrix_scramble(x)

      case ("InversivePre")

        call next_faure_inversive_pre(x)
        
      case ("InversivePost")

        call next_faure_inversive_post(x)

      case default

        write(unit=*, fmt="(A)") &
          "ERROR: unknown scrambling type for Faure sequence"
        stop

    end select

  end subroutine next_faure


  ! Return the next Faure point, with left matrix scrambling
  ! and a shift added.
  !
  subroutine next_faure_left_matrix_scramble(x)
    real(kind=qp), dimension(:), intent(out) :: x

    integer(kind=i4b) :: i, j

    do i = 1, s

      temp_digits = 0
      n_digits = 0

      ! TODO: check waarom n_digits UNDEFINED is als we dit er in laten.
      call get_digits(fau_n(i), qs, n_digits)

      ! Multiply digit-vector with effective generator matrix
      ! and apply effective shift.
      temp_digits = modulo(modulo(matmul(c_eff(i,:,:), n_digits), qs) &
                              + shift_eff(i,:), qs)

      ! Compute radical inverse value of resulting digit-vector.
      x(i) = dot_product(temp_digits, (/ (1.0_qp/qs**(j), j=1,max_digits) /))

    end do

    fau_n = fau_n + step

  end subroutine next_faure_left_matrix_scramble


  ! Return the next inversive scrambled Faure point, non-Graycode ordering.
  !
  subroutine next_faure_inversive_pre(x)
    real(kind=qp), dimension(:), intent(out) :: x

    integer(kind=i4b) :: i, j

    do i = 1, s

      temp_digits = 0
      n_digits = 0

      call get_digits(fau_n(i), qs, n_digits)

      ! Multiply with the effective generator matrices and add the effective
      ! shift = construct ordinary Faure points.
      !print *, "Effective generator matrix for dimension ", i
      !call print_matrix(c_eff(i,:,:))
      !print *, "Effective shift vector for dimension ", i
      !print *, shift_eff(i,:)
      temp_digits = modulo(modulo(matmul(c_eff(i,:,:), n_digits), qs) &
                                     + shift_eff(i,:), qs)

      ! FIRST: compute the inverse of the non-zero digits.
      where (temp_digits /= 0) 
        temp_digits = inverse_mod(temp_digits, qs)
      end where

      ! SECOND: apply a random linear digit scrambling to the inversed digits.
      !print *, "Scramble matrix for dimension ", i
      !call print_matrix(owen_matrix(i,:,:))
      !print *, "Random shift vector for dimension ", i
      !print *, owen_shift(i,:)
      temp_digits = modulo(modulo(matmul(owen_matrix(i,:,:), temp_digits), qs) &
                                     + owen_shift(i,:), qs)

      ! Compute radical inverse value of resulting digit-vector.
      x(i) = dot_product(temp_digits, (/ (1.0_qp/qs**(j), j=1,max_digits) /))

    end do

    fau_n = fau_n + step

  end subroutine next_faure_inversive_pre


  subroutine next_faure_inversive_post(x)
    real(kind=qp), dimension(:), intent(out) :: x

    integer(kind=i4b) :: i, j

    do i = 1, s

      temp_digits = 0
      n_digits = 0

      call get_digits(fau_n(i), qs, n_digits)

      ! Multiply with the effective generator matrices and add the effective
      ! shift = construct ordinary Faure points.
      temp_digits = modulo(modulo(matmul(c_eff(i,:,:), n_digits), qs) &
                                     + shift_eff(i,:), qs)

      ! FIRST: apply a random linear digit scrambling to the digits.
      temp_digits = modulo(modulo(matmul(owen_matrix(i,:,:), temp_digits), qs) &
                                     + owen_shift(i,:), qs)

      ! SECOND: compute the inverse of the non-zero digits.
      where (temp_digits /= 0)
        temp_digits = inverse_mod(temp_digits, qs)
      end where

      ! Compute radical inverse value of resulting digit-vector.
      x(i) = dot_product(temp_digits, (/ (1.0_qp/qs**(j), j=1,max_digits) /))

    end do

    fau_n = fau_n + step

  end subroutine next_faure_inversive_post


    ! Initialize the data that is necessary to generate a scrambled Faure
    ! sequence in Graycode-ordering.
    !
    ! Possible scrambling types are:
    !
    !   * "None"                        Do not apply any scrambling at all.
    !   * "GFaure"                      Construct the GFaure sequence.
    !   * "RandomLinearScrambling"      The Random Linear Scrambling from Matousek.
    !   * "RandomLinearDigitScrambling" The Random Linear Digit Scrambling form Matousek.
    !   * "FaureTezuka"                 Apply Faure-Tezuka scrambling (Ref???)
    !
    subroutine init_faure_gray(maxpoints, init_s, init_scrambletype, init_startindex, init_step)
      integer(kind=i4b), intent(in)                         :: maxpoints
      integer(kind=i4b), intent(in)                         :: init_s
      character(len=*), intent(in)                          :: init_scrambletype
      integer(kind=i4b), dimension(:), intent(in), optional :: init_startindex
      integer(kind=i4b), dimension(:), intent(in), optional :: init_step

      integer(kind=i4b)                              :: i
      integer(kind=i4b), dimension(init_s)           :: n_gray, temp

      s = init_s

      ! Determine the first prime larger or equal to s
      qs = prime_ge(s)

      ! Set the n-value for which we start generating points (this allows
      ! to skip an initial segment of the (scrambled) Faure sequence).
      if (allocated(fau_n)) then
        deallocate(fau_n)
      end if
      allocate(fau_n(s))
      if (present(init_startindex)) then
        fau_n = init_startindex
      else
        ! Use 0 as default.  This is the same value as is used in
        ! Hickernell's 823 implementation.
        fau_n = 0
      end if

      ! Set a possible step value
      if (allocated(step)) then
        deallocate(step)
      end if
      allocate(step(s))
      if (present(init_step)) then
        step = init_step
      else
        step = 1
      end if

      ! Calculate an upper bound on the number of digits we are ever gonna need.
      ! Add 1 because the graycode-value of the maximum of fau_n can have one more digit?
      max_digits = get_nb_digits(maxval(fau_n+(maxpoints-1)*step), qs)+1

      ! Allocate and construct the `Pascal Power matrices'
      if (allocated(pascal_powers)) then
        deallocate(pascal_powers)
      end if
      allocate(pascal_powers(0:s-1, 0:max_digits-1, 0:max_digits-1))
      pascal_powers = 0
      call construct_pascal_power_matrices(pascal_powers)

      ! Construct effective generator matrix and shift vector
      ! depending on the type of scrambling used

      if (allocated(c_eff)) then
        deallocate(c_eff)
      end if
      allocate(c_eff(s, 0:max_digits-1, 0:max_digits-1))

      if (allocated(shift_eff)) then
        deallocate(shift_eff)
      end if
      allocate(shift_eff(s, 0:max_digits-1))

      call generate_scrambling(init_scrambletype)
      call apply_scrambling(init_scrambletype)

      deallocate(pascal_powers)

      ! Calculate the graycode for the previous n-values.
      n_gray = b_ary_graycode(fau_n-1, qs)

      ! Get graycode bits of n_gray for all dimensions.  Initially, the graycode
      ! is the same for each dimension, but from the moment we start
      ! generating next points, the bits will become different due to
      ! the differences in the powers of the Pascal matrix.
      allocate(graycode_bits(s, max_digits))
      graycode_bits = 0
      temp = n_gray
      i = 1
      do
        ! TODO: HIER NOG ERGENS DIE SHIFT INBRENGEN???
        ! OOK JE GENERATOR MATRIX ER AL IN BETREKKEN???
        graycode_bits(:,i) =  modulo(temp, qs)
        temp = temp/qs
        if (all(temp == 0)) then
          exit
        end if
        i = i+1
      end do

      inverse_qs = 1.0_qp/real(qs, kind=qp)


    end subroutine init_faure_gray


    ! Initialize the data that is necessary to generate the Faure
    ! sequence.  The initialization consists of two things:
    !
    !   1) Determine qs, the first prime that is larger or equal to
    !      the dimension.
    !
    !   2) Set up the matrix C (not the powers of C, because we follow
    !      the technique used in Bratley and Fox for the non-graycode version).
    !   
    subroutine init_faure_bratley_fox(maxpoints, init_s, init_startindex, init_step)

      integer(kind=i4b), intent(in)                         :: maxpoints
      integer(kind=i4b), intent(in)                         :: init_s
      integer(kind=i4b), dimension(:), intent(in), optional :: init_startindex
      integer(kind=i4b), dimension(:), intent(in), optional :: init_step

      integer(kind=i4b) :: i, j

      s = init_s

      ! Determine the first prime larger or equal to s
      qs = prime_ge(s)

      if (allocated(fau_n)) then
        deallocate(fau_n)
      end if
      allocate(fau_n(s))
      if (present(init_startindex)) then
        fau_n = init_startindex
      else
        ! Use qs**4-1 as default.  This is the same value as is used in
        ! the Bratley and Fox implementation.
        fau_n = qs**4-1
      end if

      if (allocated(step)) then
        deallocate(step)
      end if
      allocate(step(s))
      if (present(init_step)) then
        step = init_step
      else
        step = 1
      end if

      ! Calculate an upper bound on the number of digits we are ever gonna need.
      max_digits = get_nb_digits(maxval(fau_n+(maxpoints-1)*step), qs)

      ! Build Pascal matrix (modulo qs) for MAX_DIGITS digits by finding
      ! binomial coefficients modulo qs in a lower-triangular matrix
      ! COEF using the recursion
      !
      !   binom(i, j) = binom(i-1, j) + binom(i-1, j-1)
      !
      ! and where a = b + c implies
      !
      !   mod(a, qs) = mod( mod(b, qs) + mod(c, qs), qs )
      !
      ! Note that if C is the matrix to be used in the generation of the Faure
      ! sequence, then we have that COEF=transpose(C)!!!  This is probably
      ! better because of the column-major ordering of Fortran.  So we have
      !
      !         1     1     1     1     1
      !         0     1     2     3     4
      !    C =  0     0     1     3     6
      !         0     0     0     1     4
      !         0     0     0     0     1 ... and so on ...
      !
      ! but
      !         1     0     0     0     0
      !         1     1     0     0     0
      ! COEF =  1     2     1     0     0
      !         1     3     3     1     0
      !         1     4     6     4     1 ... and so on ...

      ! TODO: check if coef should be initialized to zero.
      if (allocated(coef)) then
        deallocate(coef)
      end if
      allocate(coef(0:max_digits-1, 0:max_digits-1))
      ! Set first column and diagonal to 1...
      coef(0, 0) = 1
      do j=1, max_digits-1
        coef(j, 0) = 1
        coef(j, j) = 1
      end do
      ! ... and fill up the rest of the elements using the recursion cited
      ! above.
      do j=1, max_digits-1
        do i=j+1, max_digits-1
          coef(i, j) = modulo( coef(i-1, j) + coef(i-1, j-1), qs )
        end do
      end do

      inverse_qs = 1.0_qp/real(qs, kind=qp)

    end subroutine init_faure_bratley_fox


    subroutine next_faure_bratley_fox(x)

      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b), dimension(:), allocatable :: ytemp
      integer(kind=i4b)                            :: ztemp
      integer(kind=i4b)                            :: temp
      integer(kind=i4b)                            :: i, j, k
      integer(kind=i4b)                            :: nb_digits
      real(kind=qp)                                :: radinv
      real(kind=qp)                                :: b_factor


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Find the component of the first dimension. !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! fau_n has a representation in base qs of the form
      !
      !       nb_digits - 1
      !          -----
      !           \               j
      !   fau_n =  )   ytemp(j) qs
      !           /
      !          -----
      !          j = 0

      nb_digits = get_nb_digits(fau_n(1), qs)
      allocate(ytemp(0:nb_digits-1))

      radinv = 0.0_qp
      temp = fau_n(1) ! TODO: decide on scalar or array for fau_n
      b_factor = 1.0_qp
      i = 0
      do
        b_factor = b_factor/qs
        ytemp(i) = modulo(temp, qs)
        radinv = radinv + ytemp(i)*b_factor
        temp = temp/qs  ! we rely on integer division here!
        if (temp==0) then
          exit
        end if
        i = i+1
      end do
      x(1) = radinv


      ! Find the other s-1 components.  Do this by using the
      ! recursion
      !                       k      {k - 1}
      !                      x  = C x
      !                       n      n
      ! instead of
      !
      !                       k    {k - 1}  {1}
      !                      x  = C        x
      !                       n             n
      !
      ! Notes:
      !   * we do transpose(COEF)*transpose([y_0, y_1, ..., y_{nb_digits-1}]) here!
      !   * due to the fact that we make use of this recursion, we cannot use different
      !     startindices for the points, all coordinates are based on the first one.
      !     
      ! 

      do k=2, s

        x(k) = 0.0_qp
        b_factor = inverse_qs

        do j=0, nb_digits-1

          !             nb_digits-1
          !                -----
          !                 \
          ! ytemp_new(j) =   )   modulo(ytemp_old(i)*binom(i, j), qs)
          !                 /
          !                -----
          !                i = j

          ztemp = 0

          do i=j, nb_digits-1
            ztemp = ztemp + coef(i,j)*ytemp(i)
          end do

          ytemp(j) = modulo(ztemp, qs)
          x(k) = x(k) + ytemp(j)*b_factor
          b_factor = b_factor*inverse_qs

        end do

      end do

      fau_n = fau_n + step

    end subroutine next_faure_bratley_fox


    subroutine next_faure_gray(x)

      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b) :: d
      integer(kind=i4b) :: j, k
      integer(kind=i4b) :: temp

      do d = 1, s

        ! Find the smallest index k for which the digit is different from
        ! qs-1 in the previous n.  We assume that indices start at 0.
        temp = fau_n(d)-1
        k = 0
        do
          if (modulo(temp, qs) /= qs-1) then
            exit
          end if
          k = k+1
          temp = temp/qs
        end do

        ! Add the k-th column of the generator matrix for this dimension to the
        ! vector with graycode-bits and calculate the radical inverse value.
        x(d) = 0.0_qp
        do j = max_digits, 1, -1
          graycode_bits(d, j) = modulo( graycode_bits(d,j) + c_eff(d,j-1,k), qs )
          x(d) = x(d)/qs + graycode_bits(d, j)
        end do
        x(d) = x(d)/qs

      end do

      fau_n = fau_n + step

    end subroutine next_faure_gray


  ! Construct the powers of the upper triangular Pascal matrix (modulo qs).
  !
  subroutine construct_pascal_power_matrices(p)
    integer(kind=i4b), dimension(0:,0:,0:), intent(inout) :: p

    integer(kind=i4b) :: row, col, mymax
    integer(kind=i4b) :: d, l, s
    integer(kind=i4b) :: temp

    s = size(p, dim=1)

    p = 0
    mymax = ubound(p, dim=2)

    ! Zero'th power is unity matrix.
    do row = 0, mymax
      p(0, row, row) = 1
    end do


    if (s > 1) then

      ! First power has first row and diagonal equal to 1...
      do col = 0, mymax
        p(1, 0, col) = 1
        p(1, col, col) = 1
      end do
      ! ... and the rest of the elements are found using the recursion.
      do row = 1, mymax
        do col = row+1, mymax
          p(1, row, col) = modulo( p(1, row-1, col-1) + p(1, row, col-1), qs )
        end do
      end do

    end if

    if (s > 2) then

      ! Now do the other powers of the Pascal matrix.
      do d = 2, ubound(p, dim=1)

        p(d, :, :) = 0
        do col = 0, mymax
          do row = 0, mymax
            if (col < row) then
              exit
            else
              temp = 0
              !temp = temp + dot_product(p(d-1, row, row:col), p(1, row:col, row))
              do l = 0, col
                temp = temp + p(d-1, row, l)*p(1, l, col)
              end do
            end if
            p(d, row, col) = modulo(temp, qs)
          end do
        end do

      end do

    end if

  end subroutine construct_pascal_power_matrices


  ! Setup the scrambling matrices for the specified type of scrambling.
  !
  subroutine generate_scrambling(scrambletype)
    character(len=*), intent(in) :: scrambletype

    integer(kind=i4b) :: i

    select case (scrambletype)

      case ("None")


      case ("GFaure")

        ! GFaure does not add a shift.
        allocate(gfaure_a(s, max_digits, max_digits))

        do i = 1,s
          call random_linear_scrambling(qs, gfaure_a(i,:,:))
        end do


      case ("RandomLinearScrambling")

        ! MATRIX + SHIFT
        allocate(owen_matrix(s, max_digits, max_digits), owen_shift(s, max_digits))

        do i = 1,s
          call random_linear_scrambling(qs, owen_matrix(i,:,:), owen_shift(i,:))
        end do


      case ("RandomLinearDigitScrambling", "InversivePre", "InversivePost")

        ! MATRIX + SHIFT
        allocate(owen_matrix(s, max_digits, max_digits), &
                 owen_shift(s, max_digits))

        do i = 1,s
          call random_linear_digit_scrambling(qs, owen_matrix(i,:,:), &
                                                  owen_shift(i,:))
        end do


      case ("FaureTezuka")
      
        ! TODO: see 'Faure, Tezuka, Another random scrambling of digital
        ! (t,s)-sequences' and also I-binomial paper page 748.
        ! Also called 'reordered Faure sequences'.


      case ("LeftIbinomial")

        ! MATRIX + SHIFT
        allocate(owen_matrix(s, max_digits, max_digits), &
                 owen_shift(s, max_digits))

        do i = 1,s
          call left_i_binomial_scrambling(qs, owen_matrix(i,:,:), &
                                              owen_shift(i,:))
        end do


      case ("RightIbinomial")


      case ("StripedMatrix")

        ! MATRIX + SHIFT
        allocate(owen_matrix(s, max_digits, max_digits), owen_shift(s, max_digits))

        do i = 1,s
          call striped_matrix_scrambling(qs, owen_matrix(i,:,:), owen_shift(i,:))
        end do

      case default

        write(unit=*, fmt="(A)") &
          "ERROR: Unknown scrambling method requested for Faure sequence!"
        stop

    end select

  end subroutine generate_scrambling


  ! Apply the scrambling.  For the Faure sequence, this means that we
  ! compute the effective generator matrices and shifts.
  !
  subroutine apply_scrambling(scrambletype)
    character(len=*), intent(in) :: scrambletype

    integer(kind=i4b) :: i

    select case (scrambletype)

      case ("None", "InversivePre", "InversivePost")

        c_eff = pascal_powers
        shift_eff = 0

        !do i = 1,s
        !  call print_matrix(c_eff(i,:,:))
        !end do

      case ("GFaure")

        do i = 1,s
          c_eff(i,:,:) = modulo(matmul(gfaure_a(i,:,:), &
                                       pascal_powers(i-1,:,:)), qs)
        end do
        shift_eff = 0


      case ("RandomLinearScrambling",      &
            "RandomLinearDigitScrambling", &
            "LeftIbinomial",               &
            "StripedMatrix")

        do i = 1,s
          c_eff(i,:,:) = modulo(matmul(owen_matrix(i,:,:), &
                                       pascal_powers(i-1,:,:)), qs)
          shift_eff(i,:) = owen_shift(i,:)
        end do


      case default

        write(unit=*, fmt="(A)") &
          "ERROR: applying an unknown scrambling method to the Faure sequence!"
        stop

    end select

  end subroutine apply_scrambling


    ! Free up memory that is consumed by certain datastructures
    ! related to the Faure sequence.
    !
    subroutine free_faure()
      
      if (allocated(coef)) then
        deallocate(coef)
      end if
      if (allocated(pascal_powers)) then
        deallocate(pascal_powers)
      end if
      if (allocated(graycode_bits)) then
        deallocate(graycode_bits)
      end if
      if (allocated(fau_n)) then
        deallocate(fau_n)
      end if
      if (allocated(step)) then
        deallocate(step)
      end if
      if (allocated(c_eff)) then
        deallocate(c_eff)
      end if
      if (allocated(shift_eff)) then
        deallocate(shift_eff)
      end if
      if (allocated(n_digits)) then
        deallocate(n_digits)
      end if
      if (allocated(temp_digits)) then
        deallocate(temp_digits)
      end if
      if (allocated(gfaure_a)) then
        deallocate(gfaure_a)
      end if
      if (allocated(owen_matrix)) then
        deallocate(owen_matrix)
      end if
      if (allocated(owen_shift)) then
        deallocate(owen_shift)
      end if

    end subroutine free_faure


    ! Print the given matrices.
    ! Only for debugging purposes.
    !
    subroutine print_matrices(c)

      integer(kind=i4b), dimension(:,:,:), intent(in) :: c

      integer(kind=i4b) :: d

      do d = 1, size(c, dim=1)
        write(unit=*, fmt="(A, I0.0, A)") "Matrix ", d, ":"
        call print_matrix(c(d,:,:))
      end do

    end subroutine print_matrices


    ! Print the given vectors.
    ! Only for debugging purposes.
    !
    subroutine print_vectors(v)

      integer(kind=i4b), dimension(:,:), intent(in) :: v

      integer(kind=i4b) :: d
      character(len=50) :: formatstring

      write(unit=formatstring, fmt="(A, I0.0, A)") "(A, I0.0, A, ", size(v,2), "I5)"

      do d = 1, size(v,1)
        write(unit=*, fmt=formatstring) "Vector ", d, ":", v(d,:)
      end do

    end subroutine print_vectors


end module mod_faure
