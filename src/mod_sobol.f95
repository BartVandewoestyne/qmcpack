! Module implementing the Sobol sequence.
!
! References:
!
!  [1] `On the systematic search in a hypercube', I. M. Sobol, SIAM
!      Journal of Numerical Analysis, Vol. 16, No. 5, October 1979,
!      page 790-793.
!
!  [2] `ALGORITHM 659: Implementing Sobol's Quasirandom Sequence
!      Generator', Bratley, Paul and Fox, L. Bennet, ACM Transactions
!      on Mathematical Software, Vol. 14, No. 1, March 1988, Pages 88-100.
!
!  [3] `Remark on Algorithm 659: Implementing Sobol's Quasirandom Sequence
!      Generator', Joe, Stephen and Kuo, Frances Y., ACM Transactions on
!      Mathematical Software, Vol. 29, No. 1, March 2003, Pages 49-57.
!
!  [4] Hee Sun Hong, Fred J. Hickernell, `Algorithm 823: Implementing
!      Scrambled Digital Sequences', ACM Transactions on Mathematical
!      Software, Vol. 29, No. 2, June 2003, Pages 95--109.
!
! TODO: 
!
!    * fix this little problem when generating sequences consisting
!      of just 1 dimension.  For now, sequences have to have at least
!      dimension 2.

module mod_sobol

  use numeric_kinds
  use mod_utilities
  use mod_radical_inverse
  use mod_primpoly

  private

  type, public :: soboltype
    character(len=20) :: poly_order
    character(len=20) :: dirnumbers
    logical :: use_ones
    logical :: use_antonov_saleev
  end type soboltype

  public :: init_sobol
  public :: next_sobol
  public :: sobol
  public :: free_sobol
  public :: print_direction_matrix
  public :: print_sobol_type


  private :: initialize_m_values
  private :: fill_direction_matrix
  private :: init_m_values_joe_kuo
  private :: init_m_values_bratley_fox
  private :: init_m_values_numerical_recipes

  integer(kind=i4b), private                              :: sob_s
  integer(kind=i4b), dimension(:), allocatable, private   :: sob_startindex
  integer(kind=i4b), dimension(:), allocatable, private   :: sob_n
  integer(kind=i4b), dimension(:), allocatable, private   :: sob_lastq
  integer(kind=i4b), private                              :: maxcol
  integer(kind=i4b), dimension(:,:), allocatable, private :: v
  real(kind=qp), private                                  :: recipd

contains


  subroutine init_sobol(maxpoints, s, startindex, mysoboltype)

    integer(kind=i4b), intent(in)               :: maxpoints
    integer(kind=i4b), intent(in)               :: s
    integer(kind=i4b), dimension(:), intent(in) :: startindex
    type(soboltype), intent(in)                 :: mysoboltype

    integer(kind=i4b)               :: i, j, temp
    integer(kind=i4b)               :: mydim
    integer(kind=i4b), dimension(s) :: n_gray

    if (any(startindex<1)) then
      write(unit=*, fmt="(A)") &
        "ERROR: array STARTINDEX for sobol-sequence should only contain "
      write(unit=*, fmt="(A)") &
        "       strict positive numbers!"
    end if

    sob_s = s
    if (allocated(sob_startindex)) then
      deallocate(sob_startindex)
    end if
    allocate(sob_startindex(s))
    sob_startindex = startindex

   ! Determine the number of bits in the largest number n we will use.
    maxcol = get_nb_digits(maxpoints+maxval(startindex), 2)

    ! Get the necessary initial m-values
    if (allocated(v)) then
      deallocate(v)
    end if
    allocate (v(maxcol, s))
    ! Values that are undefined yet are -1
    v = -1
    call initialize_m_values(maxcol, s, mysoboltype, v)

    ! Use the recurrence to fill up the rest of the direction matrix.
    call fill_direction_matrix(v, mysoboltype)

    ! Multiply the columns of v by the appropriate power of 2
    v = v*spread(2**(/ (i, i=maxcol-1,0,-1) /), 2, s)
    recipd = 1.0_qp/2**maxcol
    
    ! Set sob_lastq to the appropriate value belonging to graycode(startindex-1)
    ! so that the first point generated will be the correct value for
    ! graycode(startindex).
    if (allocated(sob_lastq)) then
      deallocate(sob_lastq)
    end if
    allocate(sob_lastq(sob_s))
    n_gray = graycode(sob_startindex-1)

    do mydim=1,sob_s

      temp = n_gray(mydim)
      j = 1
      sob_lastq(mydim) = modulo(temp, 2)*v(1, mydim)
      do
        temp = temp/2
        if (temp == 0) then
          exit
        end if
        j = j+1
        sob_lastq(mydim) = ieor(sob_lastq(mydim), modulo(temp,2)*v(j, mydim))
      end do

    end do

    ! At entry of next_sobol(), sob_n is the index of the previous point.
    if (allocated(sob_n)) then
      deallocate(sob_n)
    end if
    allocate(sob_n(sob_s))
    ! No points were generated yet...
    sob_n = sob_startindex-1

  end subroutine init_sobol


  ! Return the next Sobol point.
  !
  ! TODO:
  !   * Do not only support the Antonov-Saleev variant but also the slower one
  !     because it generates the points in another ordering.
  !
  subroutine next_sobol(x)

    real(kind=qp), dimension(:), intent(out) :: x

    integer(kind=i4b) :: mydim
    integer(kind=i4b) :: c

    do mydim=1,sob_s

        ! Search for the index of the rightmost zero bit of the previous
        ! (non-graycode) number.
        c = find_rightmost_zerobitpos(sob_n(mydim))

        sob_lastq(mydim) = ieor(sob_lastq(mydim), v(c, mydim))

    end do

    x = sob_lastq*recipd

    sob_n = sob_n+1

  end subroutine next_sobol


  subroutine sobol(n, s, startindex, x, mysoboltype)

    integer(kind=i4b), intent(in)                :: n, s
    integer(kind=i4b), dimension(:), intent(in)  :: startindex
    real(kind=qp), dimension(:,:), intent(out)   :: x
    type(soboltype), intent(in)                  :: mysoboltype

    integer(kind=i4b), dimension(n,s)            :: res
    !integer(kind=i4b)                 :: maxcol
    !integer(kind=i4b), allocatable    :: v(:,:)
    integer(kind=i4b), dimension(:), allocatable :: b
    integer(kind=i4b)                            :: i, j, w, mydim, c
    real(kind=qp)                                :: recipd


    ! Determine the number of bits in n
    maxcol = get_nb_digits(n+maxval(startindex), 2)

    ! Get the necessary initial m-values
    allocate (v(maxcol, s))
    ! Values that are undefined yet are -1
    v = -1
    call initialize_m_values(maxcol, s, mysoboltype, v)
    call fill_direction_matrix(v, mysoboltype)

    ! Multiply the columns of v by the appropriate power of 2
    v = v*spread(2**(/ (i, i=maxcol-1,0,-1) /), 2, s)
    recipd = 1.0_qp/2**maxcol

    ! TODO: optimize this for speed!  This is the bottleneck in our Sobol
    ! implementation
    if (.not. mysoboltype%use_antonov_saleev) then
      do mydim=1,s
         do j=startindex(mydim),startindex(mydim)+n-1
           allocate(b(get_nb_digits(j, 2)))
           call get_digits(j, 2, b)
           res(j-startindex(mydim)+1, mydim) = b(1)*v(1,mydim)
           do w=2,size(b)
             res(j-startindex(mydim)+1, mydim) = ieor(res(j-startindex(mydim)+1, mydim), b(w)*v(w,mydim))
           end do
           deallocate(b)
         end do
      end do
    else
      ! TODO: check with the plots of Numerical Recipes... they seem to
      ! differ a bit... but I'm not sure if those plots used the Antonov
      ! Saleev variant or not...
      do mydim=1,s
        ! Calculate first point using the traditional method:
        allocate(b(get_nb_digits(startindex(mydim), 2)))
        call get_digits(startindex(mydim), 2, b)
        res(1, mydim) = b(1)*v(1,mydim)
        do w=2,size(b)
          res(1, mydim) = ieor(res(1, mydim), b(w)*v(w,mydim))
        end do
        deallocate(b)
        ! Calculate the rest using the Antonov-Saleev method
        do j=2, n
          ! Find the position of the rightmost zero bit of the
          ! previous n
          c = find_rightmost_zerobitpos(startindex(mydim)+j-2)
          res(j, mydim) = ieor(res(j-1, mydim), v(c, mydim))
        end do
      end do
    end if

    x = res*recipd

    deallocate(v)

  end subroutine sobol


  ! Initialize the direction numbers.
  !
  ! TODO:
  !   * In the user manual of RandQMC, on page 6, another method to set up
  !     initial m-values is described...
  !
  subroutine initialize_m_values(w, s, mysoboltype, v)

    integer(kind=i4b), intent(in)                  :: w
    integer(kind=i4b), intent(in)                  :: s
    type(soboltype), intent(in)                    :: mysoboltype
    integer(kind=i4b), dimension(:,:), intent(out) :: v

    select case (mysoboltype%dirnumbers)

      case ("NumericalRecipes")
        call init_m_values_numerical_recipes(w, s, mysoboltype, v)

      case ("BratleyFox")
        call init_m_values_bratley_fox(w, s, mysoboltype, v)

      case ("JoeKuo")
        call init_m_values_joe_kuo(w, s, mysoboltype, v)

    end select

  end subroutine initialize_m_values


  ! Fill up the rest of the direction matrix.
  !
  subroutine fill_direction_matrix(v, mysoboltype)

    type(soboltype), intent(in)                      :: mysoboltype
    integer(kind=i4b), dimension(:,:), intent(inout) :: v

    integer(kind=i4b)                            :: skip
    integer(kind=i4b)                            :: maxcol, s
    integer(kind=i4b), dimension(:), allocatable :: p
    integer(kind=i4b), dimension(:), allocatable :: a
    integer(kind=i4b)                            :: i, j, k
    integer(kind=i4b)                            :: startdim, mydim
    integer(kind=i4b)                            :: newv, d

    maxcol = size(v, 1)
    s = size(v, 2)

    ! Check for what cases we have to skip the first dimension
    if ( ((mysoboltype%use_ones)                            &
          .OR. (mysoboltype%dirnumbers == "BratleyFox")     &
          .OR. (mysoboltype%dirnumbers == "JoeKuo")         &
          .OR. (mysoboltype%dirnumbers == "Jaeckel")        &
          .OR. (mysoboltype%dirnumbers == "Sobol76"))       &
         .AND. (mysoboltype%dirnumbers /= "PressTeukolsky") &
         .AND. (mysoboltype%dirnumbers /= "Winiarski")) then
      skip = 1
      allocate(p(s-1))
      call get_primpolys(s-1, p, mysoboltype%poly_order)
    else
      skip = 0
      allocate(p(s))
      call get_primpolys(s, p, mysoboltype%poly_order)
    end if

    ! Find first dimension that needs to be computed (where we still have
    ! a -1 value)
    !
    !   TODO:
    !     * check if there's a better way to do this...
    !     * we are in trouble here for 1-dimensinal sequences.  Check how and
    !       how to solve.
    startdim = s
    do j=1,s
      do i=1,maxcol
         if (v(i,j) == -1) then
           startdim = j
         end if
         if (startdim /= s) then
           exit
         end if
      end do
      if (startdim /= s) then
        exit
      end if
    end do


    do mydim = startdim,s

      ! Get the degree of the primpoly to see where we have to start filling up
      ! for this dimension
      !print *, "size(p) = ", size(p)
      !print *, "mydim-skip = ", mydim-skip
      d = get_degree(p(mydim-skip))

      ! Get the bit representation of the polynomial so we know what a_i
      ! we have to include in the recurrence relation
      allocate(a(d+1))
      call get_digits(p(mydim-skip), 2, a)
      
      ! Start filling up...
      do i=d+1,maxcol
        newv = v(i-d, mydim)
        do k=1,d
          if (a(d-k+1)==1) then
            newv = ieor(newv, 2**k*v(i-k, mydim))
          end if
        end do
        v(i, mydim) = newv
      end do

      deallocate(a)

    end do
    
    deallocate(p)

  end subroutine fill_direction_matrix


  ! Initialize with the direction numbers from Numerical Recipes.
  !
  subroutine init_m_values_numerical_recipes(w, s, mysoboltype, v)

    integer(kind=i4b), intent(in)                  :: w
    integer(kind=i4b), intent(in)                  :: s
    type(soboltype), intent(in)                    :: mysoboltype
    integer(kind=i4b), dimension(:,:), intent(out) :: v

    integer(kind=i4b), dimension(4,6) :: nrdata

    if (mysoboltype%use_ones) then
      stop "The Numerical Recipes implementation does *NOT* set the initial direction numbers to one!"
    end if

    nrdata = -1
    nrdata(1:1,1) = (/ 1 /)
    nrdata(1:2,2) = (/ 1, 1 /)
    nrdata(1:3,3) = (/ 1, 3, 7 /)
    nrdata(1:3,4) = (/ 1, 3, 3 /)
    nrdata(1:4,5) = (/ 1, 1, 3, 13 /)
    nrdata(1:4,6) = (/ 1, 1, 5,  9 /)

    ! Take only the part we need (e.g. if w or dim is smaller)
    v(1:min(4,w), 1:min(6,s)) = nrdata(1:min(4,w), 1:min(6,s))

  end subroutine init_m_values_numerical_recipes


  ! For Bratley and Fox, the initialization of this array is from the paper
  ! 'THE PRODUCTION OF POINTS UNIFORMLY DISTRIBUTED IN A MULTIDIMENSIONAL CUBE
  ! (IN RUSSIAN), PREPRINT IPM AKAD. NAUK SSSR, NO. 40, MOSCOW 1976.'
  !
  subroutine init_m_values_bratley_fox(w, s, mysoboltype, v)

    integer(kind=i4b), intent(in)                  :: w
    integer(kind=i4b), intent(in)                  :: s
    type(soboltype), intent(in)                    :: mysoboltype
    integer(kind=i4b), dimension(:,:), intent(out) :: v

    integer(kind=i4b), dimension(8,40) :: bfdata

    if  (.NOT. mysoboltype%use_ones) then
      print *, "The Bratley and Fox implementation *does* set the initial"
      print *, "direction numbers to all ones for the first dimension,"
      print *, "so you should do this too!"
      stop
    end if

    ! Initialize this to all -1 (we need this to determine where to fill up
    ! the direction numbers table)
    bfdata = -1

    ! Bratley and Fox initialize the first dimension to all ones and start
    ! using the recurrence relation from dimension 2 on.
    bfdata(1:8,1) = 1

    ! From the second dimension on, they use the recurrence relation, so
    ! now we need 'degree' initial values for each dimension.
    bfdata(1,2:40) = 1
    bfdata(2,3:40) = (/1,3,1,3,1,3,3,1,           &
                3,1,3,1,3,1,1,3,1,3,              &
                1,3,1,3,3,1,3,1,3,1,              &
                3,1,1,3,1,3,1,3,1,3/)
    bfdata(3,4:40) = (/7,5,1,3,3,7,5,             &
                5,7,7,1,3,3,7,5,1,1,              &
                5,3,3,1,7,5,1,3,3,7,              &
                5,1,1,5,7,7,5,1,3,3/)
    bfdata(4,6:40) = (/1,7,9,13,11,               &
                1,3,7,9,5,13,13,11,3,15,          &
                5,3,15,7,9,13,9,1,11,7,           &
                5,15,1,15,11,5,3,1,7,9/)
    bfdata(5,8:40) = (/9,3,27,                    &
                15,29,21,23,19,11,25,7,13,17,     &
                1,25,29,3,31,11,5,23,27,19,       &
                21,5,1,17,13,7,15,9,31,9/)
    bfdata(6,14:40) = (/37,33,7,5,11,39,63,       &
                27,17,15,23,29,3,21,13,31,25,     &
                9,49,33,19,29,11,19,27,15,25/)
    bfdata(7,20:40) = (/13,                       &
                33,115,41,79,17,29,119,75,73,105, &
                7,59,65,21,3,113,61,89,45,107/)
    bfdata(8,38:40) = (/7,23,39/)

    v(1:w, 1) = 1
    v(1:min(8,w), 1:min(40,s)) = bfdata(1:min(8,w), 1:min(40,s))
    
  end subroutine init_m_values_bratley_fox


  subroutine init_m_values_joe_kuo(w, s, mysoboltype, v)

    integer(kind=i4b), intent(in)                  :: w
    integer(kind=i4b), intent(in)                  :: s
    type(soboltype), intent(in)                    :: mysoboltype
    integer(kind=i4b), dimension(:,:), intent(out) :: v

    integer(kind=i4b), dimension(13,1111) :: jkdata

    if  (.NOT. mysoboltype%use_ones) then
      print *, "Joe and Kuos Sobol implementation *does* set the initial"
      print *, "direction numbers to all ones for the first dimension,"
      print *, "so you should do this too!"
      stop
    end if

    ! Initialize direction matrix to all -1 (we need this to determine where
    ! to start to fill up the direction numbers table)
    jkdata = -1

    ! Joe and Kuo initialize the first dimension to all ones and start
    ! using the recurrence relation from dimension 2 on.
    jkdata(1:13,1) = 1

   ! From the second dimension on, they use the recurrence relation, so
   ! now we need 'degree' initial values for each dimension.
    jkdata(1:1,2) = (/1/)
    jkdata(1:2,3) = (/1, 1/)
    jkdata(1:3,4) = (/1, 3, 7/)
    jkdata(1:3,5) = (/1, 1, 5/)
    jkdata(1:4,6) = (/1, 3, 1,  1/)
    jkdata(1:4,7) = (/1, 1, 3,  7/)
    jkdata(1:5,8) = (/1, 3, 3,  9,  9/)
    jkdata(1:5,9) = (/1, 3, 7, 13,  3/)
    jkdata(1:5,10) = (/1, 1, 5, 11, 27/)
    jkdata(1:5,11) = (/1, 3, 5,  1, 15/)
    jkdata(1:5,12) = (/1, 1, 7,  3, 29/)
    jkdata(1:5,13) = (/1, 3, 7,  7, 21/)
    jkdata(1:6,14) = (/1, 1, 1,  9, 23, 37/)
    jkdata(1:6,15) = (/1, 3, 3,  5, 19, 33/)
    jkdata(1:6,16) = (/1, 1, 3, 13, 11,  7/)
    jkdata(1:6,17) = (/1, 1, 7, 13, 25,  5/)
    jkdata(1:6,18) = (/1, 3, 5, 11,  7, 11/)
    jkdata(1:6,19) = (/1, 1, 1,  3, 13, 39/)
    jkdata(1:7,20) = (/1, 3, 1, 15, 17, 63,  13/)
    jkdata(1:7,21) = (/1, 1, 5,  5,  1, 59,  33/)
    jkdata(1:7,22) = (/1, 3, 3,  3, 25, 17, 115/)
    jkdata(1:7,23) = (/1, 1, 7, 15, 29, 15,  41/)
    jkdata(1:7,24) = (/1, 3, 1,  7,  3, 23,  79/)
    jkdata(1:7,25) = (/1, 3, 7,  9, 31, 29,  17/)
    jkdata(1:7,26) = (/1, 1, 5, 13, 11,  3,  29/)
    jkdata(1:7,27) = (/1, 1, 1,  9,  5, 21, 119/)
    jkdata(1:7,28) = (/1, 1, 3,  1, 23, 13,  75/)
    jkdata(1:7,29) = (/1, 3, 7, 11, 27, 31,  73/)
    jkdata(1:7,30) = (/1, 1, 7,  7, 19, 25, 105/)
    jkdata(1:7,31) = (/1, 3, 1,  5, 21,  9,   7/)
    jkdata(1:7,32) = (/1, 1, 1, 15,  5, 49,  59/)
    jkdata(1:7,33) = (/1, 3, 1,  1,  1, 33,  65/)
    jkdata(1:7,34) = (/1, 3, 5, 15, 17, 19,  21/)
    jkdata(1:7,35) = (/1, 1, 7, 11, 13, 29,   3/)
    jkdata(1:7,36) = (/1, 3, 7,  5,  7, 11, 113/)
    jkdata(1:7,37) = (/1, 1, 5, 11, 15, 19,  61/)
    jkdata(1:8,38) = (/1, 1, 1,  1,  9, 27,  89,   7/)
    jkdata(1:8,39) = (/1, 1, 3,  7, 31, 15,  45,  23/)
    jkdata(1:8,40) = (/1, 3, 3,  9, 25, 25, 107,  39/)
    jkdata(1:8,41) = (/1, 1, 7,  7,  3, 63,  21, 217/)
    jkdata(1:8,42) = (/1, 3, 5,  7,  5, 55,  71, 141/)
    jkdata(1:8,43) = (/1, 1, 5,  1, 23, 17,  79,  27/)
    jkdata(1:8,44) = (/1, 1, 5, 15,  7, 63,  19,  53/)
    jkdata(1:8,45) = (/1, 1, 3, 15,  3, 49,  71, 181/)
    jkdata(1:8,46) = (/1, 3, 3, 15, 17, 19,  61, 169/)
    jkdata(1:8,47) = (/1, 3, 3, 13, 23, 41,  41,  35/)
    jkdata(1:8,48) = (/1, 1, 1,  3,  3, 59,  57,  15/)
    jkdata(1:8,49) = (/1, 3, 1,  3,  3,  3, 121, 207/)
    jkdata(1:8,50) = (/1, 3, 5, 15, 21, 57,  87,  45/)
    jkdata(1:8,51) = (/1, 1, 1,  5, 25, 33, 119, 247/)
    jkdata(1:8,52) = (/1, 1, 1,  9, 25, 49,  55, 185/)
    jkdata(1:8,53) = (/1, 3, 5,  7, 23, 53,  85, 117/)
    jkdata(1:9,54) = (/1, 3, 3, 13, 11, 57, 121,  41, 235/)
    jkdata(1:9,55) = (/1, 1, 3,  3, 19, 57, 119,  81, 307/)
    jkdata(1:9,56) = (/1, 3, 3,  7,  3, 39,  11, 223, 495/)
    jkdata(1:9,57) = (/1, 3, 3,  5, 11, 21,  23, 151, 417/)
    jkdata(1:9,58) = (/1, 3, 1, 11, 31,  7,  61,  81,  57/)
    jkdata(1:9,59) = (/1, 1, 3,  9,  7, 53,  11, 189, 151/)
    jkdata(1:9,60) = (/1, 3, 7,  1,  9,  9,  35,  61,  19/)
    jkdata(1:9,61) = (/1, 1, 5,  9,  5, 55,  33,  95, 119/)
    jkdata(1:9,62) = (/1, 3, 7,  1, 17, 15,  43, 185, 375/)
    jkdata(1:9,63) = (/1, 1, 3,  5, 23, 59, 107,  23, 451/)
    jkdata(1:9,64) = (/1, 1, 7,  7, 17, 19, 113,  73,  55/)
    jkdata(1:9,65) = (/1, 3, 1, 13, 17, 49, 101, 113, 449/)
    jkdata(1:9,66) = (/1, 3, 3,  9, 25, 31,  29, 239, 501/)
    jkdata(1:9,67) = (/1, 1, 3,  9, 13,  3,  87,  85,  53/)
    jkdata(1:9,68) = (/1, 1, 5,  1, 11, 39, 119,   9, 185/)
    jkdata(1:9,69) = (/1, 1, 1,  7, 31,  5,  97, 201, 317/)
    jkdata(1:9,70) = (/1, 1, 3,  3, 27,  5,  29,  83,  17/)
    jkdata(1:9,71) = (/1, 3, 5,  5, 19, 41,  17,  53,  21/)
    jkdata(1:9,72) = (/1, 1, 5,  1, 17,  9,  89, 183, 487/)
    jkdata(1:9,73) = (/1, 1, 7, 11, 23, 19,   5, 203,  13/)
    jkdata(1:9,74) = (/1, 3, 7, 11,  7,  9, 127,  91, 347/)
    jkdata(1:9,75) = (/1, 1, 7, 13,  5, 57,  89, 149, 393/)
    jkdata(1:9,76) = (/1, 1, 1,  7, 11, 25, 119, 101,  15/)
    jkdata(1:9,77) = (/1, 1, 1,  7, 19,  1, 117,  13, 391/)
    jkdata(1:9,78) = (/1, 3, 3,  9, 19, 15, 103, 111, 307/)
    jkdata(1:9,79) = (/1, 3, 3,  9,  7, 51, 105, 239, 189/)
    jkdata(1:9,80) = (/1, 1, 1,  1, 13, 11,  41,   3, 381/)
    jkdata(1:9,81) = (/1, 3, 1,  1, 21, 19,  83, 205,  71/)
    jkdata(1:9,82) = (/1, 3, 5,  3, 21, 61,  25, 253, 163/)
    jkdata(1:9,83) = (/1, 1, 1,  9,  7, 53,  41, 247,  99/)
    jkdata(1:9,84) = (/1, 3, 5, 15,  9, 29,  55, 121, 467/)
    jkdata(1:9,85) = (/1, 3, 7,  1, 11, 19,  69, 189, 167/)
    jkdata(1:9,86) = (/1, 3, 5,  5,  1, 11, 117, 169, 433/)
    jkdata(1:9,87) = (/1, 1, 1, 13,  5,  9,  49, 179, 337/)
    jkdata(1:9,88) = (/1, 3, 7,  1, 21, 21, 127, 197, 257/)
    jkdata(1:9,89) = (/1, 3, 5,  9, 11, 19,  29, 175, 179/)
    jkdata(1:9,90) = (/1, 3, 3,  9, 13, 43,   1, 217,  47/)
    jkdata(1:9,91) = (/1, 1, 3,  9, 25, 13,  99, 249, 385/)
    jkdata(1:9,92) = (/1, 3, 1,  9,  9, 13,  53, 195,  23/)
    jkdata(1:9,93) = (/1, 3, 5,  9,  7, 41,  83,  95, 117/)
    jkdata(1:9,94) = (/1, 1, 7, 13,  7, 25,  15,  63, 369/)
    jkdata(1:9,95) = (/1, 3, 1, 11, 27, 31,  31,  19, 425/)
    jkdata(1:9,96) = (/1, 3, 7,  3, 15,  9,  73,   7, 207/)
    jkdata(1:9,97) = (/1, 3, 5,  5, 25, 11, 115,   5, 433/)
    jkdata(1:9,98) = (/1, 1, 1, 11, 15, 19,  35,  75, 301/)
    jkdata(1:9,99) = (/1, 3, 7, 11, 21,  5,  21, 217, 147/)
    jkdata(1:9,100) = (/1, 1, 3, 13, 17, 53,  89, 245, 333/)
    jkdata(1:9,101) = (/1, 3, 1,  5, 19, 37,   5, 111,  85/)
    jkdata(1:10,102) = (/1, 1, 7,  3, 19,  7,   1, 189, 221,  519/)
    jkdata(1:10,103) = (/1, 1, 1, 15, 21, 51,  91, 165, 423,  307/)
    jkdata(1:10,104) = (/1, 3, 7,  1,  5, 45,  53, 169,  49,  931/)
    jkdata(1:10,105) = (/1, 3, 3, 11, 11,  7,  35, 141,   3, 1023/)
    jkdata(1:10,106) = (/1, 1, 3, 11,  3,  7,  95, 221,  43,  517/)
    jkdata(1:10,107) = (/1, 3, 5,  7,  5, 61,  83, 249, 229,  771/)
    jkdata(1:10,108) = (/1, 3, 7, 13, 29, 23,  19, 159, 227,  151/)
    jkdata(1:10,109) = (/1, 1, 3, 15, 31, 45,  85, 253, 201, 1023/)
    jkdata(1:10,110) = (/1, 1, 3, 11, 29,  7,  55, 207, 383,  539/)
    jkdata(1:10,111) = (/1, 1, 5, 13,  5, 59,  51, 249, 281,  725/)
    jkdata(1:10,112) = (/1, 3, 1,  9,  5, 41, 101, 219, 229,   45/)
    jkdata(1:10,113) = (/1, 3, 3, 11,  1,  1,  33,  23, 207,  927/)
    jkdata(1:10,114) = (/1, 1, 3, 15, 31, 29,  41,  49,  21,  707/)
    jkdata(1:10,115) = (/1, 3, 1, 15, 27, 61,  55, 127, 343,   29/)
    jkdata(1:10,116) = (/1, 3, 3, 13, 11, 37,  45, 237, 251,  125/)
    jkdata(1:10,117) = (/1, 1, 5,  3, 13, 27,  95,   5, 397,  371/)
    jkdata(1:10,118) = (/1, 3, 1, 15,  1, 47,  61,  25, 173,  275/)
    jkdata(1:10,119) = (/1, 1, 3,  7,  3, 15,  27, 177, 507,  279/)
    jkdata(1:10,120) = (/1, 1, 3,  9,  7, 31,  37,  37, 421,  817/)
    jkdata(1:10,121) = (/1, 3, 3, 11, 11, 35,  89, 103, 443,  389/)
    jkdata(1:10,122) = (/1, 3, 7, 13,  7, 31,  75,  65, 399,  453/)
    jkdata(1:10,123) = (/1, 3, 1, 11,  3, 17,  57, 167,  53,  989/)
    jkdata(1:10,124) = (/1, 1, 1,  9, 23, 51,  61,  81, 345, 1015/)
    jkdata(1:10,125) = (/1, 1, 7,  9, 13, 13,  15,  87,  77,   29/)
    jkdata(1:10,126) = (/1, 1, 3,  5, 31, 25, 117, 119, 385,  169/)
    jkdata(1:10,127) = (/1, 3, 1, 13, 17, 45,  15,  45, 317,  743/)
    jkdata(1:10,128) = (/1, 1, 3,  9,  1,  5,  21,  79, 155,   99/)
    jkdata(1:10,129) = (/1, 1, 7,  1, 27,  5,  27, 143, 187,  923/)
    jkdata(1:10,130) = (/1, 3, 5, 13, 11, 33,  25,  57, 269,  981/)
    jkdata(1:10,131) = (/1, 1, 5,  7, 25, 39,  27,  79, 501,  181/)
    jkdata(1:10,132) = (/1, 1, 7,  7,  1,  5, 123, 187,  19,  693/)
    jkdata(1:10,133) = (/1, 3, 5,  7, 23, 47,  39, 143, 169,  309/)
    jkdata(1:10,134) = (/1, 3, 5,  7, 29, 29, 109, 183, 235,  227/)
    jkdata(1:10,135) = (/1, 1, 3,  7, 17, 35,  93,  75, 415,  111/)
    jkdata(1:10,136) = (/1, 3, 1,  5, 25, 47,  51,  97,  61,  219/)
    jkdata(1:10,137) = (/1, 1, 3,  9,  7, 63,  21, 211, 247,  897/)
    jkdata(1:10,138) = (/1, 3, 3,  7, 25, 45,  91, 149, 183,  377/)
    jkdata(1:10,139) = (/1, 3, 3, 13, 27, 37, 109, 175,   5,  425/)
    jkdata(1:10,140) = (/1, 3, 1, 11, 17, 47, 107,  37, 257,  609/)
    jkdata(1:10,141) = (/1, 3, 3,  9, 13, 59,  45, 135, 401,  227/)
    jkdata(1:10,142) = (/1, 1, 3, 11, 17, 21,  15, 189, 451,   19/)
    jkdata(1:10,143) = (/1, 1, 7, 15, 23, 59,  93, 225,  95,  221/)
    jkdata(1:10,144) = (/1, 1, 3,  3,  5, 33, 127, 241, 455,  143/)
    jkdata(1:10,145) = (/1, 3, 3, 13, 17, 51,   3,  63,  49,  581/)
    jkdata(1:10,146) = (/1, 3, 1, 11,  5,  9,  53,  33, 489,  147/)
    jkdata(1:10,147) = (/1, 1, 7,  1, 13, 27,  81,  43,  75,  919/)
    jkdata(1:10,148) = (/1, 1, 5, 11, 11, 13,  79,  13, 459,  127/)
    jkdata(1:10,149) = (/1, 3, 1,  3, 21, 25, 107,  73, 377,  725/)
    jkdata(1:10,150) = (/1, 1, 7,  3,  5, 43,  79, 213,  87,  793/)
    jkdata(1:10,151) = (/1, 1, 7,  9, 11,  3,  87,  57, 463,  289/)
    jkdata(1:10,152) = (/1, 1, 5, 11,  5, 17,  35, 239, 155,  411/)
    jkdata(1:10,153) = (/1, 1, 7,  1,  9, 21, 109, 183, 233,  835/)
    jkdata(1:10,154) = (/1, 1, 5,  7, 31, 59,  73, 117, 115,  921/)
    jkdata(1:10,155) = (/1, 1, 1,  1, 19, 61,  35,  21, 429,  957/)
    jkdata(1:10,156) = (/1, 3, 3, 15, 17, 27,  83,  29, 211,  443/)
    jkdata(1:10,157) = (/1, 1, 1, 15,  9, 47, 107, 115, 419,  349/)
    jkdata(1:10,158) = (/1, 3, 7,  3,  9, 57,   1,  43, 143,  813/)
    jkdata(1:10,159) = (/1, 1, 3,  1, 27, 11,  51, 205, 487,    5/)
    jkdata(1:10,160) = (/1, 1, 7,  9, 21, 17,   7, 223, 195,  105/)
    jkdata(1:10,161) = (/1, 1, 3,  1, 15, 39,  59,  15, 209,  457/)
    jkdata(1:11,162) = (/1, 3, 5,  7, 15,  1,  33,   3, 461,  393,    7/)
    jkdata(1:11,163) = (/1, 1, 7, 13,  1, 63, 115, 159, 193,  539, 2011/)
    jkdata(1:11,164) = (/1, 3, 3, 11,  1, 21,  43,  51, 157,  101, 1001/)
    jkdata(1:11,165) = (/1, 1, 1,  3, 29, 59, 111, 101, 193,  197,   49/)
    jkdata(1:11,166) = (/1, 3, 3, 13,  5, 17,  45, 127, 363,  697,  825/)
    jkdata(1:11,167) = (/1, 3, 3, 11, 31, 13, 121,  99, 181,   27,  415/)
    jkdata(1:11,168) = (/1, 3, 3,  7, 11, 31, 105, 239, 271,  343, 1441/)
    jkdata(1:11,169) = (/1, 1, 1,  3, 17,  3, 125, 171, 445,  515,  383/)
    jkdata(1:11,170) = (/1, 1, 5,  3, 23, 31,  87, 113, 381,   69, 1581/)
    jkdata(1:11,171) = (/1, 3, 7,  5, 19,  7, 101, 171, 231,  485,  623/)
    jkdata(1:11,172) = (/1, 3, 3, 13, 21,  9,  41, 119, 135,  383, 1621/)
    jkdata(1:11,173) = (/1, 1, 3, 11, 25, 27,  95, 189, 327,  855, 1319/)
    jkdata(1:11,174) = (/1, 3, 7,  5, 15, 37,  75, 245, 403,  693, 1387/)
    jkdata(1:11,175) = (/1, 1, 7, 11, 11, 23,   1, 201, 171,  133,  619/)
    jkdata(1:11,176) = (/1, 3, 7,  1,  5, 31,  57,  27, 197,   87,  839/)
    jkdata(1:11,177) = (/1, 1, 5,  3,  5,  9, 117, 185, 181,  743,  217/)
    jkdata(1:11,178) = (/1, 1, 3,  9,  1, 45,  21, 229, 343,  747,   75/)
    jkdata(1:11,179) = (/1, 3, 1,  7, 19, 43,  27, 105, 113,  475, 1955/)
    jkdata(1:11,180) = (/1, 1, 7, 15, 19, 31,  67, 153, 313,   87,  505/)
    jkdata(1:11,181) = (/1, 3, 1,  7, 19, 63,  29, 189, 393,  469,  281/)
    jkdata(1:11,182) = (/1, 1, 3,  5,  7, 21,  53,  33, 311,  763, 1629/)
    jkdata(1:11,183) = (/1, 3, 7, 13, 13, 39, 117,  35, 415,  721, 1379/)
    jkdata(1:11,184) = (/1, 1, 5,  7, 21, 51,  63, 137, 267,  345,   53/)
    jkdata(1:11,185) = (/1, 3, 3,  9, 17, 27,   1,  77, 247,  479, 1111/)
    jkdata(1:11,186) = (/1, 1, 3, 13, 17,  7,  77,  97, 425,  965, 1399/)
    jkdata(1:11,187) = (/1, 1, 3, 15, 25, 53,  89,  17, 233,  527,  301/)
    jkdata(1:11,188) = (/1, 1, 7, 13, 23, 11, 115, 181, 289,  121,  209/)
    jkdata(1:11,189) = (/1, 3, 1,  9, 19,  1,  49,  55,  55,  271,   49/)
    jkdata(1:11,190) = (/1, 3, 1,  7, 23, 59, 127, 197,  39,  353,  155/)
    jkdata(1:11,191) = (/1, 1, 3, 15, 15, 39,  15, 201, 247,  467, 1647/)
    jkdata(1:11,192) = (/1, 3, 1,  7, 13, 23,  79, 155, 327,  177,  631/)
    jkdata(1:11,193) = (/1, 3, 5,  9,  5, 49,  81,  37, 141,  245,  129/)
    jkdata(1:11,194) = (/1, 1, 7,  5, 19, 23,  29, 197,   5,  627, 1569/)
    jkdata(1:11,195) = (/1, 3, 1, 11, 25,  7,  65, 137, 189,  113,  335/)
    jkdata(1:11,196) = (/1, 1, 3, 11,  9, 55, 103, 223, 183,  357,   67/)
    jkdata(1:11,197) = (/1, 1, 5, 13,  7, 59,  33,  25,  27,    7, 1955/)
    jkdata(1:11,198) = (/1, 1, 3, 13,  3,  3,  73, 179, 337,  691, 1611/)
    jkdata(1:11,199) = (/1, 3, 5,  9, 21, 19,  79,  91, 341,  725, 2021/)
    jkdata(1:11,200) = (/1, 1, 3,  3, 17, 35,  29,  23, 327,  355, 1305/)
    jkdata(1:11,201) = (/1, 3, 3,  5, 25, 13,  21, 235,  87,  889,  121/)
    jkdata(1:11,202) = (/1, 1, 7, 13,  1,  9, 113,  53, 429,  635,   37/)
    jkdata(1:11,203) = (/1, 1, 5,  9, 27, 13,  31, 253, 357,  737,  877/)
    jkdata(1:11,204) = (/1, 3, 5, 11, 25, 15,  33,  49, 265,  429,  835/)
    jkdata(1:11,205) = (/1, 1, 3, 15, 27, 23, 107, 181, 251,  545, 1457/)
    jkdata(1:11,206) = (/1, 1, 3, 11, 25,  9,  95, 249, 437,  925,  669/)
    jkdata(1:11,207) = (/1, 3, 1,  7,  9,  7, 111,  53, 201,  357, 1405/)
    jkdata(1:11,208) = (/1, 3, 3,  1, 13, 43,  59, 173,  29,  873,  935/)
    jkdata(1:11,209) = (/1, 1, 7,  7,  3, 55,  99,  97, 339,  187, 1735/)
    jkdata(1:11,210) = (/1, 1, 7, 13, 17,  3, 117, 247, 257,  351,  665/)
    jkdata(1:11,211) = (/1, 3, 7,  3, 25, 19,  63,  67, 377,  677,  551/)
    jkdata(1:11,212) = (/1, 3, 1, 13, 23,  9,  63, 115,  17,  999,  789/)
    jkdata(1:11,213) = (/1, 3, 5,  3,  9, 27,  99, 103,  53,  921, 1543/)
    jkdata(1:11,214) = (/1, 1, 7, 13, 25, 33,  39, 159, 327,  477, 1267/)
    jkdata(1:11,215) = (/1, 3, 1,  9,  9, 27,   9, 239,  47,  233, 1027/)
    jkdata(1:11,216) = (/1, 3, 3, 15, 13, 49,  35,  69, 375,  765,    1/)
    jkdata(1:11,217) = (/1, 3, 1,  7, 17, 23,  63, 173, 393,  495, 1911/)
    jkdata(1:11,218) = (/1, 1, 1, 13, 17, 47, 125, 217, 369,   81,  163/)
    jkdata(1:11,219) = (/1, 3, 7, 13,  3, 19,  99,  95, 403,  953, 1929/)
    jkdata(1:11,220) = (/1, 1, 1,  3, 15,  7,  45, 221, 125,  479,   67/)
    jkdata(1:11,221) = (/1, 3, 3, 13,  7, 11,  93, 247, 429,   89, 1975/)
    jkdata(1:11,222) = (/1, 1, 1, 15,  7, 55,  33,  97, 257,  173, 1681/)
    jkdata(1:11,223) = (/1, 1, 7, 15, 29, 27,  93,  91, 157,  473, 1413/)
    jkdata(1:11,224) = (/1, 1, 1, 11,  3, 35,   9, 123, 217,  131,  191/)
    jkdata(1:11,225) = (/1, 3, 5,  9, 19,  5, 105, 223,  85,  961, 1711/)
    jkdata(1:11,226) = (/1, 1, 3, 13, 29,  5,  75, 213, 267,  411, 1307/)
    jkdata(1:11,227) = (/1, 1, 5,  9, 29, 55,  51, 129, 117,  291,  401/)
    jkdata(1:11,228) = (/1, 1, 3, 15, 19, 35, 115, 181, 337,  967,  725/)
    jkdata(1:11,229) = (/1, 3, 1,  1, 29, 37,  11,  87, 447,   65, 1229/)
    jkdata(1:11,230) = (/1, 1, 1,  1, 13,  9,  37, 239, 219,  511, 1403/)
    jkdata(1:11,231) = (/1, 1, 5, 15, 15, 33,  17,  85, 501,   13, 1609/)
    jkdata(1:11,232) = (/1, 1, 5, 11, 25, 29,  41,  89,  41,  805, 2035/)
    jkdata(1:11,233) = (/1, 1, 3, 11, 27, 47,  21, 249,  41,  945,  917/)
    jkdata(1:11,234) = (/1, 1, 3,  7,  1, 25,  43, 141, 193,  369,  921/)
    jkdata(1:11,235) = (/1, 3, 5,  1,  3, 11,  73,  39, 509,  827, 1789/)
    jkdata(1:11,236) = (/1, 3, 7, 11,  9, 47,  19,  57, 131,  295,   41/)
    jkdata(1:11,237) = (/1, 3, 1, 13,  9, 53,  93, 249, 207,  163, 2003/)
    jkdata(1:11,238) = (/1, 1, 5,  9, 13, 61,   7,  71, 505,  835,  187/)
    jkdata(1:11,239) = (/1, 1, 3, 13, 31, 59,  95, 101, 421,  259,   67/)
    jkdata(1:11,240) = (/1, 1, 7,  3, 29,  3,  81, 159, 149,  207, 1635/)
    jkdata(1:11,241) = (/1, 1, 7,  5, 31, 53,  93,  33, 111,  331,  717/)
    jkdata(1:11,242) = (/1, 3, 3, 11,  5, 47,  79, 137, 177,   29, 1449/)
    jkdata(1:11,243) = (/1, 3, 5, 13, 15,  5,  81, 189, 167,  315,  277/)
    jkdata(1:11,244) = (/1, 3, 3,  9, 29, 19,  55,  71, 223,  999, 1903/)
    jkdata(1:11,245) = (/1, 1, 3,  9,  1, 59,   9, 253, 291,  133, 1179/)
    jkdata(1:11,246) = (/1, 3, 1, 13, 19,  5,  51, 205,  91,  967,  363/)
    jkdata(1:11,247) = (/1, 3, 7,  1,  5, 47,  63, 171,  29,   41, 1211/)
    jkdata(1:11,248) = (/1, 1, 3, 11,  9, 23,  45,  13, 305,  117, 1231/)
    jkdata(1:11,249) = (/1, 1, 1, 15, 19, 45,  89, 249, 151,  677,  647/)
    jkdata(1:11,250) = (/1, 1, 3, 13,  5, 53,  73, 109, 177,  471, 1261/)
    jkdata(1:11,251) = (/1, 1, 5,  3, 15,  3,  19, 131, 337,  717, 1029/)
    jkdata(1:11,252) = (/1, 3, 7, 13,  3, 49, 115, 199, 183,  881, 1485/)
    jkdata(1:11,253) = (/1, 1, 1,  7,  5, 61,  39, 189, 361,  755, 1309/)
    jkdata(1:11,254) = (/1, 1, 3, 15,  7, 47,  47, 179, 435,  351, 1149/)
    jkdata(1:11,255) = (/1, 3, 7,  1, 15, 39,  81,  31, 307,  723,  317/)
    jkdata(1:11,256) = (/1, 1, 1, 15, 17, 29,  39,  99, 507,  259, 1335/)
    jkdata(1:11,257) = (/1, 3, 5,  3, 17, 17,   5, 113,  77,  879,  171/)
    jkdata(1:11,258) = (/1, 3, 1,  3, 23, 57,   5,  41, 181,  455,  243/)
    jkdata(1:11,259) = (/1, 1, 3, 11, 11,  5,  45, 173, 507,  721,  271/)
    jkdata(1:11,260) = (/1, 1, 1,  7,  9, 17,  53,  23, 315,  289, 1055/)
    jkdata(1:11,261) = (/1, 3, 5, 13, 23, 31,  65, 189, 145,  149, 1601/)
    jkdata(1:11,262) = (/1, 3, 3,  7, 19, 23,  49, 197, 423,  199, 1129/)
    jkdata(1:11,263) = (/1, 1, 1,  7,  3, 41,  17,   3,  71,  805, 1653/)
    jkdata(1:11,264) = (/1, 1, 7,  9, 17, 39, 105, 135, 103,  987,  205/)
    jkdata(1:11,265) = (/1, 1, 1,  7,  1,  5,  13,   9, 493,  851, 1463/)
    jkdata(1:11,266) = (/1, 1, 5,  5, 27, 27, 107,  95, 271,  423, 1681/)
    jkdata(1:11,267) = (/1, 3, 5, 15,  9,  7,   5, 195, 469,  597, 1621/)
    jkdata(1:11,268) = (/1, 1, 5,  9,  9, 29,   5,  27, 339,  129,  197/)
    jkdata(1:11,269) = (/1, 3, 3,  5, 17, 29,  19, 183, 237,   11,  951/)
    jkdata(1:11,270) = (/1, 3, 7,  5, 13, 33,  73,   1, 437,  733,  573/)
    jkdata(1:11,271) = (/1, 1, 1,  7, 25, 31,  59, 123, 483,  549, 1697/)
    jkdata(1:11,272) = (/1, 3, 1, 15, 29, 41,  43,  73,  31,  153, 1265/)
    jkdata(1:11,273) = (/1, 3, 7, 13, 23, 31,  83,  53, 219,  285, 1321/)
    jkdata(1:11,274) = (/1, 1, 3, 15, 29, 29,  97,  99,  61,  451, 1805/)
    jkdata(1:11,275) = (/1, 1, 1,  5, 11, 17, 115, 197, 131,  559, 1235/)
    jkdata(1:11,276) = (/1, 1, 1, 15, 31, 29,  27,  59, 391,  377, 1853/)
    jkdata(1:11,277) = (/1, 3, 7,  5, 25, 29,   1,  27, 233,  109, 1307/)
    jkdata(1:11,278) = (/1, 3, 5,  3, 21,  9,  69, 101, 219,  357,  945/)
    jkdata(1:11,279) = (/1, 3, 7,  1, 29,  9, 103,  55,  69,  143, 1197/)
    jkdata(1:11,280) = (/1, 1, 5, 11, 19, 31,   3, 193,  57,  693, 1411/)
    jkdata(1:11,281) = (/1, 3, 7,  7, 27, 27,  99,  31, 459,  615,  833/)
    jkdata(1:11,282) = (/1, 3, 7,  1, 31, 53, 103,  61, 225,  677,  273/)
    jkdata(1:11,283) = (/1, 1, 3,  5,  3, 35,  63, 119, 421,  701, 1517/)
    jkdata(1:11,284) = (/1, 3, 7,  7,  5,  5,  67,  11,   7,  475, 1747/)
    jkdata(1:11,285) = (/1, 3, 1,  9,  3, 61,  25,   7, 461,  767, 1095/)
    jkdata(1:11,286) = (/1, 1, 3,  3,  3,  1, 121, 255, 111,   85, 1345/)
    jkdata(1:11,287) = (/1, 3, 7, 11, 13, 49,  97, 233, 451,  229,  869/)
    jkdata(1:11,288) = (/1, 1, 7,  1, 21, 13,  77,  53, 277,  509,   57/)
    jkdata(1:11,289) = (/1, 3, 3, 15,  9, 57,  13, 157, 185,  547, 1383/)
    jkdata(1:11,290) = (/1, 3, 5,  1, 29, 29,  83, 193, 193,  151,  221/)
    jkdata(1:11,291) = (/1, 3, 1,  3,  3,  5, 103,  97, 125,  389, 1713/)
    jkdata(1:11,292) = (/1, 1, 1, 15, 17, 21,  41,  83, 251,  711,  335/)
    jkdata(1:11,293) = (/1, 3, 7, 11, 11, 43,  11,  65, 199,  785, 1751/)
    jkdata(1:11,294) = (/1, 1, 1, 13, 11, 25,  27,  81,  73,  657, 1141/)
    jkdata(1:11,295) = (/1, 1, 5,  5,  9, 57,  81, 239,  71,  319,  839/)
    jkdata(1:11,296) = (/1, 3, 5, 13, 21, 49,  37, 167,   7,  509,  523/)
    jkdata(1:11,297) = (/1, 1, 5,  1, 19, 37,  33,  69, 409,   99, 1861/)
    jkdata(1:11,298) = (/1, 3, 1,  7,  7, 27, 125,  71, 417, 1007, 1105/)
    jkdata(1:11,299) = (/1, 1, 5,  1, 17, 11,  71, 109, 149,  775,  389/)
    jkdata(1:11,300) = (/1, 1, 1, 15, 31, 61,  41,  97, 193,  359, 1177/)
    jkdata(1:11,301) = (/1, 1, 7,  7, 25, 37,  41, 137,  53,  697, 1877/)
    jkdata(1:11,302) = (/1, 3, 5,  5,  1, 49,  59,  71, 437,  677,  805/)
    jkdata(1:11,303) = (/1, 3, 5,  1, 27,  5,  41, 193,  29,   85,   93/)
    jkdata(1:11,304) = (/1, 3, 7,  1,  5, 63,  87, 189, 467,  497, 1591/)
    jkdata(1:11,305) = (/1, 1, 1, 15, 15, 63, 123, 115, 229,  105,  423/)
    jkdata(1:11,306) = (/1, 1, 1, 13, 27,  3,  43,  79,  31,  615, 1835/)
    jkdata(1:11,307) = (/1, 3, 7, 11, 29, 45, 101, 205,  35,  891,   99/)
    jkdata(1:11,308) = (/1, 1, 1, 11, 29, 37,  63,  37,  75,   71, 1781/)
    jkdata(1:11,309) = (/1, 3, 7, 13, 29, 63,  45, 227, 105,  449, 1515/)
    jkdata(1:11,310) = (/1, 1, 7,  5, 25, 21,  39,  53, 503,  835, 1909/)
    jkdata(1:11,311) = (/1, 1, 1, 11, 27, 21,  21,  33,  75,  609, 1011/)
    jkdata(1:11,312) = (/1, 1, 1,  7, 25, 19,  97,  91, 317,  377,  303/)
    jkdata(1:11,313) = (/1, 1, 3,  9,  3, 27,  15, 229, 401,  693,  385/)
    jkdata(1:11,314) = (/1, 1, 3,  7, 21, 59,  97, 245, 367,  665, 1635/)
    jkdata(1:11,315) = (/1, 1, 3,  1, 17, 21, 111, 105, 131,  627,  357/)
    jkdata(1:11,316) = (/1, 3, 7,  5, 25, 45,  21,  77, 365,  215,  973/)
    jkdata(1:11,317) = (/1, 1, 7,  3, 13, 23,  49, 229, 441,  911, 1781/)
    jkdata(1:11,318) = (/1, 1, 5,  9, 15, 13,  13, 161, 433,  503, 1707/)
    jkdata(1:11,319) = (/1, 3, 3,  5, 17, 15,  17, 103,  93,  729, 1363/)
    jkdata(1:11,320) = (/1, 1, 7,  5, 13,  3,  79,  93, 377,  131, 1053/)
    jkdata(1:11,321) = (/1, 3, 3, 11, 23, 43,  91,  13, 405,   19,  649/)
    jkdata(1:11,322) = (/1, 3, 1,  5,  9, 63,  65, 161, 465,  895, 1469/)
    jkdata(1:11,323) = (/1, 1, 3,  1,  3, 39, 105, 229, 259,  199,  623/)
    jkdata(1:11,324) = (/1, 1, 7,  7, 11, 19,  75, 223, 283,  161, 1429/)
    jkdata(1:11,325) = (/1, 1, 5,  1,  7, 63,   1,  69, 443,  239, 1241/)
    jkdata(1:11,326) = (/1, 1, 3, 11,  9, 31,  45,  15, 143,  633, 1151/)
    jkdata(1:11,327) = (/1, 3, 3,  7,  9, 41,  67,  25, 445, 1013, 1055/)
    jkdata(1:11,328) = (/1, 1, 5,  9,  7, 41,  83,  23,   3,  537,  503/)
    jkdata(1:11,329) = (/1, 3, 7, 13, 17, 15, 107, 233, 461,  255,  921/)
    jkdata(1:11,330) = (/1, 1, 1, 15,  7, 43, 125,  93, 329,   23,    3/)
    jkdata(1:11,331) = (/1, 3, 1, 13,  1, 63,  87,  25, 309,  149,  349/)
    jkdata(1:11,332) = (/1, 1, 5,  3, 27, 53,  15, 217,  77,  679, 1149/)
    jkdata(1:11,333) = (/1, 1, 5,  1,  1,  1,  81, 247, 323, 1021,  293/)
    jkdata(1:11,334) = (/1, 1, 7, 11,  9, 63,  95,  61, 155,  595,   45/)
    jkdata(1:11,335) = (/1, 1, 7, 13,  5, 31, 105,  75, 347,  199,  303/)
    jkdata(1:11,336) = (/1, 3, 1, 15, 31,  7,  65,  27,  45,  557,  877/)
    jkdata(1:11,337) = (/1, 3, 1,  1, 21, 17,  45,   9, 381,  659, 1565/)
    jkdata(1:12,338) = (/1, 1, 1,  1, 25, 11,  59, 223, 315,  251, 1583, 3915/)
    jkdata(1:12,339) = (/1, 1, 1, 11, 25, 61, 103, 213, 463,  829, 1001,   97/)
    jkdata(1:12,340) = (/1, 1, 5,  9, 21, 31,  23,  55, 207,  727,  663, 3047/)
    jkdata(1:12,341) = (/1, 1, 5, 13, 11, 51, 103, 197, 321,  439, 1535,  937/)
    jkdata(1:12,342) = (/1, 1, 5,  3,  1, 37,  99, 145, 157,  495,  395, 2897/)
    jkdata(1:12,343) = (/1, 3, 7, 13, 23, 29,  67,  89, 109,  647, 1141,  953/)
    jkdata(1:12,344) = (/1, 3, 5, 11, 19, 59,  99, 199, 479,  223, 1481,  127/)
    jkdata(1:12,345) = (/1, 3, 7, 15, 27, 25,  47,  41, 313,  949, 1797, 1201/)
    jkdata(1:12,346) = (/1, 1, 1, 13, 15, 63, 117, 201, 345,  625,  643, 3819/)
    jkdata(1:12,347) = (/1, 1, 1,  9,  3, 59,  71,   5, 167,   87, 1507,  193/)
    jkdata(1:12,348) = (/1, 3, 3,  9,  5, 47,  89, 149, 439,  481,  465, 2053/)
    jkdata(1:12,349) = (/1, 3, 5,  9, 23, 15,  35,  35, 307,   85, 2027, 3061/)
    jkdata(1:12,350) = (/1, 3, 1,  5,  9, 27,  53, 119, 235,  799, 1695, 3759/)
    jkdata(1:12,351) = (/1, 3, 3,  5, 25, 19,  73, 183, 473,  917,  367, 1553/)
    jkdata(1:12,352) = (/1, 3, 3,  5,  7, 29,   9,  53,  79,  769,  937, 2007/)
    jkdata(1:12,353) = (/1, 1, 7,  5, 29, 45, 115,  11, 101,  949,  719, 2493/)
    jkdata(1:12,354) = (/1, 3, 3,  1, 11, 35,  49,  13, 245,  739,  545,  603/)
    jkdata(1:12,355) = (/1, 3, 7, 15,  9, 55,  37,   3,  19,  115, 1991, 3343/)
    jkdata(1:12,356) = (/1, 1, 5,  5, 13, 39,   1, 179, 381,  499,   83, 3751/)
    jkdata(1:12,357) = (/1, 3, 3,  9,  5, 19,  35, 229, 251,  945,  819, 1059/)
    jkdata(1:12,358) = (/1, 3, 5, 11, 11, 43,   9,  43,  35,  547,  239,  783/)
    jkdata(1:12,359) = (/1, 3, 3,  7,  1, 21,  45,  55,  25,  225, 1791, 1789/)
    jkdata(1:12,360) = (/1, 3, 1, 15,  3, 19,  81, 187, 107, 1015, 1461, 1589/)
    jkdata(1:12,361) = (/1, 1, 7,  5, 31, 13,  19, 233, 187,  469, 1647,  283/)
    jkdata(1:12,362) = (/1, 1, 1,  3, 27, 17, 127,  47, 115,  737, 1501, 1093/)
    jkdata(1:12,363) = (/1, 1, 7, 13,  3, 51,  17, 133, 113,  495, 1161, 3919/)
    jkdata(1:12,364) = (/1, 1, 7,  5, 17, 37,  17,  91, 321,  353, 1629, 2747/)
    jkdata(1:12,365) = (/1, 1, 1,  3, 27,  5, 105,  47, 115,  103,  139,  277/)
    jkdata(1:12,366) = (/1, 1, 1, 11, 11, 33,  89,  71, 445,   17, 1595, 2605/)
    jkdata(1:12,367) = (/1, 3, 7,  5, 13, 35,  49,  93,  61,  665, 1921, 2169/)
    jkdata(1:12,368) = (/1, 1, 7,  1, 15, 49, 101, 105,  77,  639, 1267, 2905/)
    jkdata(1:12,369) = (/1, 1, 7, 11, 29, 25,   7, 145, 293,  525, 1415,  721/)
    jkdata(1:12,370) = (/1, 3, 5, 13, 15, 45,  37,  45, 405,   75,  509, 4069/)
    jkdata(1:12,371) = (/1, 1, 5,  9,  1,  1,  33, 255,  13,  447,  347,  233/)
    jkdata(1:12,372) = (/1, 1, 1, 11, 15, 63,  11, 221,  53,  185,  777,  261/)
    jkdata(1:12,373) = (/1, 1, 1,  3, 23, 47,  95, 115,  17,   43, 1083, 1137/)
    jkdata(1:12,374) = (/1, 3, 7,  7, 25,  9,  95, 175, 171,  729,  363, 3993/)
    jkdata(1:12,375) = (/1, 1, 5, 13, 13, 63,  17,  19, 299,  577,  269, 3619/)
    jkdata(1:12,376) = (/1, 1, 5, 15, 21, 15, 111, 129,  41,  863, 1015, 2881/)
    jkdata(1:12,377) = (/1, 1, 7,  1, 15, 25, 105,   5,  79,  735, 1809, 1275/)
    jkdata(1:12,378) = (/1, 3, 5,  7,  3, 25,  41, 209,   3,  317, 1105, 3865/)
    jkdata(1:12,379) = (/1, 3, 1, 11, 29, 15, 115, 197, 485,   99, 1429, 1299/)
    jkdata(1:12,380) = (/1, 3, 1,  1, 29, 41,   5,  57, 331,   17, 1471, 3757/)
    jkdata(1:12,381) = (/1, 1, 5, 13,  5, 13,  69, 177,  13,  477, 2019, 1193/)
    jkdata(1:12,382) = (/1, 3, 5,  1, 25,  3, 101, 115, 257,  893,  381,  733/)
    jkdata(1:12,383) = (/1, 1, 5, 15, 17, 19,  27, 187,  59,  537, 2025,  993/)
    jkdata(1:12,384) = (/1, 1, 5,  1, 11, 51,  27, 119, 201,  519, 1223, 1153/)
    jkdata(1:12,385) = (/1, 3, 5,  9,  7, 49, 101,  77, 497, 1017,  827, 2945/)
    jkdata(1:12,386) = (/1, 3, 5,  7, 15, 37, 103, 211,  81,  375, 1733, 3163/)
    jkdata(1:12,387) = (/1, 3, 1,  3,  5, 25,  53, 111, 451,  297,  887, 3179/)
    jkdata(1:12,388) = (/1, 1, 3,  9, 21, 49,   9,  33, 199,  325, 1321,  437/)
    jkdata(1:12,389) = (/1, 3, 1, 11,  7, 13,  21, 113, 171,  999,  803,  271/)
    jkdata(1:12,390) = (/1, 3, 5,  1, 31, 53,  43,  23,  81,  353, 1951, 3493/)
    jkdata(1:12,391) = (/1, 1, 7,  9, 13, 47,  79,  87, 253,  343, 1297, 3971/)
    jkdata(1:12,392) = (/1, 3, 3, 13, 11, 23,  91, 137, 365,  729, 1995, 1005/)
    jkdata(1:12,393) = (/1, 1, 3, 13, 23, 35,  65,  41,  75,  135,  833, 2615/)
    jkdata(1:12,394) = (/1, 3, 5,  3,  5, 29, 117,   7, 451,  489, 1107, 2253/)
    jkdata(1:12,395) = (/1, 3, 7, 11,  7, 33,  87,  83, 149,  859, 1135, 1131/)
    jkdata(1:12,396) = (/1, 1, 3,  7, 23, 21, 125,  43, 483,  267, 1181,  585/)
    jkdata(1:12,397) = (/1, 3, 7,  9, 27, 35,  55, 121,  81,  141, 1251, 2775/)
    jkdata(1:12,398) = (/1, 3, 1,  1, 21, 23,  45, 145, 453,  831,  983, 2171/)
    jkdata(1:12,399) = (/1, 3, 7,  7, 29,  3,  63,   5, 469,  141, 1389, 2383/)
    jkdata(1:12,400) = (/1, 1, 7, 15, 15, 43,  85, 219, 485,  893, 1565, 2937/)
    jkdata(1:12,401) = (/1, 1, 1,  9,  7, 31,  83,  27, 305,  249,  273, 2447/)
    jkdata(1:12,402) = (/1, 3, 3,  1, 27, 63,  97,  11, 163,  807,  137, 1745/)
    jkdata(1:12,403) = (/1, 3, 5,  5, 27,  9,  45, 111, 401,   53,   71,  663/)
    jkdata(1:12,404) = (/1, 1, 1, 13, 19,  1,  83, 207,  15,  613,  735, 1515/)
    jkdata(1:12,405) = (/1, 3, 5,  5,  7, 61,  87,  55,  91,  131, 1005, 3767/)
    jkdata(1:12,406) = (/1, 1, 5, 11, 15, 43, 113,  97,   3,  547,  933, 2709/)
    jkdata(1:12,407) = (/1, 3, 3,  3, 27,  3,  93,  63, 129,  977,   67, 1767/)
    jkdata(1:12,408) = (/1, 1, 7,  9, 27, 11,  95, 229,  35,  131, 1471, 3185/)
    jkdata(1:12,409) = (/1, 1, 3, 15, 19, 55,   5,  53, 239,  999,  551, 3017/)
    jkdata(1:12,410) = (/1, 1, 7, 11, 19, 11,  17,  33, 355,  175,  457, 2815/)
    jkdata(1:12,411) = (/1, 3, 7, 13,  9, 35,  77, 149, 211,   31, 1667, 1829/)
    jkdata(1:12,412) = (/1, 3, 5,  5, 15,  1,  77,  23, 387,  341, 1729,   87/)
    jkdata(1:12,413) = (/1, 3, 7,  1,  1, 63, 127, 187, 101,  739,  919, 3341/)
    jkdata(1:12,414) = (/1, 3, 5,  7,  3, 35, 123, 153, 299,  467,  285,  793/)
    jkdata(1:12,415) = (/1, 1, 7,  7, 29, 49,  45,  91,  67,  675, 1629, 2627/)
    jkdata(1:12,416) = (/1, 3, 1,  5, 29, 19,  81, 193, 375,  241, 1815, 2169/)
    jkdata(1:12,417) = (/1, 1, 1, 13,  5, 45,  85, 183, 405,  645,  653, 1875/)
    jkdata(1:12,418) = (/1, 1, 5,  7, 27,  9, 121,  59, 357,  247, 1919, 3745/)
    jkdata(1:12,419) = (/1, 3, 3,  7, 31, 57, 119, 211, 267,  391, 1039,  367/)
    jkdata(1:12,420) = (/1, 1, 5,  9,  9, 51,  27,  93, 363,  583,  531, 3783/)
    jkdata(1:12,421) = (/1, 3, 1,  5,  1,  1,  85, 139,  79,  183,  393,  783/)
    jkdata(1:12,422) = (/1, 1, 5, 11,  7, 47,  41,  59,  83,  973, 1411,  827/)
    jkdata(1:12,423) = (/1, 1, 3, 11,  3, 41,  49, 179, 437,  433,  359, 3253/)
    jkdata(1:12,424) = (/1, 1, 7,  1, 19,  9,  15, 163, 457,  367,  221, 2639/)
    jkdata(1:12,425) = (/1, 3, 1,  1, 19, 11, 107, 209,  39,  131,  699, 2955/)
    jkdata(1:12,426) = (/1, 1, 5, 15, 29, 37,  21,  77,  97,  467, 1485, 3539/)
    jkdata(1:12,427) = (/1, 3, 7,  3,  9, 19,  51,  39, 473,  571,  471, 1579/)
    jkdata(1:12,428) = (/1, 1, 7, 13,  3, 55, 119, 111, 289,  309, 1357, 2109/)
    jkdata(1:12,429) = (/1, 3, 3,  9, 21, 23,  11,  79, 179,  385, 1715,  379/)
    jkdata(1:12,430) = (/1, 1, 5, 13, 31, 55,  87, 229,  57,  977,  595, 2939/)
    jkdata(1:12,431) = (/1, 3, 1,  9, 29, 55, 101,  85,  23,  111, 1677, 3019/)
    jkdata(1:12,432) = (/1, 3, 3,  9, 25, 13, 115, 237,  49,  917,  153, 1999/)
    jkdata(1:12,433) = (/1, 3, 5, 11,  1,  7,  63, 199,  79,  935, 1903, 2253/)
    jkdata(1:12,434) = (/1, 3, 1,  5,  3, 47,  63, 137,  71,  473, 1281, 2911/)
    jkdata(1:12,435) = (/1, 3, 5,  5,  9, 37,  37, 147, 341,  345,  215, 3733/)
    jkdata(1:12,436) = (/1, 3, 3, 13, 27, 11, 121,  25, 287,  411,  781,  481/)
    jkdata(1:12,437) = (/1, 3, 3, 15,  5, 43, 109,  73,  95,  313,  543, 1767/)
    jkdata(1:12,438) = (/1, 3, 3,  3, 27, 17,   7, 121, 229,   97,  293, 1055/)
    jkdata(1:12,439) = (/1, 1, 7,  9, 25,  3,  43, 129, 271,  149, 1807, 4019/)
    jkdata(1:12,440) = (/1, 3, 3, 15, 21, 25,  69,  83, 475,  959,  965, 4085/)
    jkdata(1:12,441) = (/1, 3, 5,  3, 11, 19,  19,  87,  49,  841, 1695,  105/)
    jkdata(1:12,442) = (/1, 3, 1, 11, 29, 55,  77,  93, 241,  839,  443, 1829/)
    jkdata(1:12,443) = (/1, 3, 3, 11, 31, 59,  49, 205, 261,  669, 1985, 2097/)
    jkdata(1:12,444) = (/1, 3, 7, 15, 27, 37,  71, 167, 495,  431,  321, 2379/)
    jkdata(1:12,445) = (/1, 1, 7, 15, 21, 33,  59,  53, 353,   51,  879, 1567/)
    jkdata(1:12,446) = (/1, 3, 3,  3, 29, 43,  35, 107, 381,   41, 1227, 2713/)
    jkdata(1:12,447) = (/1, 1, 7, 11, 17,  1,   7, 229,  13,  301, 1915,  737/)
    jkdata(1:12,448) = (/1, 3, 5, 15,  9,  5,  13, 213, 291,  247,  839, 3423/)
    jkdata(1:12,449) = (/1, 3, 3, 15, 17, 21,  55,  95,  37, 1015, 1945, 3941/)
    jkdata(1:12,450) = (/1, 3, 3,  3, 13,  5, 101, 219, 251,  377, 1993, 2659/)
    jkdata(1:12,451) = (/1, 1, 1,  1, 11, 63, 127, 109, 105,  329, 1165, 3961/)
    jkdata(1:12,452) = (/1, 3, 7,  3, 25, 49, 103, 175, 399,  945,   51, 1755/)
    jkdata(1:12,453) = (/1, 1, 5,  1, 15, 61,  85,  13,  81,  269,  557, 3613/)
    jkdata(1:12,454) = (/1, 3, 1,  3, 21, 21, 109, 209,  89,   67,  723, 1937/)
    jkdata(1:12,455) = (/1, 1, 1,  3, 11, 51,  29,  97, 265,  979, 1491, 1559/)
    jkdata(1:12,456) = (/1, 3, 3,  1, 19, 15,  61,  61, 507,  581,  817, 2287/)
    jkdata(1:12,457) = (/1, 3, 7,  3, 31, 19,  67, 147, 205,  643, 1237, 2743/)
    jkdata(1:12,458) = (/1, 1, 1, 13,  3, 43,  21,  19, 145,  823,  947,   67/)
    jkdata(1:12,459) = (/1, 3, 7,  1, 19, 47, 111,  13, 331,  557, 1215, 2859/)
    jkdata(1:12,460) = (/1, 3, 1, 11,  5, 17,  67, 123, 129,   91, 1911,  325/)
    jkdata(1:12,461) = (/1, 3, 7,  5,  3,  9,  23,  73, 119,  405, 1225, 2601/)
    jkdata(1:12,462) = (/1, 3, 3, 15,  3, 53,  57,  35, 503,  117, 1965, 1149/)
    jkdata(1:12,463) = (/1, 3, 7,  7,  9, 45,  75, 141, 249,  801, 1889, 3259/)
    jkdata(1:12,464) = (/1, 3, 3, 15, 13, 11,  71,  81,   1,  509, 1503, 2403/)
    jkdata(1:12,465) = (/1, 3, 5,  9, 13, 51, 101,  19, 289,  347, 1177, 3947/)
    jkdata(1:12,466) = (/1, 3, 7,  1,  3, 25, 123, 171, 463,  893,   73, 2011/)
    jkdata(1:12,467) = (/1, 3, 3,  7, 29, 11,  41, 255, 163,  303, 1767,  175/)
    jkdata(1:12,468) = (/1, 1, 5,  1,  7, 25, 107, 111, 443,  227,  303, 3389/)
    jkdata(1:12,469) = (/1, 1, 3,  9,  5, 47, 101, 107,  63,  783,  177, 3915/)
    jkdata(1:12,470) = (/1, 1, 1, 11,  9, 47, 107, 233, 123,  555, 1897, 1315/)
    jkdata(1:12,471) = (/1, 1, 1, 15, 23,  1, 125, 113, 361,  867, 1401, 2447/)
    jkdata(1:12,472) = (/1, 1, 1,  1, 13, 43,  27, 133, 261,   99,  321,  141/)
    jkdata(1:12,473) = (/1, 1, 5, 13, 21, 29,  47,  89,  49,  703,  921,  359/)
    jkdata(1:12,474) = (/1, 3, 7,  9, 23, 17, 119,   9, 429,  111,  217, 3609/)
    jkdata(1:12,475) = (/1, 3, 7, 13, 21, 31,  41, 231, 137,  797, 1779, 3933/)
    jkdata(1:12,476) = (/1, 1, 3, 11, 31, 15,  19,  95, 355,  873,  327,  729/)
    jkdata(1:12,477) = (/1, 1, 3,  7, 11, 59, 127,  69, 175,  541, 1889, 2051/)
    jkdata(1:12,478) = (/1, 3, 1,  3,  7, 27,  33,  33, 507,  919,  333, 1755/)
    jkdata(1:12,479) = (/1, 3, 1,  7,  7, 63,  31,   1,  59,  513,  615, 2149/)
    jkdata(1:12,480) = (/1, 1, 1,  3,  3, 11, 109, 253, 277,  343, 1665, 2107/)
    jkdata(1:12,481) = (/1, 1, 5, 13, 23, 41,   7, 219, 391,  319, 1825, 1741/)
    jkdata(1:13,482) = (/1, 1, 5,  7,  1, 51,  91, 253,  25,  517, 1639, 1051, 2319/)
    jkdata(1:13,483) = (/1, 3, 7,  9, 23, 29,  91, 247, 185,  135,  237, 3681,  653/)
    jkdata(1:13,484) = (/1, 3, 3,  7,  5,  7,  39, 129, 381,  871, 1205,  471, 1379/)
    jkdata(1:13,485) = (/1, 1, 1,  7,  9, 27, 125,  11, 197,  917,  361, 1055, 1675/)
    jkdata(1:13,486) = (/1, 1, 1,  3, 17, 63, 105, 251,  39,  285,  129,  845, 1951/)
    jkdata(1:13,487) = (/1, 3, 3,  3, 21, 31,  47, 221,   5,  663, 1655,  257, 7075/)
    jkdata(1:13,488) = (/1, 3, 3,  9,  1, 43, 125, 153, 429,  301,  983, 1559, 2087/)
    jkdata(1:13,489) = (/1, 3, 7,  9, 17,  3, 123,  35, 119,   15, 1089, 1061, 7147/)
    jkdata(1:13,490) = (/1, 3, 3,  7, 29, 29,  91, 103, 247,  763, 1171, 2803, 1427/)
    jkdata(1:13,491) = (/1, 1, 3,  5,  7, 39,   9, 239, 177,   89,  401, 2219,  893/)
    jkdata(1:13,492) = (/1, 1, 5, 11,  5,  3, 103,   7, 329,  323,  677, 1315,  171/)
    jkdata(1:13,493) = (/1, 3, 1, 13, 17, 59,  45,  27, 465,  757,  643, 1369, 2019/)
    jkdata(1:13,494) = (/1, 1, 3, 13, 13, 59,  23, 235, 421,  317,  749, 3211, 7235/)
    jkdata(1:13,495) = (/1, 3, 7,  7, 25,  1, 117, 181, 271,  807,  303, 4027, 5697/)
    jkdata(1:13,496) = (/1, 3, 3,  7, 17, 53,   9,   5, 467,  309, 1407,  105, 3615/)
    jkdata(1:13,497) = (/1, 1, 3, 15,  9, 63, 125, 207, 151, 1013, 1873,   11, 1961/)
    jkdata(1:13,498) = (/1, 3, 7,  9, 19, 23,  73,  53,  45,  345, 1579, 1077, 7517/)
    jkdata(1:13,499) = (/1, 3, 3,  5,  9, 63,  11, 149, 429,  499, 1491, 2857, 6849/)
    jkdata(1:13,500) = (/1, 1, 5,  5,  5, 47,  37, 155, 137,  279, 1393,  337, 2893/)
    jkdata(1:13,501) = (/1, 1, 7,  3,  7, 51,  61, 225, 471,  711, 1247, 3553, 1883/)
    jkdata(1:13,502) = (/1, 1, 5,  3, 21, 23,  79, 165,  11,  915,  789, 3503, 2863/)
    jkdata(1:13,503) = (/1, 3, 7, 13, 19, 61,  21, 137,  17,  411,  763, 3917, 2173/)
    jkdata(1:13,504) = (/1, 3, 7,  3, 13, 39,   5, 155, 409,  281,   49, 2665, 4543/)
    jkdata(1:13,505) = (/1, 3, 3,  9,  9, 47,  47, 201, 347,  193,    5, 3823,   73/)
    jkdata(1:13,506) = (/1, 1, 3,  3,  7, 21, 117,  97, 199,  739, 1607, 3403,  381/)
    jkdata(1:13,507) = (/1, 1, 5,  1,  3, 39,  67, 245, 463,  365, 1891, 3711, 3893/)
    jkdata(1:13,508) = (/1, 3, 1, 11,  9, 15,  53, 203, 177,  315,  735, 2085, 6045/)
    jkdata(1:13,509) = (/1, 3, 3,  1,  3,  3,  85,  47,  11,  375, 1557, 1103, 1643/)
    jkdata(1:13,510) = (/1, 3, 5,  3, 15,  9,  33,  39,  51,  809, 1909, 1641, 7669/)
    jkdata(1:13,511) = (/1, 3, 3, 11, 31, 57,  81,  35, 361,  469, 1765,  701, 1027/)
    jkdata(1:13,512) = (/1, 3, 1, 15, 29, 61, 121, 105,  95,  487, 1777, 4095, 1549/)
    jkdata(1:13,513) = (/1, 1, 3, 11, 29, 39,  47, 239, 497,  621, 1127, 2883, 3983/)
    jkdata(1:13,514) = (/1, 1, 5, 11, 25, 37,  61,  49, 163,  857,  813, 1435, 1985/)
    jkdata(1:13,515) = (/1, 1, 1, 11, 13, 21,  51,  15, 351,  975,  695,  653, 6589/)
    jkdata(1:13,516) = (/1, 3, 1,  9,  9, 51, 127, 253, 127,  537,   97, 2363, 7497/)
    jkdata(1:13,517) = (/1, 1, 3, 13, 21,  1,  29,   7, 395,  939,  731, 1597, 2745/)
    jkdata(1:13,518) = (/1, 3, 7,  7,  9, 23,  65, 237, 511,  585, 1503,  767, 2375/)
    jkdata(1:13,519) = (/1, 3, 7,  9, 31, 43,  45, 213, 327,  129, 1751,  869, 7047/)
    jkdata(1:13,520) = (/1, 1, 1, 15,  7, 27,  41,  55, 353,  625,  333, 1825, 1117/)
    jkdata(1:13,521) = (/1, 3, 5,  9, 15, 25,  95,  87,  49,  447,  769, 1117, 1171/)
    jkdata(1:13,522) = (/1, 3, 1, 11,  5, 11,  57, 199, 105,  129,  865, 1297, 1975/)
    jkdata(1:13,523) = (/1, 3, 3,  1, 31, 13,  73,  27, 151, 1017,  693,  501, 5199/)
    jkdata(1:13,524) = (/1, 3, 7,  3,  7, 21,  33, 175, 321,  133,  377,  505, 3915/)
    jkdata(1:13,525) = (/1, 1, 3,  3, 15, 43, 117,  49, 331,   83, 1919,  149, 3695/)
    jkdata(1:13,526) = (/1, 1, 7,  9, 27,  7,  61,  41, 329,    3,  957,  873, 8113/)
    jkdata(1:13,527) = (/1, 3, 3,  7, 25, 11, 111, 229, 509,  415, 1359, 2673, 4303/)
    jkdata(1:13,528) = (/1, 1, 5, 15, 19, 33,  59,  85, 107,  661, 1627,  551, 3773/)
    jkdata(1:13,529) = (/1, 1, 1, 13,  9, 55, 123,   3, 109,   53, 1039, 1499, 7705/)
    jkdata(1:13,530) = (/1, 3, 7, 13,  9,  1,  65, 149, 303,  115, 1783, 2793, 6855/)
    jkdata(1:13,531) = (/1, 1, 1,  7, 25, 37,  47, 179, 467,  903, 1065, 3277, 1675/)
    jkdata(1:13,532) = (/1, 3, 1, 15, 25, 35, 105, 129, 287,   49, 1665, 2143, 2245/)
    jkdata(1:13,533) = (/1, 1, 3,  9, 23, 27,  23, 185, 161,   79, 1917, 3663, 2817/)
    jkdata(1:13,534) = (/1, 3, 5, 13,  1, 61,  29, 249,  45,   55, 1947,  533, 1719/)
    jkdata(1:13,535) = (/1, 1, 3,  9,  9, 39, 107, 197, 385,  385,  991, 3991,  569/)
    jkdata(1:13,536) = (/1, 3, 7, 15,  7,  5,  37,  15, 289,  261, 1997,  575, 1021/)
    jkdata(1:13,537) = (/1, 3, 1, 13, 11, 19,  81,  97, 363,  345,  841, 1877, 2077/)
    jkdata(1:13,538) = (/1, 1, 5, 15, 15, 61,  67, 197, 331,  297,  459, 1009, 5945/)
    jkdata(1:13,539) = (/1, 1, 5,  9, 19, 61,  29, 139, 265,  199,  221, 3929, 1833/)
    jkdata(1:13,540) = (/1, 3, 1, 13, 15, 57, 115, 203, 407,  385,  327,  473, 2631/)
    jkdata(1:13,541) = (/1, 3, 1,  1, 27, 59, 119,  63,  37,  617, 1595, 3009, 4851/)
    jkdata(1:13,542) = (/1, 1, 3, 11, 17, 21,  75,  33, 433,   25, 1881, 2595, 6371/)
    jkdata(1:13,543) = (/1, 3, 1,  7, 11, 59,  73, 251, 315,  515, 1269, 3249,  833/)
    jkdata(1:13,544) = (/1, 3, 3, 11, 11, 61,  99, 217, 343,  275, 1007,  675, 7987/)
    jkdata(1:13,545) = (/1, 1, 3,  3, 31, 57, 103, 199,  63,  849,  129, 3593,  331/)
    jkdata(1:13,546) = (/1, 3, 7, 13, 13, 25,   7, 199,  51,  401, 1413, 2453, 1899/)
    jkdata(1:13,547) = (/1, 3, 1,  5, 25, 55,  57,  99, 185,  471,  475, 1567, 8093/)
    jkdata(1:13,548) = (/1, 1, 7,  1, 25, 27,  45, 249,  71,  377, 1105,  973, 6719/)
    jkdata(1:13,549) = (/1, 1, 3,  7,  9, 31,  61,  33,  27,  661,  791,  595, 6903/)
    jkdata(1:13,550) = (/1, 3, 1, 15,  7, 41,  95, 229, 267,  535, 1983, 1335, 5903/)
    jkdata(1:13,551) = (/1, 1, 7,  3, 13, 33,  49, 177, 503,  505, 1359, 1715, 5657/)
    jkdata(1:13,552) = (/1, 3, 3, 13, 29, 63, 101,  13, 239,  939,  503,  589, 5007/)
    jkdata(1:13,553) = (/1, 3, 1,  7, 19, 19, 101, 209, 293,  465,  691,   85, 2689/)
    jkdata(1:13,554) = (/1, 1, 7, 13,  5, 57,  35, 147, 245,  225,  659, 2265, 6637/)
    jkdata(1:13,555) = (/1, 1, 3, 13, 19, 35,  47,  97, 281,  929,  691, 3069, 2675/)
    jkdata(1:13,556) = (/1, 3, 5, 11, 31, 13, 119,  31, 297,  219,  343,  461, 1645/)
    jkdata(1:13,557) = (/1, 1, 3,  3, 25, 63,  39, 125,  75,  955, 1375, 1659, 1819/)
    jkdata(1:13,558) = (/1, 3, 5,  5, 13, 35,  67, 177, 461,  659, 1919, 2627,  689/)
    jkdata(1:13,559) = (/1, 1, 7,  3, 25, 17,  31, 137, 371,  441,  263, 1307, 6709/)
    jkdata(1:13,560) = (/1, 3, 3, 13, 15, 11, 103, 187, 129,  117, 1373, 1731, 7717/)
    jkdata(1:13,561) = (/1, 1, 3, 11,  5, 11,   7,  11, 189,  527,  603, 1501, 6295/)
    jkdata(1:13,562) = (/1, 1, 3,  9,  9, 49,  61,  91, 189,  427, 1383, 1699, 7013/)
    jkdata(1:13,563) = (/1, 3, 5,  9, 29, 41, 127, 223, 339,  515,  297, 3545, 7695/)
    jkdata(1:13,564) = (/1, 3, 1,  3, 31, 55,  87,  29, 287,  287,  781, 3803, 3705/)
    jkdata(1:13,565) = (/1, 1, 7, 11,  9,  5,   3, 169, 111,  191,  145, 2157, 7069/)
    jkdata(1:13,566) = (/1, 1, 7, 11, 29, 45,  35, 231, 111,   33,  285,  453, 2621/)
    jkdata(1:13,567) = (/1, 1, 1,  7, 27, 17,  29,  59, 379,  389,  767, 2813, 3631/)
    jkdata(1:13,568) = (/1, 3, 3,  9, 25, 35,  73,  31,  93,  197, 1739, 2047, 6571/)
    jkdata(1:13,569) = (/1, 3, 1, 13, 27,  5,  95, 163,  27,  825, 1715, 2999, 6259/)
    jkdata(1:13,570) = (/1, 1, 3, 11, 11, 31, 103,  41, 185,   63,  715, 3841, 7261/)
    jkdata(1:13,571) = (/1, 3, 7,  7, 17, 31,  71,  57, 347,  417,  317, 2361, 3397/)
    jkdata(1:13,572) = (/1, 1, 7, 15,  5, 37,  75,  87, 337,  949, 1333, 1079, 7645/)
    jkdata(1:13,573) = (/1, 1, 1, 13, 17, 17,  51, 247, 247,   35,   85,  573, 1115/)
    jkdata(1:13,574) = (/1, 3, 3,  7,  3, 45,  87,  25, 507,  571,  831,   69, 4753/)
    jkdata(1:13,575) = (/1, 3, 7,  5, 23, 51,  57, 127, 161,    9, 1615, 1363, 2047/)
    jkdata(1:13,576) = (/1, 1, 3,  3, 15,  1,  97, 101, 231,  131,   81, 1597, 7579/)
    jkdata(1:13,577) = (/1, 1, 1,  1,  9, 39,  11, 207,  43,  609, 1667, 3427, 2271/)
    jkdata(1:13,578) = (/1, 3, 5,  5,  9, 49, 105, 187, 499,  439, 1467, 2899, 5403/)
    jkdata(1:13,579) = (/1, 1, 3, 15, 17, 55,  87,  73,  73,   95, 1457, 2771, 4911/)
    jkdata(1:13,580) = (/1, 3, 1, 15, 17, 19,  41,  61, 327,   19, 1453, 1327, 7629/)
    jkdata(1:13,581) = (/1, 1, 1,  3, 31, 41,  73, 105, 263,  569, 1825, 1117, 4225/)
    jkdata(1:13,582) = (/1, 1, 1, 11, 11, 13, 109,  27, 331,  893,  109, 1523, 1209/)
    jkdata(1:13,583) = (/1, 1, 5,  1, 19,  5,  69,  91, 249,  451,  387, 3521, 6955/)
    jkdata(1:13,584) = (/1, 1, 3,  7, 25, 51,  35, 171, 493,  397, 1207, 2393, 6951/)
    jkdata(1:13,585) = (/1, 1, 3,  3, 13,  5, 121, 243,  37,  971, 2039, 2537, 1829/)
    jkdata(1:13,586) = (/1, 3, 7, 15, 23, 49,  39,  33,  25,  801,  213, 1979, 5579/)
    jkdata(1:13,587) = (/1, 1, 1, 11, 15,  1, 111,   3, 115,  125, 1351, 3179, 5231/)
    jkdata(1:13,588) = (/1, 1, 5,  5, 25, 21,   1,   1,   3,  471, 1329,  683, 1783/)
    jkdata(1:13,589) = (/1, 1, 3,  5, 21, 13,  77,  21, 167,  187, 1173, 2453, 4285/)
    jkdata(1:13,590) = (/1, 1, 5,  3, 31, 17,  39, 229, 197,  257,   57,  453, 7425/)
    jkdata(1:13,591) = (/1, 3, 1,  5, 19, 59,  47,  93, 127,   67, 1769, 1227,  599/)
    jkdata(1:13,592) = (/1, 1, 3,  5,  3, 51,  53,  71, 357,  949,  951,  779, 5785/)
    jkdata(1:13,593) = (/1, 3, 1,  1, 11, 11,  91,  61, 497,  621,  183,  671, 3275/)
    jkdata(1:13,594) = (/1, 1, 3, 15, 25,  3,   3,  37, 103,  453,   23, 3483, 5643/)
    jkdata(1:13,595) = (/1, 1, 1,  5,  7, 61,  17, 183, 125,  411,  451, 2135, 2263/)
    jkdata(1:13,596) = (/1, 3, 5,  1, 15,  1,  51,  65, 191,  621, 1155, 3139,  657/)
    jkdata(1:13,597) = (/1, 3, 7,  5, 19, 33,  83, 211, 165,  955, 1551, 3381, 6769/)
    jkdata(1:13,598) = (/1, 1, 7,  3,  7, 37,  39,  53,  55,  309, 2037, 3945, 6261/)
    jkdata(1:13,599) = (/1, 1, 1,  7,  5, 33, 125,  11, 101,  783,  811,   57, 1251/)
    jkdata(1:13,600) = (/1, 3, 1,  5,  3, 61,  85, 151,  95,  893,  635, 1541, 3249/)
    jkdata(1:13,601) = (/1, 1, 5, 11, 13, 25, 111, 165,  79,  597, 1671, 3405, 4447/)
    jkdata(1:13,602) = (/1, 3, 3,  3, 13, 27,  21,  47, 351,  377, 1451, 3381, 4111/)
    jkdata(1:13,603) = (/1, 1, 1, 13,  1, 59,  69,   5, 341,  753,  863, 2371, 3991/)
    jkdata(1:13,604) = (/1, 3, 5,  9, 23,  7,  85, 129,  43,  145, 1499, 2879, 1215/)
    jkdata(1:13,605) = (/1, 3, 1, 13,  5, 49,  29,  79, 125,  637, 1673, 1985,  131/)
    jkdata(1:13,606) = (/1, 3, 1, 15, 25, 13,  55, 101, 135,  941,  363,  987, 4397/)
    jkdata(1:13,607) = (/1, 1, 7,  5, 11, 63,  11, 147, 173,  593, 1029, 3017, 3487/)
    jkdata(1:13,608) = (/1, 3, 7,  3, 25,  3, 117, 169, 289,  317, 1077, 3031, 7585/)
    jkdata(1:13,609) = (/1, 3, 3,  5, 15, 33,   1, 181, 373,  555, 1525, 3839, 5565/)
    jkdata(1:13,610) = (/1, 3, 5,  9, 13,  3,  47,  19, 133,  375,  277, 1401, 7199/)
    jkdata(1:13,611) = (/1, 1, 5,  5, 21, 15,  17,  95, 421,  575, 1023, 3749, 3573/)
    jkdata(1:13,612) = (/1, 1, 1,  3, 11,  9,  65,  77, 241,  175,  655, 2977, 7105/)
    jkdata(1:13,613) = (/1, 3, 7, 11, 23, 13,  63, 139, 281,  403,  665,  681, 7409/)
    jkdata(1:13,614) = (/1, 3, 1,  1, 29, 35,  47, 197, 213,  571, 1869, 1175, 1671/)
    jkdata(1:13,615) = (/1, 3, 5, 13,  5, 39, 117, 219, 177,  555, 1255, 1519,  949/)
    jkdata(1:13,616) = (/1, 1, 1,  9, 17, 11,  17,  97, 363,  109,  965, 3355, 3889/)
    jkdata(1:13,617) = (/1, 1, 1, 15, 27, 59, 115, 239, 151,  377,  277,  907, 5971/)
    jkdata(1:13,618) = (/1, 1, 3,  3,  9, 59,  51, 183, 227,  931, 1601,  117, 3333/)
    jkdata(1:13,619) = (/1, 1, 1,  5, 19,  1,  25, 143, 145,  499,  329,  771,  225/)
    jkdata(1:13,620) = (/1, 3, 5, 11, 15, 57,  33,   9, 363,  649, 1603, 3741, 3647/)
    jkdata(1:13,621) = (/1, 1, 7,  9,  5, 11, 123,  13, 239,  653, 1901, 3337, 5403/)
    jkdata(1:13,622) = (/1, 3, 5,  1, 29,  5, 123, 209, 431,  329,  395, 1743, 3409/)
    jkdata(1:13,623) = (/1, 1, 7,  3, 23, 57,  83,  23,  81,  279,   65, 1227, 7459/)
    jkdata(1:13,624) = (/1, 3, 7, 15, 19, 13,  51, 215, 397,  271, 1307, 3335, 6879/)
    jkdata(1:13,625) = (/1, 1, 1,  9,  1, 31, 113,  53, 241,  647, 2029, 2755, 5789/)
    jkdata(1:13,626) = (/1, 1, 5,  9, 27, 13,  95, 137,  67,  721,   21, 1909, 6567/)
    jkdata(1:13,627) = (/1, 3, 1,  9,  3, 11, 121, 203, 291,  665, 1321, 3603, 5581/)
    jkdata(1:13,628) = (/1, 3, 1, 11, 23, 55,  51,  19, 255,  429,  543, 2397, 4919/)
    jkdata(1:13,629) = (/1, 1, 3,  7, 21, 45,  91, 151, 405,  957, 1569,  653, 1927/)
    jkdata(1:13,630) = (/1, 1, 5,  5, 19,  9, 109, 171, 421,  803, 1185,   87, 4407/)
    jkdata(1:13,631) = (/1, 1, 1, 13, 27, 55,  43, 133, 399,  767, 1905, 2025, 8085/)
    jkdata(1:13,632) = (/1, 3, 5,  1, 11, 55,  55, 219,  75,  425, 1701, 2617, 4691/)
    jkdata(1:13,633) = (/1, 3, 5, 15, 17, 19,  35, 231, 399,  477,  413, 3257,  611/)
    jkdata(1:13,634) = (/1, 1, 3,  3, 13, 25,  55,   3, 105,  995, 2041,  287, 3005/)
    jkdata(1:13,635) = (/1, 3, 1, 13, 27, 41,  87,  15, 329,  105, 1697, 3051,  591/)
    jkdata(1:13,636) = (/1, 1, 3,  9, 11, 23,  33, 253,  41,  495,  725, 3809,  753/)
    jkdata(1:13,637) = (/1, 3, 1, 13, 31, 45,  37, 225, 425,  575, 1417,  897,  589/)
    jkdata(1:13,638) = (/1, 1, 5,  5, 23, 29,   5,  33,   7,  687, 1847, 2215,  171/)
    jkdata(1:13,639) = (/1, 1, 5,  1,  5, 63,   3, 111, 283,  385,  411,   63, 5729/)
    jkdata(1:13,640) = (/1, 1, 3,  5,  9, 59,  45, 183, 375,  227,  211, 2043, 5891/)
    jkdata(1:13,641) = (/1, 1, 3,  1, 21, 27,  21, 213, 475,  923,  915, 1757, 1033/)
    jkdata(1:13,642) = (/1, 1, 3, 13, 31, 39, 105, 169, 427,  563, 1891, 3671, 3049/)
    jkdata(1:13,643) = (/1, 1, 3, 13, 29, 21, 127, 119, 277,  723,   17,  297, 6567/)
    jkdata(1:13,644) = (/1, 3, 1,  7, 11, 37,  35, 111, 209,  481, 1877, 3131, 5257/)
    jkdata(1:13,645) = (/1, 1, 1,  7, 21,  7,  17,  15, 411,  717, 1699, 1305, 8003/)
    jkdata(1:13,646) = (/1, 3, 3,  1, 17, 61,  35, 201,   3,  111,  687,  293, 1757/)
    jkdata(1:13,647) = (/1, 3, 1,  9, 15, 49,  37, 123, 137,  633, 1089, 3865, 4489/)
    jkdata(1:13,648) = (/1, 1, 3,  5,  7, 35,  97, 121, 195,  113, 1973, 3173, 4923/)
    jkdata(1:13,649) = (/1, 3, 5, 11, 15, 39,  97, 225, 289,  369, 1809, 3397, 6379/)
    jkdata(1:13,650) = (/1, 3, 5,  9,  7,  9,  21, 113, 509,  955,  851, 2269, 5171/)
    jkdata(1:13,651) = (/1, 3, 7, 11,  9, 29,  77, 113, 121,  253, 1495, 3673, 1757/)
    jkdata(1:13,652) = (/1, 1, 5, 13, 21,  7, 123, 225,  55,  321, 1257,  717,  689/)
    jkdata(1:13,653) = (/1, 3, 5,  3, 27, 25,  17, 161, 147,  409,   63, 3041, 3081/)
    jkdata(1:13,654) = (/1, 1, 7, 15, 25, 23,  89, 165, 275,  909, 1323, 3341, 1389/)
    jkdata(1:13,655) = (/1, 1, 5, 15, 29, 57,  53,   1, 251,  367, 1307, 3595, 4113/)
    jkdata(1:13,656) = (/1, 3, 7, 13, 11,  5, 105, 139,  19,   33,  609, 3819,  455/)
    jkdata(1:13,657) = (/1, 3, 1, 15,  3, 19,  75,  55, 129,  967,  881, 2871, 2761/)
    jkdata(1:13,658) = (/1, 1, 3,  7, 21, 15,  25,   3, 285,  453, 1543, 3973,  847/)
    jkdata(1:13,659) = (/1, 1, 7,  5, 13, 33, 125,  93, 415,  863,  177, 1129, 7575/)
    jkdata(1:13,660) = (/1, 3, 7,  7, 23, 49,  13, 217, 487,  449,  617,  513, 5829/)
    jkdata(1:13,661) = (/1, 3, 3,  9, 19, 37,  47, 193, 491,  539, 1505,  871,  633/)
    jkdata(1:13,662) = (/1, 1, 5,  7, 27, 25,  21,  97, 193,  781, 1747, 1485, 6629/)
    jkdata(1:13,663) = (/1, 1, 5,  9, 17, 17, 125,  29, 219,  911, 1537, 3977, 1103/)
    jkdata(1:13,664) = (/1, 1, 7,  9, 29, 45,  23,  69, 403,  113,  925, 2473, 7635/)
    jkdata(1:13,665) = (/1, 3, 5,  9, 25, 29,  55, 231,  23,    7,  183, 1171,  803/)
    jkdata(1:13,666) = (/1, 1, 5, 11, 17, 15,  63, 161,  97,  219,   77, 1143, 6175/)
    jkdata(1:13,667) = (/1, 3, 3,  9,  9, 25,  61,  93,  65,  725, 1723, 3063, 6587/)
    jkdata(1:13,668) = (/1, 3, 3,  3,  1,  3,   5,  69, 285, 1015, 1877, 3547, 2711/)
    jkdata(1:13,669) = (/1, 1, 3, 11, 19,  3,  17, 143,  75,  971, 1703, 2183, 3879/)
    jkdata(1:13,670) = (/1, 1, 1, 15, 23, 49,  93, 137,  21, 1021,  397, 3993,   67/)
    jkdata(1:13,671) = (/1, 3, 7, 13,  5, 11,  57,   9, 373,  525,  459,  133, 1179/)
    jkdata(1:13,672) = (/1, 1, 1, 13, 23, 39, 121,  87, 261,  785,  521, 2529, 4761/)
    jkdata(1:13,673) = (/1, 1, 5,  5,  1, 15,  69, 183, 339,  873,  257, 2699, 7281/)
    jkdata(1:13,674) = (/1, 3, 5,  9, 17, 19,  73, 113, 239,  191, 1177,  233, 1557/)
    jkdata(1:13,675) = (/1, 1, 5, 15, 17, 57,  93, 183, 495,  893,  389, 2355, 3379/)
    jkdata(1:13,676) = (/1, 3, 3,  1, 13, 39, 121,  73, 415,  297, 1947,  231, 2459/)
    jkdata(1:13,677) = (/1, 1, 3,  1, 27, 15, 105, 215, 333,  507, 1553, 3241, 4273/)
    jkdata(1:13,678) = (/1, 1, 5,  9, 23, 11,  75, 137, 107,  215, 1583,  611, 4127/)
    jkdata(1:13,679) = (/1, 1, 1,  5,  7,  3,  91,  89, 435,   21, 1831, 1309, 7147/)
    jkdata(1:13,680) = (/1, 3, 3, 13,  7, 57,  67, 251, 297,  153,  261, 3829,   35/)
    jkdata(1:13,681) = (/1, 3, 1,  3, 11, 31,  95, 163, 213,  645,  485, 1839, 3549/)
    jkdata(1:13,682) = (/1, 3, 3, 13, 13, 55,  75,  41, 149,  913,  289, 1495,  395/)
    jkdata(1:13,683) = (/1, 3, 3, 15, 17, 61,   9, 227, 463,  755, 1281,  301, 3735/)
    jkdata(1:13,684) = (/1, 1, 3,  3, 13, 19,  69, 145, 199,  371, 1543, 1169, 5787/)
    jkdata(1:13,685) = (/1, 1, 7,  1, 11,  5,  97,  57, 323,  881, 1591, 1613, 4179/)
    jkdata(1:13,686) = (/1, 3, 1,  3, 21, 41,  99,  81,  45,  113, 1123, 2673, 5889/)
    jkdata(1:13,687) = (/1, 3, 7, 11, 13, 35,  93,  57,  19,  903,  573,  243, 5057/)
    jkdata(1:13,688) = (/1, 1, 7, 13, 23, 59,  11,  11, 301,  225,  821, 3601, 7473/)
    jkdata(1:13,689) = (/1, 1, 3,  1,  1, 61,  53, 135, 121,   49, 1065, 3669, 4713/)
    jkdata(1:13,690) = (/1, 1, 7, 15, 27, 39,  19, 145, 499,  587, 1933, 2813, 2133/)
    jkdata(1:13,691) = (/1, 1, 1,  9, 13, 41,  73, 161, 187,  201, 1373, 2671, 2897/)
    jkdata(1:13,692) = (/1, 3, 1,  9,  9, 53,   5, 175, 229,  927, 2005, 2679, 1841/)
    jkdata(1:13,693) = (/1, 1, 5,  3,  7, 53,  33, 159,  63,  429,  905, 3463, 2125/)
    jkdata(1:13,694) = (/1, 1, 7,  1,  1, 63,  79,  25, 425,  599,  207, 2477, 1029/)
    jkdata(1:13,695) = (/1, 3, 1,  9, 27, 31, 107,  55,  99,  513,  173, 1795, 1695/)
    jkdata(1:13,696) = (/1, 3, 7,  1, 29,  9,  65, 167, 281,   97, 1573,  617, 6523/)
    jkdata(1:13,697) = (/1, 3, 1,  9,  5, 59,  69, 157,  35,  319, 1597, 2317, 1143/)
    jkdata(1:13,698) = (/1, 1, 7,  1, 13, 13,  79, 211, 125,  331,  573, 1855, 5105/)
    jkdata(1:13,699) = (/1, 1, 7, 13, 25, 35, 125,  97, 349,  833, 1883, 1057, 7133/)
    jkdata(1:13,700) = (/1, 3, 1, 11, 21, 55,  25, 247,  87,  325, 1795, 1703, 3351/)
    jkdata(1:13,701) = (/1, 3, 3, 15,  3, 41,  93, 249, 101,  887, 1499, 1761, 2775/)
    jkdata(1:13,702) = (/1, 1, 7,  7, 31, 49,  55,  23,  59,  139, 1743, 2515, 3971/)
    jkdata(1:13,703) = (/1, 3, 5, 11, 15,  5,  61, 129, 195,  927,  553,  801, 4503/)
    jkdata(1:13,704) = (/1, 3, 1, 15, 13, 41,  17, 159, 511,  399,  335, 1205, 7589/)
    jkdata(1:13,705) = (/1, 1, 3, 13,  3, 25, 117,  71, 355,  163,  333, 1311, 5155/)
    jkdata(1:13,706) = (/1, 1, 5, 15, 19, 27,  69, 197,  73,  307, 1645,  473, 4305/)
    jkdata(1:13,707) = (/1, 3, 5,  1, 13, 43,  97, 127, 263,  803,  791, 3963, 1641/)
    jkdata(1:13,708) = (/1, 3, 5,  9,  1,  5,  87, 141, 243,  169,  871,  697, 4717/)
    jkdata(1:13,709) = (/1, 3, 1,  9, 27,  5, 111, 219, 101, 1019, 1157, 1221, 2427/)
    jkdata(1:13,710) = (/1, 3, 1,  7, 15, 43,  37,   5, 165,  869,  969,  251, 5617/)
    jkdata(1:13,711) = (/1, 3, 7,  3, 17,  5,  93, 233, 141,  537,  557,  381, 1267/)
    jkdata(1:13,712) = (/1, 3, 1,  5,  1,  5,  59, 131,  11,  907,  141, 3887,  399/)
    jkdata(1:13,713) = (/1, 3, 7, 11,  3, 17,  79, 217, 389,  479,  223, 1761, 5831/)
    jkdata(1:13,714) = (/1, 1, 1,  7, 13,  5,  95, 101, 219,  335, 1129, 3093, 4305/)
    jkdata(1:13,715) = (/1, 3, 7,  3, 13, 15,  53, 131, 187,  697, 1685, 3721, 4241/)
    jkdata(1:13,716) = (/1, 3, 7,  9, 13, 27, 115,  33, 449,  479,  423, 2079, 3395/)
    jkdata(1:13,717) = (/1, 1, 3,  5, 31, 29,  53, 157, 447,  353, 1069, 4085, 3045/)
    jkdata(1:13,718) = (/1, 3, 1, 15, 29, 17,  85, 173, 393,  769,  391,  379, 4899/)
    jkdata(1:13,719) = (/1, 1, 1,  7, 27,  9,  85,  69, 477,  787,   99, 3601, 1713/)
    jkdata(1:13,720) = (/1, 3, 5,  5,  7,  3,  65, 207, 305, 1023,   95, 3845,  171/)
    jkdata(1:13,721) = (/1, 1, 1,  3,  7, 55,  59, 239, 221,  855, 1847,  433,  411/)
    jkdata(1:13,722) = (/1, 1, 5, 13, 21, 31,  23,  81,  51,  493,  531, 1781, 7099/)
    jkdata(1:13,723) = (/1, 3, 1,  7, 29,  1,  75, 205, 355,  883, 1859,   29, 5473/)
    jkdata(1:13,724) = (/1, 3, 5,  1, 15, 45,  21,  11, 209,  521, 1833, 1897, 5209/)
    jkdata(1:13,725) = (/1, 1, 3,  1, 17, 45,  67,  41, 499,  735, 1833, 1599, 1195/)
    jkdata(1:13,726) = (/1, 1, 5,  9, 17, 13,  27, 169, 479,  297,  341, 2163, 1077/)
    jkdata(1:13,727) = (/1, 1, 5, 15, 21, 57,  99,  65, 265, 1011,  237,   75, 1309/)
    jkdata(1:13,728) = (/1, 3, 5, 15, 19, 17,  79, 193, 377,  991, 1997, 3475, 2953/)
    jkdata(1:13,729) = (/1, 1, 5, 15, 17,  3,  27,  77, 145,  879, 1799, 3957, 7343/)
    jkdata(1:13,730) = (/1, 3, 5, 11,  3, 61,   3, 201, 411,  855,  409, 1641, 4887/)
    jkdata(1:13,731) = (/1, 3, 3,  3, 15, 15,  95, 173, 173,  591,  431, 3911, 3229/)
    jkdata(1:13,732) = (/1, 1, 3,  5,  5, 49,  27,   1,  11,  415, 1917, 2959, 6759/)
    jkdata(1:13,733) = (/1, 3, 7, 15, 27, 15,  69, 221, 433,  917,  363, 2833, 6721/)
    jkdata(1:13,734) = (/1, 3, 3, 13, 27, 47,  19, 157, 483,  375,  335, 1279, 6775/)
    jkdata(1:13,735) = (/1, 1, 3,  7,  3,  9,  75,   1, 135,  453, 1039, 1099,  675/)
    jkdata(1:13,736) = (/1, 3, 5, 15, 31, 37,  47,  15, 385,  553, 1085,  403, 4039/)
    jkdata(1:13,737) = (/1, 1, 5, 15, 31, 45,  59, 113, 341,  189, 1657,  799, 2493/)
    jkdata(1:13,738) = (/1, 1, 3, 11,  7,  9,  41, 147,  89,  841, 1975, 2183, 7511/)
    jkdata(1:13,739) = (/1, 3, 7, 11, 21, 51,  85, 137, 209,  339, 1527, 2699, 3269/)
    jkdata(1:13,740) = (/1, 3, 1,  9,  3, 61,  77, 205, 391,  211, 1111, 1711, 4199/)
    jkdata(1:13,741) = (/1, 3, 5,  5, 13, 21,  99, 225,  33,  601,  659, 2037, 6625/)
    jkdata(1:13,742) = (/1, 1, 7, 15, 11, 33,  55,  73, 395,   57,  389,  727, 7943/)
    jkdata(1:13,743) = (/1, 1, 5,  9, 17, 11,  49,  45, 319,  765,  899,  289, 2013/)
    jkdata(1:13,744) = (/1, 1, 1,  7, 27, 21,  93,  49, 451,  745,  595, 1785, 4145/)
    jkdata(1:13,745) = (/1, 3, 5,  3, 25, 63,  93, 149, 119,  621, 1439, 1575,  667/)
    jkdata(1:13,746) = (/1, 1, 5, 13,  1, 63, 119, 113, 341,  209, 1861, 3633,  513/)
    jkdata(1:13,747) = (/1, 1, 3,  1,  9, 47,  51, 253, 227,  875, 1979, 2367, 2303/)
    jkdata(1:13,748) = (/1, 1, 5,  1,  7, 57, 125,  99, 375,  639, 1569, 1261, 4591/)
    jkdata(1:13,749) = (/1, 3, 5,  5, 29, 61,  63,  17,  61,    7, 1087, 3953, 7941/)
    jkdata(1:13,750) = (/1, 3, 7,  1, 27, 49,  13, 119, 331,  595, 1009, 1735, 2741/)
    jkdata(1:13,751) = (/1, 3, 5,  3, 21,  9,  15, 105, 493,  971,  165,  171,  987/)
    jkdata(1:13,752) = (/1, 1, 3,  1, 23, 59,  45, 117, 411,  263, 1895, 1959, 8061/)
    jkdata(1:13,753) = (/1, 3, 5,  7, 13, 19,  61, 129, 293, 1009, 1481, 2867, 3161/)
    jkdata(1:13,754) = (/1, 3, 5,  1, 25, 29,  19, 243,  47,  201, 1583,  859, 5951/)
    jkdata(1:13,755) = (/1, 1, 5,  1, 29, 21, 105,  75, 203,   23,   29, 2951, 1431/)
    jkdata(1:13,756) = (/1, 3, 1,  5, 15, 23, 115, 203, 375,   77, 1193, 3211,  831/)
    jkdata(1:13,757) = (/1, 1, 5,  1, 17, 55,  17,  53, 167,  621, 1673,   15, 5559/)
    jkdata(1:13,758) = (/1, 1, 5, 11, 29, 23,  83,  29, 395,   33, 1075, 1279, 7405/)
    jkdata(1:13,759) = (/1, 3, 5, 11,  9, 43,   7, 247, 155,  535,  301, 1323, 1357/)
    jkdata(1:13,760) = (/1, 3, 5,  9, 15, 41,   7,  35,   5,  963, 1081,  599, 4319/)
    jkdata(1:13,761) = (/1, 3, 1,  9,  3, 57,  11, 247, 237,  661, 1377, 1651, 4235/)
    jkdata(1:13,762) = (/1, 1, 3,  5, 21,  9,  61, 171, 361,  523, 1747, 3951, 5421/)
    jkdata(1:13,763) = (/1, 3, 5, 13, 15, 39,  37,  31, 489,  263, 1497, 1011, 2559/)
    jkdata(1:13,764) = (/1, 3, 3,  7, 17, 27,  63, 199, 127,  917, 1103,  315, 4415/)
    jkdata(1:13,765) = (/1, 1, 1,  7, 17, 41,  89, 213,  21,  103, 1789, 3513, 2439/)
    jkdata(1:13,766) = (/1, 1, 7,  7, 31, 35,  95,  29, 345,  623,  887, 3351,  823/)
    jkdata(1:13,767) = (/1, 1, 5,  1,  9, 61, 119, 251, 101,  231,  739, 1725, 1725/)
    jkdata(1:13,768) = (/1, 3, 5,  1,  9, 29, 113,   7, 371,   47, 1577, 3793, 6219/)
    jkdata(1:13,769) = (/1, 1, 7,  9, 23, 57,  67, 251, 233,  301,  313, 2399, 4903/)
    jkdata(1:13,770) = (/1, 3, 1,  9, 19, 63, 123, 187, 431,  549, 1367,  287, 6699/)
    jkdata(1:13,771) = (/1, 3, 5, 11, 25, 21,  91,  91, 109,  337, 1299, 4017, 5451/)
    jkdata(1:13,772) = (/1, 3, 3, 11,  3, 31,  33,  11, 119,  675, 1801, 3571,  349/)
    jkdata(1:13,773) = (/1, 3, 3, 15,  1, 59,  37, 149, 277,  189, 1131, 1007, 7703/)
    jkdata(1:13,774) = (/1, 3, 1,  7, 11, 35,  99,  13, 125,  357, 1837,  541, 2927/)
    jkdata(1:13,775) = (/1, 3, 5,  5, 27, 49,  43, 205, 263, 1005,   73, 3115, 7809/)
    jkdata(1:13,776) = (/1, 3, 3,  5, 29,  3,  11,  37,  73,  789, 1865,  429, 6179/)
    jkdata(1:13,777) = (/1, 3, 7,  3,  1, 49,  33, 249, 135,  189, 1065, 1585, 1417/)
    jkdata(1:13,778) = (/1, 1, 1, 11, 31, 47,  65, 137, 123,  319,  843, 1285, 5987/)
    jkdata(1:13,779) = (/1, 3, 7,  1, 29, 49,  81, 139,  83,  721,  635,  755, 3017/)
    jkdata(1:13,780) = (/1, 3, 5,  3, 25, 33,  79,   9, 123, 1005,   55, 1211, 4983/)
    jkdata(1:13,781) = (/1, 1, 1,  7, 29, 21,  81,   7, 405,  525, 1655, 3047, 3479/)
    jkdata(1:13,782) = (/1, 3, 1, 13,  1, 19, 107, 113,  69,  675,  913,  915, 4525/)
    jkdata(1:13,783) = (/1, 1, 3,  7, 23, 21,  63, 183,  75,  539, 1037, 3611, 4643/)
    jkdata(1:13,784) = (/1, 1, 1,  7, 29, 35,  63, 205, 287,  191,  223, 2697, 4911/)
    jkdata(1:13,785) = (/1, 3, 1,  7, 25, 11,  55, 187, 401,  813, 1871, 2129,  227/)
    jkdata(1:13,786) = (/1, 3, 7,  3, 13, 17,  89,  39,  23,  917, 1161, 3669, 5475/)
    jkdata(1:13,787) = (/1, 3, 1, 15,  3, 37,  91,   3, 283,   51,  461,   81, 2287/)
    jkdata(1:13,788) = (/1, 1, 5, 15, 31, 23,  25,  79, 393,  167,  479, 3939, 5581/)
    jkdata(1:13,789) = (/1, 3, 5, 11, 25, 59,  93, 155,  41,  415,  511, 2437, 6817/)
    jkdata(1:13,790) = (/1, 3, 3,  9,  5, 13, 101, 227, 379,  579, 1721,  915, 1937/)
    jkdata(1:13,791) = (/1, 3, 7,  3,  5, 37,  27,  89, 431,  755, 1107,  779, 1421/)
    jkdata(1:13,792) = (/1, 3, 3,  9, 11, 35,  55, 185,  11,  605,  389, 3567, 4415/)
    jkdata(1:13,793) = (/1, 3, 7,  3,  3, 55,  75,  51, 475,  721,  151, 3701, 7977/)
    jkdata(1:13,794) = (/1, 1, 5, 15, 21, 57, 121, 127, 505,  837,   35, 2479, 1789/)
    jkdata(1:13,795) = (/1, 3, 3, 13,  9,  1,  79,  63,  19,  529,  375, 3807, 3907/)
    jkdata(1:13,796) = (/1, 3, 1,  5, 23, 29,  43,  83, 365,   31, 1099, 1893, 6815/)
    jkdata(1:13,797) = (/1, 3, 1,  3,  7, 45, 125,  41, 265,  327,  937, 3927, 6789/)
    jkdata(1:13,798) = (/1, 1, 3,  3, 11, 11,  73, 133, 271,  799, 1185, 2619, 6003/)
    jkdata(1:13,799) = (/1, 1, 1,  3, 23,  1,  27, 183, 499,  961, 1701, 2543, 5609/)
    jkdata(1:13,800) = (/1, 1, 3,  5, 11, 15, 109, 181, 489,  279,  769, 3633, 4507/)
    jkdata(1:13,801) = (/1, 3, 5,  9,  1,  9,  35, 127, 443,  409,  639, 2007,  337/)
    jkdata(1:13,802) = (/1, 3, 5, 15,  1, 33,  21,  19, 165,  847, 1633, 3857, 7427/)
    jkdata(1:13,803) = (/1, 1, 7,  9,  3, 19,  71, 255,  91,  649, 1609, 3837, 7943/)
    jkdata(1:13,804) = (/1, 3, 5,  9, 23, 53, 113, 219,  83,  241,  379,  487, 3075/)
    jkdata(1:13,805) = (/1, 3, 3,  1, 25, 43,  89,  59, 291,  285, 1613, 1769, 6427/)
    jkdata(1:13,806) = (/1, 1, 7,  5, 23, 39,  59, 251, 319,  545, 2031, 3759, 1019/)
    jkdata(1:13,807) = (/1, 3, 7,  9,  1, 23,  95,   3, 199,  407,  685, 3105, 7121/)
    jkdata(1:13,808) = (/1, 1, 7,  9, 23,  7,  41, 187, 107,  161,  289, 2727, 4763/)
    jkdata(1:13,809) = (/1, 3, 3, 15,  3, 13,  45,  57, 245,  591,  975, 3155,   81/)
    jkdata(1:13,810) = (/1, 1, 7,  5, 27, 13, 113, 217, 389,   73,  671, 2479, 3587/)
    jkdata(1:13,811) = (/1, 3, 3, 15,  9,  1, 119, 115, 143,  313, 1599, 1341, 2929/)
    jkdata(1:13,812) = (/1, 1, 7,  7, 27, 19, 113, 217, 137,  811, 1447, 1657, 1795/)
    jkdata(1:13,813) = (/1, 3, 1,  9,  3, 41,  39, 229,  89,   17,  871, 2767, 8067/)
    jkdata(1:13,814) = (/1, 3, 3,  1, 23, 55,  59, 181, 125,  663,  647, 2541, 2415/)
    jkdata(1:13,815) = (/1, 3, 1,  9, 25,  1,  73, 185, 281,  269,   99,  577, 1265/)
    jkdata(1:13,816) = (/1, 3, 7,  9, 19, 13,  15, 149, 381,  261,  139, 2105, 4025/)
    jkdata(1:13,817) = (/1, 3, 7,  5, 29, 15,  13,  83, 215,   37, 1427,  799, 5599/)
    jkdata(1:13,818) = (/1, 3, 1, 11, 29, 59,  59, 115, 131,  783,  959,   17, 4771/)
    jkdata(1:13,819) = (/1, 1, 7,  5, 13, 55,  67,  11, 299,  127,   89, 2871, 3025/)
    jkdata(1:13,820) = (/1, 1, 3, 15, 27, 15, 121, 123, 249,  917,  117, 3637, 2313/)
    jkdata(1:13,821) = (/1, 3, 7, 15,  5,  3,  27,  19, 375,  231,  841,  953, 6129/)
    jkdata(1:13,822) = (/1, 1, 3, 11,  9, 57,   7, 109, 455,  577,  891,   65, 7611/)
    jkdata(1:13,823) = (/1, 3, 7,  7, 29, 37, 105, 165,  43,  975, 1959,   69, 6881/)
    jkdata(1:13,824) = (/1, 1, 3,  7, 29, 31,  15, 103,  73,  793,  223, 2897, 5253/)
    jkdata(1:13,825) = (/1, 1, 7,  7, 13, 17,  59, 123, 281,  921, 1697, 3841, 4413/)
    jkdata(1:13,826) = (/1, 1, 3,  1, 17,  1,  59, 219, 217,  343, 1145, 3559, 7869/)
    jkdata(1:13,827) = (/1, 1, 5,  1,  3,  3,  35, 129, 297,  751,  499, 4067,  105/)
    jkdata(1:13,828) = (/1, 1, 1, 11, 23, 21,  91, 155, 229,  139, 1435, 2335, 3173/)
    jkdata(1:13,829) = (/1, 3, 1, 11, 19, 29,  89, 207, 431,  221, 1809, 3409, 1629/)
    jkdata(1:13,830) = (/1, 1, 7, 13,  7, 25,  23, 177, 357,   79, 1413, 1087, 2537/)
    jkdata(1:13,831) = (/1, 1, 3, 15, 13, 55, 125,   9,  81,  817, 1445,  425, 1023/)
    jkdata(1:13,832) = (/1, 1, 1,  3,  3,  9,  97,  49, 357,  393, 1675, 2813, 4409/)
    jkdata(1:13,833) = (/1, 3, 5, 13, 19, 37,  53, 181, 171,  545,  171, 1705, 7209/)
    jkdata(1:13,834) = (/1, 1, 5,  5, 23, 33,  41, 231, 451,   11, 1073, 1701, 4413/)
    jkdata(1:13,835) = (/1, 3, 7,  1,  5, 53,  91,  33, 481,  781, 1349, 1237, 7107/)
    jkdata(1:13,836) = (/1, 1, 1,  7, 29, 41, 111, 233,  13,   71, 1545,  821, 7469/)
    jkdata(1:13,837) = (/1, 1, 5,  1, 29, 51,  29,  67, 387,    1, 2039, 1375,   33/)
    jkdata(1:13,838) = (/1, 3, 5, 11, 13, 19,  31, 155, 491,  699, 1027, 3673, 1955/)
    jkdata(1:13,839) = (/1, 3, 5,  3, 13, 57,   3,  41, 489,  767, 1563, 2693, 2881/)
    jkdata(1:13,840) = (/1, 3, 7, 13,  5, 13, 103,   9, 439,  917,  859, 3925, 5167/)
    jkdata(1:13,841) = (/1, 1, 1, 15, 19, 63,  61,  95, 385,    9,  215, 1541, 6451/)
    jkdata(1:13,842) = (/1, 3, 5,  3,  5, 43,  71, 123, 487,  107, 1673, 1871, 4211/)
    jkdata(1:13,843) = (/1, 1, 5,  5, 17, 19,  35,  65, 177,  341, 1919, 2285,  179/)
    jkdata(1:13,844) = (/1, 3, 1,  3,  9,  7,   7, 117, 393,  587, 1633,  847, 5573/)
    jkdata(1:13,845) = (/1, 1, 5,  5, 11, 13, 119, 249,  33,  903,  779, 4035, 7879/)
    jkdata(1:13,846) = (/1, 1, 5,  7, 11, 37,  29,  85,  71,  965,  411, 1101, 3387/)
    jkdata(1:13,847) = (/1, 3, 3,  3, 29, 33,  45, 169, 375,  599, 1845, 2029, 7759/)
    jkdata(1:13,848) = (/1, 1, 1,  9, 27, 19,  49, 129, 443,  507, 1477,  855, 5455/)
    jkdata(1:13,849) = (/1, 3, 3,  9, 23, 15, 111, 241, 129,  843, 1489, 2733, 7157/)
    jkdata(1:13,850) = (/1, 3, 1,  5, 19, 63,  41, 173, 407,  739,  447, 2503, 1891/)
    jkdata(1:13,851) = (/1, 1, 7,  1, 17, 51, 109, 251, 395,  579, 1545,  121, 5683/)
    jkdata(1:13,852) = (/1, 3, 3,  7, 25, 11,  59, 225, 127,  397,  351, 2855, 5689/)
    jkdata(1:13,853) = (/1, 1, 1, 11, 13, 49, 125, 147,  65,  397, 1989, 1069, 6535/)
    jkdata(1:13,854) = (/1, 3, 3,  9,  1, 23,  13, 165, 333,  325,  495, 3463, 3109/)
    jkdata(1:13,855) = (/1, 3, 5,  3, 13, 57,  27,  69, 309,  775,  183, 3505, 6555/)
    jkdata(1:13,856) = (/1, 1, 7,  5,  3, 47,  19,  81, 119,  565, 1639, 1539, 6873/)
    jkdata(1:13,857) = (/1, 3, 7, 11, 11, 51,  79, 239, 197,  925, 1385,  607, 1249/)
    jkdata(1:13,858) = (/1, 3, 7, 13,  1, 15,   9,  95, 435,   75, 1805, 1349, 4251/)
    jkdata(1:13,859) = (/1, 1, 1, 13, 17, 53,  75,  23, 497,   55, 1097,  575, 6437/)
    jkdata(1:13,860) = (/1, 3, 1, 13, 29, 41,  83,  83, 373,  979, 1249, 2301,   49/)
    jkdata(1:13,861) = (/1, 3, 7,  9,  1,  1,  81, 227,  71,  931, 1431, 2321, 2745/)
    jkdata(1:13,862) = (/1, 3, 3, 15, 13, 15,  33, 249, 379,   93, 1571, 1101, 1201/)
    jkdata(1:13,863) = (/1, 3, 1,  5, 17, 37,  91, 143, 509,  957,  591,  333, 7327/)
    jkdata(1:13,864) = (/1, 3, 5,  7,  9, 61, 109, 171, 387,  857,  697,  291, 4179/)
    jkdata(1:13,865) = (/1, 3, 5,  1, 17, 11,  33, 193, 159,  753, 1509, 2171, 6783/)
    jkdata(1:13,866) = (/1, 1, 5, 15, 21, 35,  29,   9, 265,  965,  709, 4085,  623/)
    jkdata(1:13,867) = (/1, 3, 1, 11,  1, 29, 107,  21, 477,  795,   31, 2173, 2779/)
    jkdata(1:13,868) = (/1, 1, 1,  9, 11, 33, 111,  57, 463,   67, 1563, 2541, 5963/)
    jkdata(1:13,869) = (/1, 1, 1, 15,  1, 23, 101,  73, 449,    5,  165, 1195, 2585/)
    jkdata(1:13,870) = (/1, 3, 1, 15,  1, 55, 107,  97,  47,   87,  513,  925, 6927/)
    jkdata(1:13,871) = (/1, 3, 1, 13, 25, 11, 109,  57, 353,  909, 1425, 4039, 5333/)
    jkdata(1:13,872) = (/1, 3, 5, 13,  5, 59,  65,  29, 249,   97, 1299, 1379, 4033/)
    jkdata(1:13,873) = (/1, 1, 3, 13,  7, 19,  59, 239, 335,  995, 1081,  699,  285/)
    jkdata(1:13,874) = (/1, 1, 5,  1, 29, 61,  43, 151, 505,  271,  145, 1979, 7467/)
    jkdata(1:13,875) = (/1, 3, 1, 11, 29, 61,  37, 159,  89,  875, 1841,  275, 4443/)
    jkdata(1:13,876) = (/1, 3, 3,  9, 19, 45,   1, 191, 141,  671, 1211,  953, 4917/)
    jkdata(1:13,877) = (/1, 3, 5, 15, 19, 13,   9,  47,  55,  613,  941, 1755,    3/)
    jkdata(1:13,878) = (/1, 3, 3,  9,  1, 49,  15,  51, 235,   33,  609, 1643, 4319/)
    jkdata(1:13,879) = (/1, 3, 1,  5, 29, 13, 109,   1, 187,  351,  845,  325, 5517/)
    jkdata(1:13,880) = (/1, 3, 1, 15, 13, 63,  37, 223,  87,   69, 1169,  101, 3449/)
    jkdata(1:13,881) = (/1, 3, 1,  5,  3,  5, 111, 251, 363,  811, 1865, 2263,  813/)
    jkdata(1:13,882) = (/1, 1, 1,  7,  1, 61, 113, 251,  93,  669, 1593, 3329, 5499/)
    jkdata(1:13,883) = (/1, 3, 3,  3, 31,  5, 119, 151, 363,  729,  347, 3673, 2515/)
    jkdata(1:13,884) = (/1, 3, 7, 11, 15, 31,  79,  41, 101,  401,  293, 3413, 5771/)
    jkdata(1:13,885) = (/1, 3, 3,  3, 13, 17,  73, 119,  67,  647, 1277, 1977, 3357/)
    jkdata(1:13,886) = (/1, 3, 7, 15,  3, 61,  65, 127, 215,  241,  157, 2727, 2073/)
    jkdata(1:13,887) = (/1, 1, 5,  7,  1, 63,  71, 131, 321,  435,  211, 2313, 4395/)
    jkdata(1:13,888) = (/1, 3, 7, 13, 11, 13,  93,  33, 331,  447,   93, 1419, 4925/)
    jkdata(1:13,889) = (/1, 1, 1, 11, 19, 27,  17, 209, 305,  721, 1679,  887, 2643/)
    jkdata(1:13,890) = (/1, 3, 5,  7,  5, 57, 101, 123, 261,  271, 1799,  609, 7215/)
    jkdata(1:13,891) = (/1, 3, 5,  3, 29,  1,  87,  53, 411,  745,  527, 2475, 5817/)
    jkdata(1:13,892) = (/1, 3, 7,  7, 13, 21,  97, 241, 491,   53,   41,  591, 1199/)
    jkdata(1:13,893) = (/1, 1, 5, 13, 29,  5,  43,  25, 479,  775,  473, 2613, 1597/)
    jkdata(1:13,894) = (/1, 3, 3,  5, 23, 11,  23,  31,  65,   99,  563, 2081, 1619/)
    jkdata(1:13,895) = (/1, 1, 3, 13,  3, 39,  75, 183, 307,  343,  187, 3805, 7535/)
    jkdata(1:13,896) = (/1, 3, 7, 15,  1, 57, 109, 107, 469,  451, 1525, 3435, 4833/)
    jkdata(1:13,897) = (/1, 1, 5,  5, 31, 51,  41,  25, 415,  427,  575, 2409,  609/)
    jkdata(1:13,898) = (/1, 1, 3, 13, 13, 53,  49, 115, 131,  593, 1579,  111, 4797/)
    jkdata(1:13,899) = (/1, 1, 1,  9, 19, 39,  53,  39, 315,  339,  857, 3557, 8171/)
    jkdata(1:13,900) = (/1, 3, 1,  1, 17, 25,  31,  11, 487,  845,  703, 3607, 6847/)
    jkdata(1:13,901) = (/1, 3, 3, 15,  5, 41,  97, 213,  83,  243, 1211,  903,  793/)
    jkdata(1:13,902) = (/1, 1, 1, 11,  5, 39, 105, 239, 455,  345,  647,  231, 6757/)
    jkdata(1:13,903) = (/1, 3, 3,  5,  1, 37, 109, 219,  19,   17,  709, 3059, 8165/)
    jkdata(1:13,904) = (/1, 1, 1,  5, 29, 23, 119, 109, 113,  573,  981,  473, 3371/)
    jkdata(1:13,905) = (/1, 1, 1,  1, 23, 31,  51, 185, 163,  421,  285, 2959, 2431/)
    jkdata(1:13,906) = (/1, 3, 3, 11,  3, 25,   9,  35, 503,  517,  697, 2925, 5235/)
    jkdata(1:13,907) = (/1, 3, 7,  3, 19, 33,  53, 133,  99,  971,  163, 3861, 4739/)
    jkdata(1:13,908) = (/1, 1, 1,  3, 25, 17, 113, 123, 499,  499,  981, 2043, 7703/)
    jkdata(1:13,909) = (/1, 3, 7,  7, 19, 57,  97, 185, 251,  435,  153, 3887, 7223/)
    jkdata(1:13,910) = (/1, 1, 1,  1, 27, 29,  73,  27, 239,  769, 1515,  351, 6525/)
    jkdata(1:13,911) = (/1, 1, 1,  9,  9, 27,  89,  55,  81,   75,   47, 2865, 5891/)
    jkdata(1:13,912) = (/1, 1, 5,  7, 27, 23,  79, 245, 167,  203, 1553,  369, 5605/)
    jkdata(1:13,913) = (/1, 1, 1, 15, 13, 47,  49,  61, 391,  793,  599, 1377, 4433/)
    jkdata(1:13,914) = (/1, 3, 7,  9, 15, 41,  61,  75, 255,  985,  225, 2639, 3533/)
    jkdata(1:13,915) = (/1, 1, 5,  9, 29, 29, 105, 205, 317,  343, 1147, 1261, 5267/)
    jkdata(1:13,916) = (/1, 3, 3,  3, 23, 19,  13, 213, 363,  955,  381, 3625, 5125/)
    jkdata(1:13,917) = (/1, 1, 7, 11, 13, 47,  99, 169, 359,  735,  135, 3279, 5037/)
    jkdata(1:13,918) = (/1, 1, 3, 15, 25, 41,  53, 163, 395,  523,  821, 2201,  225/)
    jkdata(1:13,919) = (/1, 3, 5,  7, 25, 25,  71,  63, 419,  659, 1965, 2949, 6717/)
    jkdata(1:13,920) = (/1, 1, 3,  1, 17,  5,   7,  55, 307,  703,  609, 3049, 1121/)
    jkdata(1:13,921) = (/1, 3, 1,  3, 19, 51,  87,  49, 251,  303, 1033,  449, 5741/)
    jkdata(1:13,922) = (/1, 1, 1,  1, 17, 43,  21,  83, 267,  421,  983, 1297, 2013/)
    jkdata(1:13,923) = (/1, 3, 5,  1, 15, 39, 101, 195, 171,  951,  503,  897, 4327/)
    jkdata(1:13,924) = (/1, 3, 5,  1, 27, 29,   5,  51, 461,  405, 1117, 1891, 4839/)
    jkdata(1:13,925) = (/1, 3, 1,  9,  3,  7,  71,  31, 183,  631,  327,  411,  569/)
    jkdata(1:13,926) = (/1, 3, 7,  1, 25, 31,  31,  41, 465,  825,  453, 2773, 5227/)
    jkdata(1:13,927) = (/1, 3, 7,  5, 17, 45, 123,  15, 165,  735, 2005,  749, 7677/)
    jkdata(1:13,928) = (/1, 3, 3, 15, 27, 51, 121, 203, 163,  433, 1257, 2753, 4315/)
    jkdata(1:13,929) = (/1, 1, 7, 15,  3, 49, 121,  41, 293,  841,  343, 1825, 2391/)
    jkdata(1:13,930) = (/1, 3, 3,  7, 27, 55,  73,  63, 477,  485, 1649,  853, 5551/)
    jkdata(1:13,931) = (/1, 3, 7,  5, 31, 17,  79, 127, 223,   49, 1199, 2775,  859/)
    jkdata(1:13,932) = (/1, 3, 1,  5, 23, 43, 115, 161, 403,  749,  599, 3547, 3627/)
    jkdata(1:13,933) = (/1, 3, 5,  7, 13, 49,  13,   5, 389,  107, 1877, 3923, 6377/)
    jkdata(1:13,934) = (/1, 1, 1,  9, 31, 45,  39, 143,  97,  669,  569, 3923, 3903/)
    jkdata(1:13,935) = (/1, 3, 5,  7, 11,  9, 101,   7, 335,  211,  695,  987, 4311/)
    jkdata(1:13,936) = (/1, 3, 3, 15, 15, 29,  19, 199, 357,  497, 1587, 3723, 6527/)
    jkdata(1:13,937) = (/1, 1, 7, 13,  7,  3,  37, 251, 297,  143, 1475, 2189, 7573/)
    jkdata(1:13,938) = (/1, 3, 3, 13, 21,  5,  51,  95,  19,   99,  187, 3877, 4905/)
    jkdata(1:13,939) = (/1, 3, 5, 11, 19, 47,  83,  75, 469,   57,  973, 3577, 7731/)
    jkdata(1:13,940) = (/1, 3, 7,  1, 27,  9,  97, 101, 501,  277,  233,  297, 1909/)
    jkdata(1:13,941) = (/1, 3, 7,  9, 19, 15,  55,  15, 249,  969,  511, 2763, 1555/)
    jkdata(1:13,942) = (/1, 3, 7, 11, 21, 19,  81,  43,  85,  107,   51, 1845, 3279/)
    jkdata(1:13,943) = (/1, 1, 3,  1, 29, 51,  91, 237, 213,  397, 1083, 3083, 1949/)
    jkdata(1:13,944) = (/1, 1, 3, 13,  7, 45, 127, 197, 311,  563,  665, 2951, 1887/)
    jkdata(1:13,945) = (/1, 1, 1,  1, 31, 57, 105, 117, 265,  551, 1321,  483, 6675/)
    jkdata(1:13,946) = (/1, 1, 1,  7, 13, 63,  89, 167, 379,  447,  531, 2169, 5509/)
    jkdata(1:13,947) = (/1, 3, 5, 15,  9,  9,  63, 155, 297,  381, 1875, 3985, 2033/)
    jkdata(1:13,948) = (/1, 3, 5, 15,  9, 21,  47,  21, 283,  187, 1939,  245, 5473/)
    jkdata(1:13,949) = (/1, 3, 3,  5,  7, 59,  49,  83, 393,   57,  859, 3655, 3539/)
    jkdata(1:13,950) = (/1, 1, 7,  5, 21,  3,  75, 205, 449,  405, 1507, 3441, 5033/)
    jkdata(1:13,951) = (/1, 3, 1,  1, 13,  9,  37, 255, 463,  731, 1979, 1023, 5935/)
    jkdata(1:13,952) = (/1, 3, 1, 11, 11, 13,  77,  49, 289,  769, 1203,  235, 6095/)
    jkdata(1:13,953) = (/1, 1, 1,  3,  9, 45,  15, 101, 159,  923, 1965,  835, 4761/)
    jkdata(1:13,954) = (/1, 1, 3,  9, 11, 23,  49, 213, 289,  955,  737, 3693, 1771/)
    jkdata(1:13,955) = (/1, 3, 5, 11, 29, 15, 107, 237, 499,  915,  921, 3585, 1271/)
    jkdata(1:13,956) = (/1, 3, 3,  9, 19, 31,  23, 135, 407,  737, 1565,  327, 1717/)
    jkdata(1:13,957) = (/1, 1, 1,  9, 11, 21,  23, 135, 129,  595, 1943, 1003, 4415/)
    jkdata(1:13,958) = (/1, 1, 1,  9, 19, 15,  35,  21, 137,  341,  819,  543, 5083/)
    jkdata(1:13,959) = (/1, 3, 3,  1, 21, 51,  19,  73, 221,  253,  223, 3059, 6277/)
    jkdata(1:13,960) = (/1, 3, 3,  9,  5, 35,  69,  93,  43,  823,  365, 2637, 3147/)
    jkdata(1:13,961) = (/1, 1, 7,  3, 29,  9,  17, 115,  89,  197,  167, 2923, 7695/)
    jkdata(1:13,962) = (/1, 3, 5,  5, 13, 11,  59,   7, 403,  321, 1705,   87, 2461/)
    jkdata(1:13,963) = (/1, 1, 1, 15,  7, 61,  63,  85, 271,  315,  413, 3617, 4783/)
    jkdata(1:13,964) = (/1, 1, 1,  1, 19, 23,  73, 223,  75,  181, 1577, 1031, 4539/)
    jkdata(1:13,965) = (/1, 3, 3,  1, 19, 53,  29, 237,  83,  885,  745, 1043, 5833/)
    jkdata(1:13,966) = (/1, 1, 7,  9, 27, 29, 125,  79, 445,  497, 1573,  903, 5583/)
    jkdata(1:13,967) = (/1, 3, 1,  7, 23, 51,  61,  89, 453,  159,  655, 2913,  651/)
    jkdata(1:13,968) = (/1, 3, 5,  3, 31, 45,  65,   5, 389,  571, 1633, 2177, 1419/)
    jkdata(1:13,969) = (/1, 3, 7,  3,  1, 31,  95,  57, 149,  981, 1003, 2641, 2605/)
    jkdata(1:13,970) = (/1, 3, 3,  1, 27, 29, 101, 239, 143,  899,   91, 3279, 5511/)
    jkdata(1:13,971) = (/1, 3, 7,  9, 21,  5,  81,  67, 423,  785, 1123,  389, 3913/)
    jkdata(1:13,972) = (/1, 1, 5,  9,  7, 35,  57,  65, 499,  947,  477, 2009, 5795/)
    jkdata(1:13,973) = (/1, 3, 5, 11,  3, 29,  69, 201, 317,  217, 1741,  525, 2333/)
    jkdata(1:13,974) = (/1, 1, 7,  9,  7, 53,  83, 155, 445,  217, 1663, 4085, 2329/)
    jkdata(1:13,975) = (/1, 1, 3,  9, 11, 35,  37,  71, 157,  135,   35, 3299, 4431/)
    jkdata(1:13,976) = (/1, 3, 5, 13, 23, 17,  11,  85, 137,  753,  715,  987, 3725/)
    jkdata(1:13,977) = (/1, 3, 3, 13, 13, 59,  37, 195, 453,  623,   37, 2409, 6069/)
    jkdata(1:13,978) = (/1, 3, 1,  3, 29, 55,  95,  89, 163,  565, 1513,  813, 2699/)
    jkdata(1:13,979) = (/1, 3, 5, 13, 11, 27,   1, 181,  87,  717,  815, 2683, 7055/)
    jkdata(1:13,980) = (/1, 1, 3, 11, 31, 51,  73, 119,  23,  903,  941,  373, 6879/)
    jkdata(1:13,981) = (/1, 3, 1, 13, 19, 59,  27, 135, 391,  581, 1379, 2695, 1017/)
    jkdata(1:13,982) = (/1, 1, 1,  5,  1, 27,  29, 147, 119,  955,  263, 3775, 3121/)
    jkdata(1:13,983) = (/1, 1, 7,  1,  5, 47,  57, 237, 427,  621, 1831, 2375, 2547/)
    jkdata(1:13,984) = (/1, 3, 5,  5,  5, 15,   7, 173, 323,  361, 1735, 1119, 4603/)
    jkdata(1:13,985) = (/1, 3, 1,  5, 11, 29,  65,  41, 173,  869, 1111, 2791, 2385/)
    jkdata(1:13,986) = (/1, 3, 7,  9,  5, 37,  83, 155,  89,   87, 1449,  223, 6915/)
    jkdata(1:13,987) = (/1, 3, 3,  9,  3,  7,  99,  67, 259,  943,  353,  325, 6103/)
    jkdata(1:13,988) = (/1, 3, 7,  3, 27, 49,  69, 113, 377,  907, 1941,  587, 5669/)
    jkdata(1:13,989) = (/1, 3, 5, 13,  5, 55,  19, 111, 511,  853, 1655, 1379, 7833/)
    jkdata(1:13,990) = (/1, 1, 1, 13,  7,  5, 103,  21, 249,  353, 1349, 2877, 2001/)
    jkdata(1:13,991) = (/1, 1, 7,  9, 11, 19,  43, 183,  31,  335,  877, 2867, 4287/)
    jkdata(1:13,992) = (/1, 3, 1, 15, 31, 45,  95,  23, 363,  197,  285, 3793, 6619/)
    jkdata(1:13,993) = (/1, 1, 7,  9,  1, 29,  25, 103, 229,  771, 1723,  655,  955/)
    jkdata(1:13,994) = (/1, 3, 7, 11, 27, 19,  19, 207, 353,  433,  125,  831, 2761/)
    jkdata(1:13,995) = (/1, 1, 1,  7, 31, 57, 103, 253, 329,  743, 1753, 3425, 5711/)
    jkdata(1:13,996) = (/1, 1, 1, 11, 31, 33,  41,  69, 493,  195,  985, 1663, 6291/)
    jkdata(1:13,997) = (/1, 3, 7,  9, 23, 53, 125, 219, 427,   91,  723, 1681, 3415/)
    jkdata(1:13,998) = (/1, 1, 1, 13,  5, 45,  97, 205,  57, 1023,  175, 2657, 3909/)
    jkdata(1:13,999) = (/1, 1, 5,  9, 21, 21,  71, 195, 205,   63,  439, 1865, 2841/)
    jkdata(1:13,1000) = (/1, 1, 5,  1, 27,  9, 105,  43, 389,  301,  791, 3943, 5627/)
    jkdata(1:13,1001) = (/1, 1, 1, 15,  9,  3,  83, 197,  91,  647, 1051, 2977, 4939/)
    jkdata(1:13,1002) = (/1, 3, 1,  9, 25, 35,  83, 229,  83,  205, 1261, 1979, 7671/)
    jkdata(1:13,1003) = (/1, 3, 7,  7,  3, 29,  61, 139,  13,  485,  717, 2271, 6059/)
    jkdata(1:13,1004) = (/1, 1, 5,  7, 15, 43,  39, 177, 219,  927, 1555, 3247, 6275/)
    jkdata(1:13,1005) = (/1, 1, 7,  1, 19, 31,   9, 129, 439, 1003, 1757, 1267, 6517/)
    jkdata(1:13,1006) = (/1, 3, 1,  7,  1, 39,  45,  69,  45,  987, 1777, 1747, 1931/)
    jkdata(1:13,1007) = (/1, 1, 5,  9, 19,  3, 117,  97,  35,  359,  577,  811, 4583/)
    jkdata(1:13,1008) = (/1, 1, 3,  9,  9, 45,  63, 201, 371,  577, 1583,  159, 7301/)
    jkdata(1:13,1009) = (/1, 1, 5, 15,  5,  1,  31, 163, 441,  147, 1957,  429, 1267/)
    jkdata(1:13,1010) = (/1, 3, 3,  1, 25, 41,   5, 189,  17,  141,  873, 2001, 7509/)
    jkdata(1:13,1011) = (/1, 1, 3, 11, 21, 29, 117,  11, 267, 1017,  331, 1195, 1435/)
    jkdata(1:13,1012) = (/1, 3, 7,  1, 15,  5,  67,  99, 501,  701, 1163, 3065, 2169/)
    jkdata(1:13,1013) = (/1, 1, 1, 13, 25, 59, 125,  91,  53,  273,  313,  553, 6939/)
    jkdata(1:13,1014) = (/1, 1, 5, 13, 29, 41,  41, 253,  25,   89,    1, 1499, 3515/)
    jkdata(1:13,1015) = (/1, 3, 1, 15, 15, 33, 117, 239, 333,  589, 1963, 3529, 2985/)
    jkdata(1:13,1016) = (/1, 3, 1,  9, 21, 35,  43,  91,  17,  487,  963, 1081, 2787/)
    jkdata(1:13,1017) = (/1, 1, 5, 13, 11, 27,  77, 145, 201,  859, 1905, 2877, 2123/)
    jkdata(1:13,1018) = (/1, 3, 5,  7, 19, 19,  97,  19, 475,  343,  821, 3077, 1969/)
    jkdata(1:13,1019) = (/1, 1, 3, 15, 15, 13,  15, 179, 257,   91, 1677,  845, 3307/)
    jkdata(1:13,1020) = (/1, 1, 3,  3,  3, 25,  29, 231, 417,  847,  185, 1793,  353/)
    jkdata(1:13,1021) = (/1, 3, 7,  9,  7, 27,   5, 121, 345,  341,  709, 2409, 4359/)
    jkdata(1:13,1022) = (/1, 3, 5,  3, 13, 43,  59,   7, 381,  173,  545, 3995, 7059/)
    jkdata(1:13,1023) = (/1, 3, 5,  1, 11, 33,  25, 225, 377,  287, 1723, 2559, 5273/)
    jkdata(1:13,1024) = (/1, 3, 1, 13, 25, 35,  63, 237,  55, 1003,  215, 4081, 5873/)
    jkdata(1:13,1025) = (/1, 3, 1,  7, 17, 17,  87, 125, 403,  289, 1885, 1195, 6657/)
    jkdata(1:13,1026) = (/1, 1, 1,  5,  1, 17,  39, 191,  77,  639, 1249, 2955, 6765/)
    jkdata(1:13,1027) = (/1, 3, 3,  9,  5, 23,  39, 119, 389,  983,  583, 1117, 6229/)
    jkdata(1:13,1028) = (/1, 1, 1,  3, 31,  7,  77,  59, 347,  685, 1803, 1409, 3179/)
    jkdata(1:13,1029) = (/1, 1, 5,  1, 13, 35,  85, 175, 363,  697,  839,  785, 1583/)
    jkdata(1:13,1030) = (/1, 1, 7,  7, 29, 15,  37, 237, 211,   35,  885,  287, 6237/)
    jkdata(1:13,1031) = (/1, 3, 7,  1, 23, 61,  81, 131, 413,  701,  485, 1521, 2155/)
    jkdata(1:13,1032) = (/1, 1, 1,  1,  9, 61,  73,  79, 419,  645,  413, 1607,  371/)
    jkdata(1:13,1033) = (/1, 1, 7, 13,  5, 53,  89,  43,   5,  911, 1767,   85,  273/)
    jkdata(1:13,1034) = (/1, 1, 5,  3, 29,  5,  29,  45, 167,  501,  425, 3055, 7491/)
    jkdata(1:13,1035) = (/1, 3, 7,  3,  7, 15, 125, 205, 219,  705,  129, 3123, 3309/)
    jkdata(1:13,1036) = (/1, 1, 3, 11, 17, 23, 109, 199, 201,  873, 1035, 2533, 6805/)
    jkdata(1:13,1037) = (/1, 1, 7,  1, 27, 11,  21, 251, 285,  763,  329, 2329, 3015/)
    jkdata(1:13,1038) = (/1, 3, 3,  7,  7, 13,  23, 153, 425,  745, 1263, 3477, 6831/)
    jkdata(1:13,1039) = (/1, 1, 1, 13, 17, 43, 119, 207,  11,  657, 1881,  799, 7819/)
    jkdata(1:13,1040) = (/1, 3, 3, 15, 31, 55, 105,  37,  77,  559, 1779, 3683,  713/)
    jkdata(1:13,1041) = (/1, 3, 7, 15,  9, 47,  43, 179, 269,  699, 1565, 3715, 4747/)
    jkdata(1:13,1042) = (/1, 3, 3,  5, 31, 25,  93, 113, 489,  315,  359,  337, 3935/)
    jkdata(1:13,1043) = (/1, 3, 1,  7,  9, 43,  97, 255, 281,  347,  367, 3139, 4109/)
    jkdata(1:13,1044) = (/1, 3, 5, 13,  9, 15,  15, 107, 403,  429,  453, 3311, 1311/)
    jkdata(1:13,1045) = (/1, 1, 5, 13,  7, 57, 125, 217,  79,  197,  707,  431,  709/)
    jkdata(1:13,1046) = (/1, 1, 3, 15, 21, 45,  29,  61, 425,  165, 1419, 3511, 3089/)
    jkdata(1:13,1047) = (/1, 1, 5, 11,  3,  1,  51,   7, 125,  955,  831, 2299, 7059/)
    jkdata(1:13,1048) = (/1, 3, 1, 13,  3, 49,  69, 181,  81,  859, 1889,  365, 4247/)
    jkdata(1:13,1049) = (/1, 3, 3,  1,  3, 63,  37, 247, 331,  167,  887, 2941, 2989/)
    jkdata(1:13,1050) = (/1, 3, 5, 13,  9, 57,  45,  31, 437,  303, 1871, 3067, 1509/)
    jkdata(1:13,1051) = (/1, 3, 5, 13, 11, 15,  31,  13, 271,  833, 1869, 1331, 4919/)
    jkdata(1:13,1052) = (/1, 1, 5,  3, 21, 31,  75, 113, 397,  531,  747, 1081, 1841/)
    jkdata(1:13,1053) = (/1, 3, 1,  9, 11, 31, 109, 145, 299,  473,  223, 1097, 3045/)
    jkdata(1:13,1054) = (/1, 3, 1, 15, 31,  7, 119, 107, 475,  635, 1547, 2853, 3821/)
    jkdata(1:13,1055) = (/1, 3, 7, 15,  9, 53,  53, 233, 271,  641, 1799, 2299, 6929/)
    jkdata(1:13,1056) = (/1, 3, 7, 11, 25, 27,   5, 233, 249,  195,  433,  495, 4655/)
    jkdata(1:13,1057) = (/1, 1, 1, 15,  5, 15, 101,  43, 413,  589, 1441, 1745, 1333/)
    jkdata(1:13,1058) = (/1, 1, 5,  9,  1, 47, 125,  79, 233,  821,  553,  749, 6429/)
    jkdata(1:13,1059) = (/1, 3, 5, 15, 31, 23, 121,  23, 261,  205, 2021, 3819, 6649/)
    jkdata(1:13,1060) = (/1, 3, 1,  1, 13,  7,  35, 169, 495,    3, 1303,  619, 2131/)
    jkdata(1:13,1061) = (/1, 3, 3, 13, 29, 29,  29, 137, 171,  635, 1505, 1059, 5265/)
    jkdata(1:13,1062) = (/1, 1, 5, 15,  9, 53,   7, 129,  69,  371, 1735, 3559, 1051/)
    jkdata(1:13,1063) = (/1, 3, 1,  1, 29, 47,  63, 183,  27,  891, 1619,  183,  261/)
    jkdata(1:13,1064) = (/1, 1, 5,  1,  1,  9,  17,  53, 409,  249, 1065, 3743, 8057/)
    jkdata(1:13,1065) = (/1, 1, 3,  5, 11, 53,  63,  91,  21,  123, 1161,  723, 3379/)
    jkdata(1:13,1066) = (/1, 3, 5, 11, 19,  3,  13,  55, 421,   77, 2047,  949, 2179/)
    jkdata(1:13,1067) = (/1, 3, 3,  5,  7, 25,  69, 103, 367,  623,  347, 3501, 1993/)
    jkdata(1:13,1068) = (/1, 1, 3,  1, 27, 55,  15, 223,  81,  993,  867,  733, 5655/)
    jkdata(1:13,1069) = (/1, 3, 7, 11, 13, 45, 105,  87, 483,  401,  881, 2599, 3063/)
    jkdata(1:13,1070) = (/1, 3, 5, 11, 31, 63,  51, 177, 255,  525, 1447, 3983, 6381/)
    jkdata(1:13,1071) = (/1, 1, 7,  5,  7, 21, 127, 157,  15,  427,  329, 3961, 3587/)
    jkdata(1:13,1072) = (/1, 1, 3,  3, 31, 17, 105,  79, 219,   71,  781,  911, 7417/)
    jkdata(1:13,1073) = (/1, 1, 7,  9,  7, 23,   9, 213, 365,  655, 1065, 1899, 1579/)
    jkdata(1:13,1074) = (/1, 1, 3,  1, 25, 31,  57, 139, 497,  951,  219,  985, 1541/)
    jkdata(1:13,1075) = (/1, 1, 1,  3, 23, 27,  95, 183, 181,  357,  589, 2493, 2107/)
    jkdata(1:13,1076) = (/1, 3, 3,  5, 21, 27,  59, 231,  75,  851,  645, 1795, 5085/)
    jkdata(1:13,1077) = (/1, 1, 7, 13, 29, 43, 109, 205, 431,  899, 1257,  653, 2873/)
    jkdata(1:13,1078) = (/1, 1, 7,  9, 11, 63,  35, 143,  99,  535, 1833,  157, 6141/)
    jkdata(1:13,1079) = (/1, 3, 3,  7, 11, 55,  49, 129, 325,  493,  749,  433,  955/)
    jkdata(1:13,1080) = (/1, 3, 3,  7, 13, 63,  23, 243, 407,  323, 1841, 2361, 3537/)
    jkdata(1:13,1081) = (/1, 1, 1,  1, 11, 45,  33, 205, 229, 1003, 1733, 3093, 2157/)
    jkdata(1:13,1082) = (/1, 1, 1,  9, 27, 51, 107,  93, 281,  343, 1179, 3119,  841/)
    jkdata(1:13,1083) = (/1, 1, 3,  9,  1, 15,  55,  59,  63,  515, 1191, 3679, 1999/)
    jkdata(1:13,1084) = (/1, 3, 3, 15, 23, 27,  33,  15,  83,  859, 1025, 2367, 1465/)
    jkdata(1:13,1085) = (/1, 1, 3,  7, 31,  5,  57,  89, 493, 1017, 1639, 1701, 5171/)
    jkdata(1:13,1086) = (/1, 1, 3,  5, 21, 37,  79,   9,   5,    5, 1955, 1445, 5651/)
    jkdata(1:13,1087) = (/1, 3, 3,  5, 23, 43,  73,  11, 113,  423, 1423, 1321, 1535/)
    jkdata(1:13,1088) = (/1, 3, 5, 15, 21, 11,  69,  47,  15,  315, 1685, 2397, 7235/)
    jkdata(1:13,1089) = (/1, 1, 5, 13, 19, 27,  59, 133, 271, 1011, 1711, 1241, 4349/)
    jkdata(1:13,1090) = (/1, 3, 3,  9, 31,  5, 107, 227,  37,  703,  493, 3305, 1263/)
    jkdata(1:13,1091) = (/1, 3, 3,  7,  5, 27,  55,  75,  87,   41,  549, 3985, 1453/)
    jkdata(1:13,1092) = (/1, 3, 3, 13, 31, 59,  11,   9, 451,  777,  783, 2349, 1005/)
    jkdata(1:13,1093) = (/1, 3, 1,  3, 25, 21,  63,  91, 299,  163, 1653, 4067, 6893/)
    jkdata(1:13,1094) = (/1, 3, 3, 13, 25,  7,  95,  19,  83,   95,  397, 3805, 2919/)
    jkdata(1:13,1095) = (/1, 3, 5, 11, 19, 39, 103, 171, 451,  831,  895, 3073, 1947/)
    jkdata(1:13,1096) = (/1, 3, 7, 13, 17, 27,  23, 163, 311,   79,  233, 2837, 1635/)
    jkdata(1:13,1097) = (/1, 3, 7,  7, 11, 63, 125,  79, 441,  975,  759, 1567, 3963/)
    jkdata(1:13,1098) = (/1, 1, 1,  9, 25, 35,  91,   7,  47,  235, 1505, 3783,  397/)
    jkdata(1:13,1099) = (/1, 1, 5, 13,  7, 47,  31, 103, 455,  633,  677,  451,  969/)
    jkdata(1:13,1100) = (/1, 3, 7, 13, 13, 55,  91,   5,  47,  723, 1449, 2441, 4569/)
    jkdata(1:13,1101) = (/1, 3, 3, 13,  1, 17,  51, 119, 253,  297, 1573, 1181,  655/)
    jkdata(1:13,1102) = (/1, 1, 7, 15, 29, 17,  65, 155,  13,  589, 1297,  487, 6737/)
    jkdata(1:13,1103) = (/1, 1, 1,  9, 17, 17,  61,  75, 109,  317, 1821,  543, 2995/)
    jkdata(1:13,1104) = (/1, 3, 1,  5, 23,  3,  75,  11, 369,  679, 1691, 1201, 7235/)
    jkdata(1:13,1105) = (/1, 1, 3,  5, 15, 19,  69,  71, 347,  981,  791, 3735, 7713/)
    jkdata(1:13,1106) = (/1, 3, 5,  3,  7, 21, 107,  95,  11,  195,  289, 2517,  973/)
    jkdata(1:13,1107) = (/1, 3, 7,  3, 29, 13,  65,  17, 409,  399, 1187,  733, 4821/)
    jkdata(1:13,1108) = (/1, 3, 5,  3, 17, 49, 101,  13, 275, 1003,  867, 1535, 2377/)
    jkdata(1:13,1109) = (/1, 3, 3,  1, 13, 61,  59, 243,  63,  121, 1535, 2175, 1673/)
    jkdata(1:13,1110) = (/1, 3, 3,  3,  3, 39,  35, 207, 441,  501,  575, 3613,    1/)
    jkdata(1:13,1111) = (/1, 1, 3, 15, 17, 15,  15, 187,  15,  155,  183, 3019, 6541/)

    v(1:w, 1) = 1
    v(1:min(13,w), 1:min(1111,s)) = jkdata(1:min(13,w), 1:min(1111,s))
    
   end subroutine init_m_values_joe_kuo


  ! Only for debugging...
  !
  subroutine print_direction_matrix(v)

    integer(kind=i4b), dimension(:,:), intent(in) :: v

    integer(kind=i4b) :: maxcol
    integer(kind=i4b) :: i

    maxcol = size(v, 1)

    do i=1,maxcol
      print *, v(i,:)
    end do 

  end subroutine print_direction_matrix


  ! Print out the type of Sobol sequence
  !
  subroutine print_sobol_type(my_soboltype)

    type(soboltype), intent(in) :: my_soboltype

    write(unit=*, fmt="(A70,  A)")                                           &
      "  Direction numbers:                                               ", &
      my_soboltype%dirnumbers
    write(unit=*, fmt="(A70,  A)")                                           &
      "  Primitive polynomials:                                           ", &
      my_soboltype%poly_order
    write(unit=*, fmt="(A70, L2)")                                           &
      "  Using all ones for initial direction numbers of first dimension? ", &
      my_soboltype%use_ones
    write(unit=*, fmt="(A70, L2)")                                           &
      "  Using Antonov-Saleev method?                                     ", &
      my_soboltype%use_antonov_saleev

  end subroutine print_sobol_type


  subroutine free_sobol()

    if (allocated(sob_startindex)) then
      deallocate(sob_startindex)
    end if
    if (allocated(sob_n)) then
      deallocate(sob_n)
    end if
    if (allocated(sob_lastq)) then
      deallocate(sob_lastq)
    end if
    if (allocated(v)) then
      deallocate(v)
    end if

  end subroutine free_sobol

end module mod_sobol

! Tests for correctness:
!   [OK] Compare with the plots on page 77 of Press and Teukolsky in
!        'Computers and Physics' on page 77.
