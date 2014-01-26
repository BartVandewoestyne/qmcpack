! Module that implements the Halton sequence and it's different scramblings.
!
! References:
!
!   [1] Halton, John H., `Algorithm 247: Radical-inverse quasi-random point
!       sequence' in `Communications of the ACM' vol. 7, nb. 12, 1964,
!       pages 701-702.
!
!   [2] Struckmeier, Jens, `Fast Generation of Low-Discrepancy Sequences',
!       in `Journal of Computational and Applied Mathematics', vol. 61, nb. 1,
!       1995, pages 29-41.
!
!   [3] Kolar, Miroslav and O'Shea, Seamus F., `Fast, portable, and reliable
!       algorithm for the calculation of Halton numbers.' in `Computers and
!       Mathematics with Applications', vol. 25, nb. 7, april 1993,
!       pages 3-13.

!   [4] Berblinger, Michael and Schlier, Christophe, `Monte Carlo integration
!       with quasi-random numbers: some experience', Computer Physics
!       Communications, vol. 66, nb 2-3, 1991, pages 157-166.
!
! TODO:
!   * add hal_step_values array so we can leap with diffent values in each
!     dimension.
!   * complete references.
!   * check out the implementation of ref [4]
!   * rename next_halton_scrambled next_halton_permutation_scrambled or so?

module mod_halton

  use numeric_kinds
  use mod_primes
  use mod_radical_inverse
  use mod_permutation
  use mod_number_theory
  use mod_utilities

  private

  public :: init_halton
  public :: random_startindex
  public :: halton
  public :: halton_bw
  public :: halton_chi_permuted
  public :: next_halton
  public :: next_halton_scrambled
  public :: next_halton_modified
  public :: next_halton_folded
  public :: next_halton_chi_permuted
  public :: next_halton_inversive_scrambled
  public :: next_halton_quadratic_scrambled
  public :: next_halton_phicf
  public :: free_halton
  public :: check_for_overflow
  public :: set_quadratic_scramble

  integer(kind=i4b), private                                    :: s
  integer(kind=i4b), dimension(:), allocatable, private         :: startindex
  integer(kind=i4b), dimension(:), allocatable, private         :: hal_primes
  integer(kind=i4b), dimension(:), allocatable, private         :: hal_n
  type(permutation_pointer), dimension(:), allocatable, private :: hal_ptable

  ! variables used for the quadratic scrambling method
  integer(kind=i4b), dimension(:), allocatable, private         :: qu, qv, qw

  ! The multipliers used by Chi and Mascagni
  integer(kind=i4b), dimension(40), parameter, private  :: w = (/ &
                                          1,   2,   3,   3,   8,  11,  12, &
                                         14,   7,  18,  12,  13,  17,  18, &
                                         29,  14,  18,  43,  41,  44,  40, &
                                         30,  47,  65,  71,  28,  40,  60, &
                                         79,  89,  56,  50,  52,  61, 108, &
                                         56,  66,  63,  60,  66 /)

  ! The k-modifiers from Atanassov's modified Halton sequence
  integer(kind=i4b), dimension(300), parameter, private :: hal_kmod1 = (/ &
       1,    1,    2,    4,    2,   11,   14,    7,   13,    6, &
       9,   25,   39,   20,   23,   36,   32,   18,   36,   40, &
      28,   69,   37,   86,   10,   67,   78,   37,   17,   81, &
      72,    6,   97,  105,   75,   53,   31,  149,   65,    1, &
     104,  105,  128,   30,  165,  137,   49,   31,   86,   12, &
     112,  176,   14,   12,   70,   83,  208,  174,   30,  156, &
     228,  271,  163,  263,   24,    1,  103,   41,  100,  289, &
     213,  318,    8,   62,  123,  166,  237,   13,  334,  186, &
     405,   85,  319,  112,  130,  376,   61,  161,  169,  308, &
     308,  426,    1,   21,   93,  221,  220,   98,  216,  336, &
     351,  367,  457,  516,  257,   76,  334,  190,  487,   21, &
     525,  595,  500,  198,  236,  402,  615,  306,  354,  312, &
      81,  658,  202,  220,  641,   49,   94,  588,  175,  434, &
     145,  538,  448,  369,  395,  187,  509,    6,   33,  454, &
     785,  718,  162,  464,  518,  225,  716,  287,  705,  680, &
     696,  734,  321,  588,  316,   26,  396,  572,  339,  303, &
     736,  189,  729,  730,  326,  895,  727,  201,  634,  673, &
     609,  615,  148,  805,  777,  117,  225,  805,   21,  674, &
     842,  414,  906,  121,  924,  458,  973,   81,  815, 1142, &
     541,  806, 1005, 1078, 1153,  329,  205, 1186, 1191,  643, &
     509,  468, 1177,  802,  988, 1235,  380, 1079, 1032,  247, &
     206, 1014,  725,  242, 1095,  370,  748, 1335,  469,  561, &
    1302,  722,  274,  465,  748,  349, 1358, 1148,  777,   52, &
    1068, 1397, 1210,  142, 1286,   95, 1012,  805, 1098, 1019, &
    1172, 1096,  263,  715,  349, 1406,  832,  415,  575,  826, &
    1118,  864,  899, 1138, 1089, 1380, 1353, 1373,  296, 1032, &
     576, 1009,  508,  925, 1223,   68, 1115,  455,    6,  505, &
     467,  336, 1004,  992,  234,  949, 1276, 1058, 1102,  991, &
    1454, 1303,  545,  180,  968, 1859, 1207, 1372,  561, 1463, &
    1632,  337,  971, 1415,  491,  142,  804,  839,  920, 1191 /)

  integer(kind=i4b), dimension(300), parameter, private :: hal_kmod2 = (/ &
     432,  101, 1643,  320,   42,  534,  521,  229, 1648,  561, &
    1569,  833,  707,  462,  526,  456,  547, 1473, 1818, 1220, &
     333,  131,  988,  992,  809, 1856,  773,  570, 2031,  252, &
    1315, 2158,   81, 1908,  588, 1939, 1051,  619, 1415,  873, &
    1761, 1923,  360,  934,  781,  725, 1175, 2135, 1914,  639, &
    1725,  443, 2127, 1259, 2173, 1368, 2091,  837,  138, 1047, &
    1394, 1277,  905, 1771, 2109, 2275, 1520, 1104, 1853,  250, &
    1637, 1866, 2005, 1088, 1767,  240, 1140,  825,  677, 2473, &
    1760, 1495, 2303, 2559, 1465, 1043, 1045, 1606, 1780, 2265, &
    1499, 1578, 2060, 1968, 2277,  739,  998,  409, 2323, 2100, &
     672,  565,  528, 2556,  132,  355,  429, 1847,   40,  441, &
    2014,   98, 2453, 1545,  962,  243, 1782,  323, 2578, 1273, &
    1887, 1399, 1613, 2045, 1871,  540, 1330, 1348, 1661, 1947, &
    2397,   50, 1985, 2275,  334, 1526, 2830,  621, 2313, 2213, &
     854, 2340, 2136, 1715,  788, 1194, 2801,  254,  100, 2479, &
    1923, 1637, 1074, 2687,  723,  127, 2468, 2904,  613, 2288, &
    1632, 2765, 1742,  279, 2878,  661, 1871,  553,  720,  726, &
     464,   38, 1412, 2545, 2396,  944, 2660, 2225, 1413,  351, &
    3064, 2323,  706,  778,  767, 1255, 1036, 2315, 1652,  194, &
     191, 2791, 1311, 2873, 1319, 1106, 1327, 1967, 2797, 1769, &
    1235, 2206, 2524, 3253,  926, 2976,  789, 3469, 1199, 2163, &
     138,  141, 1624,  813,  274,  137, 3495, 2698, 1607, 2079, &
    2720, 2559,  233, 1715,  509,  495, 3434, 2909, 3180,  246, &
     375, 3662, 1548, 1995, 2872,  258, 1871, 1326,  736,  450, &
    2440, 3228, 2899,  871, 3183, 2465,  387, 3258, 1310,  179, &
    3019, 1361, 2846, 3301, 1390, 2694,   70,  704, 3783, 1971, &
    3121, 1549, 3247, 1301, 2102, 2413,  330, 3762, 2212,  420, &
    3774,   69, 1077, 3105, 2207,  328,  189, 1693,  369, 2152, &
    3183, 1778, 4038, 3071,  154, 1817, 2475,  446,  213,  276, &
     390, 1725, 3516,  192, 2282, 1997, 2353,    5,  639,  604 /)

  integer(kind=i4b), dimension(300), parameter, private :: hal_kmod3 = (/ &
    3082, 1224, 4415, 3195, 3471,  575, 3521, 4091, 3164, 2794, &
    1738, 1634, 2572, 2276, 2286,  206, 2399, 3199, 1752,  648, &
    2957, 3251, 1984, 2066, 1862,  463,  245, 4206, 4116,  109, &
    1824, 2999, 1866, 2037, 4185, 1623,  905,  298, 3016, 2112, &
    1632, 4193, 2627, 2325,  922,  868,  605, 3484, 3804, 2294, &
    1684,  905, 1800, 4385, 1558, 1291, 2339, 1962, 4809, 2170, &
    3091,  842, 3788, 4181, 3899,  766, 2868, 2335,  666, 4147, &
    2564,   55, 3026, 2361,  728,  253, 4612, 4437,  385, 4832, &
    3385, 2719, 4906, 3165, 3651,  744, 4576, 4940, 3764, 1705, &
     649, 1747, 1599, 3308, 1303, 4976, 2346, 1229, 2236, 2734, &
    1695, 1147, 4889, 2226, 4251, 1364,   96,  183, 3048, 1885, &
     623, 3753, 5150, 4404, 5248, 1182,  716, 4943, 1674, 4093, &
    2649, 1558, 3613, 3937, 3234, 3621, 4582, 3569, 3022, 5168, &
    4439,  721, 4776, 4783, 2761, 2466,  387, 2949, 2164, 5610, &
     646, 3145, 1376, 3606, 5390, 3393, 4470, 5483, 3414, 3585, &
    2626, 4477, 3916, 2974, 4316, 1513, 1361, 5109, 2831, 5507, &
    2362,  922,   18, 1203, 4732, 2483,  661, 3110, 1819,  108, &
    2654,  570, 1493, 4440,  223, 3899, 1685,  263, 2544, 5165, &
    4715, 5865, 4420, 2038, 2464, 3446,  644, 5092, 1753, 1882, &
    5444, 3666, 4819, 5586, 2327, 3417, 1395, 1601, 2069, 2597, &
    3044,   53, 1120, 4118, 5045, 4381, 4917, 3010, 6181, 1734, &
    2504, 2907, 2520, 6242, 4896, 1165, 2364, 2317, 1779, 5191, &
     126, 3439, 5013, 6237, 2351, 3635,   98, 4073, 2875,  708, &
     279, 5358, 3309, 2325, 2742, 5646, 5027, 3375,  388, 3503, &
    6187,  520,  472, 6492,  849, 1818, 1545, 5142, 6340, 4465, &
     144, 4745, 6081,  247, 3232, 3243,  875,  656, 4053, 6258, &
    4459, 2832, 5417, 2700, 5350, 2490, 2664, 2868, 6468, 4494, &
    5839, 1066, 1692, 3100, 5851, 1613, 5230, 5428, 5958, 4623, &
    2690, 1325, 3381, 2510,  960, 1591, 4203, 2056, 5893, 2446, &
    2333, 4171, 1149, 3882, 5054, 2739, 1224,  132, 3164, 6832 /)

  integer(kind=i4b), dimension(100), parameter, private :: hal_kmod4 = (/ &
    2615, 2678,  100, 4011, 5621, 2376, 2466,  928, 2584, 6645, &
    4972, 5967, 5153, 4730, 5272, 6313, 1624, 5368, 7004,  178, &
    2666, 3469, 3272, 4281, 2578, 6912, 5537, 3621, 1716, 2467, &
    2578, 5819, 4750,  525, 5472, 5469, 6937, 6039, 4359, 7243, &
    4059, 6766, 4263,  449, 2837, 7118, 7354, 2104, 5666, 1164, &
    5593, 7230, 4247, 3853, 6740, 5248, 4392, 4465, 4460, 6547, &
    1555,  490, 5545, 5158, 4326, 2824, 2866, 1315,   77, 7538, &
    3233, 4728, 1093, 4754, 7233, 1956, 2653, 6865, 3309, 6193, &
    1090, 6912, 6606, 2395, 2506, 5859, 4204, 1506,  643,  557, &
    5590, 1968, 5130,  414, 3180,   32, 1687, 6198, 3966,  313 /)

  integer(kind=i4b), dimension(1000), parameter, public :: &
        hal_kmod = (/ hal_kmod1, hal_kmod2, hal_kmod3, hal_kmod4 /)
        
  contains


    ! Default values:
    !   startindex: 1
    !   myprimes:   first s primes
    !   perm_type:  none
    subroutine init_halton(init_s, init_startindex, myprimes, perm_type)

      integer(kind=i4b), intent(in)                         :: init_s
      integer(kind=i4b), dimension(:), intent(in), optional :: init_startindex
      integer(kind=i4b), dimension(:), intent(in), optional :: myprimes
      character(len=*), intent(in), optional                :: perm_type

      ! Initialize the dimension of the points to be generated
      s = init_s

      ! Initialize the startindices
      if (allocated(startindex)) then
        deallocate(startindex)
      end if
      allocate(startindex(s))
      if (present(init_startindex)) then
        startindex = init_startindex
      else
        ! Use (1,...,1) as default startindex
        startindex = 1
      end if

      if (allocated(hal_n)) then
        deallocate(hal_n)
      end if
      allocate(hal_n(s))
      hal_n = startindex

      ! Initialize the primes to be used
      if (allocated(hal_primes)) then
        deallocate(hal_primes)
      end if
      allocate(hal_primes(s))
      if (present(myprimes)) then
        hal_primes = myprimes
      else
        ! Use the first s primes by default
        hal_primes = primes(s)
      end if

      if (present(perm_type)) then

        if (allocated(hal_ptable)) then
          !call free_permutation_table(hal_ptable)
          deallocate(hal_ptable)
        end if

        select case(perm_type)

          case ("Braaten_Weller")

            if (s>16) then
              print *, "ERROR: Braaten and Weller's scrambled Halton sequence"
              print *, "       is only defined up to 16 dimensions!"
            else
              allocate(hal_ptable(s))
              hal_ptable = init_perm_braaten_weller(s)
            end if

          case ("BW_extended")

            if (s>64) then
              print *, "ERROR: Extended Braaten and Weller scrambling"
              print *, "       is only available up to 64 dimensions!"
            else
              allocate(hal_ptable(s))
              hal_ptable = init_perm_bw_extended(s)
            end if

          case ("Faure")

            allocate(hal_ptable(s))
            hal_ptable = init_perm_faure(s, .true.)

          case ("MCL")

            if (s>16) then
              print *, "ERROR: Tuffin's MCL scrambled Halton sequence"
              print *, "       is only defined up to 16 dimensions!"
            else
              allocate(hal_ptable(s))
              hal_ptable = init_perm_mcl(s)
            end if

          case ("MCL_star")

            if (s>16) then
              print *, "ERROR: Tuffin's MCL* scrambled Halton sequence"
              print *, "       is only defined up to 16 dimensions!"
            else
              allocate(hal_ptable(s))
              hal_ptable = init_perm_mcl_star(s)
            end if

          case ("reverse")
            
            allocate(hal_ptable(s))
            hal_ptable = init_perm_reverse(s)

          case default
            print *, "ERROR: Unknown permutation table type for the Halton!!!"

        end select
      end if

    end subroutine init_halton


    ! The standard Halton sequence
    ! The resulting x contains its dimensions in its columns.
    !
    subroutine halton(n, s, start, myprimes, x)
      integer(kind=i4b), intent(in)                         :: n
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), dimension(:), intent(in), optional :: start
      integer(kind=i4b), dimension(:), intent(in), optional :: myprimes
      real(kind=qp), dimension(:, :), intent(out)           :: x

      integer(kind=i4b) :: i, j

      call init_halton(s, start, myprimes)

      do i=1,s
        do j=0,n-1
          x(j+1,i) = radical_inverse(startindex(i) + j, hal_primes(i))
        end do
      end do

      call free_halton()

    end subroutine halton

    ! Next standard Halton point
    subroutine next_halton(x)
      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b) :: i

      call check_for_overflow()

      do i=1,s
        x(i) = radical_inverse(hal_n(i), hal_primes(i))
      end do

      hal_n = hal_n + 1

    end subroutine next_halton


    ! The Halton sequence as proposed by Braaten and Weller
    subroutine halton_bw(n, s, start, x)
      integer(kind=i4b), intent(in)                         :: n
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), dimension(:), intent(in), optional :: start
      real(kind=qp), dimension(:,:), intent(out)            :: x

      integer(kind=i4b) :: i, j

      if (s>16) then

        print *, "ERROR: Braaten and Weller's scrambled Halton sequence is"
        print *, "       only defined up to 16 dimensions."
        print *, "       No points were generated!!!"

      else

        ! Initialize the permutation table of Braaten and Weller
        allocate(hal_ptable(s))
        hal_ptable = init_perm_braaten_weller(s)

        ! Braaten and Weller use the first s primes
        call init_halton(s, start, primes(s))

        do i=1,s
          do j=0,n-1
            x(j+1,i) = radical_inverse_scrambled(startindex(i)+j, &
                                                 hal_primes(i),       &
                                                 hal_ptable(i)%p)
          end do
        end do

        call free_halton()

      end if

    end subroutine halton_bw


    ! Generate the next scrambled Halton point, with the digit-scrambling as
    ! defined in the permutation table hal_ptable
    subroutine next_halton_scrambled(x)
      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b) :: i


      call check_for_overflow()

      do i=1,s
        x(i) = radical_inverse_scrambled(hal_n(i),        &
                                         hal_primes(i),   &
                                         hal_ptable(i)%p)
      end do

      hal_n = hal_n + 1

    end subroutine next_halton_scrambled


    ! Default values:
    !   increased = .true.
    !   scrambled = .false. 
    subroutine next_halton_modified(x, increased, scrambled)
      real(kind=qp), dimension(:), intent(out) :: x
      logical, intent(in), optional            :: increased
      logical, intent(in), optional            :: scrambled

      integer(kind=i4b) :: i
    
      if (s > 1000) then
        print *, "ERROR: Atanassov's modified Halton sequence is only"
        print *, "       defined up to 1000 dimensions."
        print *, "       No points were generated!!"
      else

        call check_for_overflow()

        do i=1,s

          call radical_inverse_modified(hal_n(i),                  &
                                        hal_primes(i),             &
                                        hal_kmod(i),               &
                                        increased, scrambled, x(i))
                                        
        end do

        hal_n = hal_n + 1

      end if

    end subroutine next_halton_modified


    subroutine next_halton_folded(x)
      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b) :: i

      call check_for_overflow()

      do i=1,s
        x(i) = radical_inverse_folded(hal_n(i),     &
                                      hal_primes(i))
      end do
      
      hal_n = hal_n + 1

    end subroutine next_halton_folded


    ! References:
    !   [1] 'On the scrambled Halton sequence', Mascagni, Michael and Chi,
    !       Hongmei, Monte Carlo Methods and Applications, Vol. 10, No. 3-4,
    !       pp. 435--442 (2004)
    subroutine halton_chi_permuted(n, s, start, x)
      integer(kind=i4b), intent(in)                         :: n
      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), dimension(:), intent(in), optional :: start
      real(kind=qp), dimension(:,:), intent(out)            :: x
   
      integer(kind=i4b) :: i, j

      ! Mascagni and Chi use the first s primes
      call init_halton(s, start, primes(s))

      do i=1,s
        do j=0,n-1
          x(j+1,i) = radical_inverse_lin_scrambled(startindex(i)+j, &
                                                   hal_primes(i),       &
                                                   w(i), 0)
        end do
      end do

      call free_halton()
    
    end subroutine halton_chi_permuted


    subroutine next_halton_chi_permuted(x)
      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b) :: i

      call check_for_overflow()

      do i=1,s
        x(i) = radical_inverse_lin_scrambled(hal_n(i),      &
                                             hal_primes(i), &
                                             w(i), 0)
      end do
      hal_n = hal_n + 1

    end subroutine next_halton_chi_permuted


    subroutine next_halton_inversive_scrambled(x)
      real(kind=qp), dimension(:), intent(out) :: x

      call check_for_overflow()

      ! TODO: How to choose the w's???

      x = radical_inverse_inversive(hal_n, hal_primes, w(1:size(x)), 0)
       
      hal_n = hal_n + 1

    end subroutine next_halton_inversive_scrambled


    subroutine next_halton_quadratic_scrambled(x)
      real(kind=qp), dimension(:), intent(out) :: x

      call check_for_overflow()

      x = radical_inverse_quadratic(hal_n, hal_primes, qu, qv, qw)
       
      hal_n = hal_n + 1

    end subroutine next_halton_quadratic_scrambled


    subroutine next_halton_phicf(x)
      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b) :: i

      call check_for_overflow()

      do i=1,s
        x(i) = radical_inverse_phicf(hal_n(i), hal_primes(i))
      end do
      hal_n = hal_n + 1

    end subroutine next_halton_phicf


    ! Fill startindex with random start integers
    !
    subroutine random_startindex(startindex, imax)
      integer(kind=i4b), dimension(:), intent(out) :: startindex
      integer(kind=i4b), intent(in), optional      :: imax

      real(kind=qp), dimension(size(startindex)) :: x

      call random_number(x)

      if (present(imax)) then
        startindex = int(imax*x, kind=i4b)
      else
        startindex = int((huge(i4b)-1)*x, kind=i4b)
      end if

    end subroutine random_startindex


    subroutine free_halton()

      if (allocated(startindex)) then
        deallocate(startindex)
      end if
      if (allocated(hal_n)) then
        deallocate(hal_n)
      end if
      if (allocated(hal_primes)) then
        deallocate(hal_primes)
      end if
      if (allocated(hal_ptable)) then

        ! Deallocate the permutations
        call free_permutation_table(hal_ptable)

        ! Deallocate the table containing the pointers itself
        deallocate(hal_ptable)

      end if
      if (allocated(qv)) then
        deallocate(qu, qv, qw)
      end if

    end subroutine free_halton


    subroutine set_quadratic_scramble(u, v, w)
      integer(kind=i4b), dimension(:), intent(in) :: u, v, w

      if (allocated(qu)) then
        deallocate(qu, qv, qw)
      end if
      allocate(qu(size(u)), qv(size(v)), qw(size(w)))
      qu = u
      qv = v
      qw = w

    end subroutine set_quadratic_scramble


    ! Check if we are going to generate more points than allowed by
    ! our architecture
    subroutine check_for_overflow()

      if (huge(i4b)-maxval(hal_n) == 0) then
        write(unit=*, fmt="(A)") "ERROR: maximum number of Halton points reached!"
      end if

    end subroutine check_for_overflow


end module mod_halton
