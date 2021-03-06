! Module to implement methods concerning primitive polynomials over GF(2).
!
! References:
!   [1] TODO
!
! TODO:
!   * Add the list of primitive polynomials from http://fchabaud.free.fr/

module mod_primpoly

use numeric_kinds

private

public :: get_degree
public :: get_primpolys

private :: get_primpolys_bratley_fox
private :: get_primpolys_jaeckel
private :: get_primpolys_joe_kuo
private :: get_primpolys_numrecipes
private :: get_primpolys_press_teukolsky
private :: get_primpolys_sobol76
private :: get_primpolys_watson


contains

  ! Returns the degree of the polynomial P which is specified in binary
  ! encoding.  The binary encoding of polynomials P works as follows:
  !
  ! 1.X^3 + 0.X^2 + 1.X^1 + 1 ====> 1.2^3 + 0.2^2 + 1.2^1 + 1 = 11
  !
  function get_degree(p) result (d)

    integer(kind=i4b), intent(in) :: p
    integer(kind=i4b)             :: d 

    integer(kind=i4b) :: temp

    temp = p/2  ! we rely on integer division here!
    d = 0
    do
      if (temp == 0) then
        exit
      end if
      d = d+1
      temp = temp/2  ! We rely on integer division here!
    end do
    
  end function get_degree


  subroutine get_primpolys(s, polys, order)
    integer(kind=i4b), intent(in)                :: s
    integer(kind=i4b), dimension(:), intent(out) :: polys
    character(len=*), intent(in)                 :: order

    select case (order)
    case ("Matlab")
      stop "Program aborted: Matlab primitive polynomials not implemented yet"
    case ("NumericalRecipes")
      call get_primpolys_numrecipes(s, polys)
    case ("PressTeukolsky")
      call get_primpolys_press_teukolsky(s, polys)
    case ("JoeKuo")
      call get_primpolys_joe_kuo(s, polys)
    case ("BratleyFox")
      call get_primpolys_bratley_fox(s, polys)
    case ("Winiarski")
      stop "Program aborted: Winiarski primitive polynomials not implemented yet"
    case ("Jaeckel")
      call get_primpolys_jaeckel(s, polys)
    case ("Sobol76")
      call get_primpolys_sobol76(s, polys)
    case ("Watson")
      call get_primpolys_watson(s, polys)
    case default
      stop "Unknown order for the primitive polynomials."
    end select

  end subroutine get_primpolys


  ! The primitive polynomials used by the
  ! Numerical Recipes implementation
  subroutine get_primpolys_numrecipes(s, polys)
    integer(kind=i4b), intent(in)                :: s
    integer(kind=i4b), dimension(:), intent(out) :: polys

    integer(kind=i4b), dimension(6) :: primdata

    if (s>6) then
      stop "Numerical Recipes only specifies the first 6 primitive polynomials!"
    end if
    primdata = (/ 3, 7, 11, 13, 19, 25 /)
    polys = primdata(1:s)

  end subroutine get_primpolys_numrecipes


  ! Return the polynomials in the order as Press and Teukolsky specify them
  subroutine get_primpolys_press_teukolsky(s, polys)
    integer(kind=i4b), intent(in)                :: s
    integer(kind=i4b), dimension(:), intent(out) :: polys

    integer(kind=i4b), dimension(160) :: primdata

    if (s>160) then
      stop "Press and Teukolsky only specify the first 160 primitive polynomials!"
    end if

    ! These are the representations *with* the higher order and the unit
    ! bit!!!
    primdata(1) = 3
    primdata(2) = 7
    primdata(3:4) = (/ 11, 13 /) 
    primdata(5:6) = (/ 19, 25 /)
    primdata(7:12) = (/ 37, 41, 47, 55, 59, 61 /)
    primdata(13:18) = (/ 67, 91, 97, 103, 109, 115 /)
    primdata(19:36) = (/ 131, 137, 143, 145, 157, 167, 171, 185, 191, 193, &
                         203, 211, 213, 229, 239, 241, 247, 253 /)
    primdata(37:100) = (/ 285, 299,  301,  333, 351, 355, 357, 361, 369, 391, &
                          397, 425,  451,  463, 487, 501, 529, 539, 545, 557, &
                          563, 601,  607,  617, 623, 631, 637, 647, 661, 675, &
                          677, 687,  695,  701, 719, 721, 731, 757, 761, 787, &
                          789, 799,  803,  817, 827, 847, 859, 865, 875, 877, &
                          883, 895,  901,  911, 949, 953, 967, 971, 973, 981, &
                          985, 995, 1001, 1019 /)
    primdata(101:160) = (/ 1033, 1051, 1063, 1069, 1125, 1135, 1153, 1163, &
                           1221, 1239, 1255, 1267, 1279, 1293, 1305, 1315, &
                           1329, 1341, 1347, 1367, 1387, 1413, 1423, 1431, &
                           1441, 1479, 1509, 1527, 1531, 1555, 1557, 1573, &
                           1591, 1603, 1615, 1627, 1657, 1663, 1673, 1717, &
                           1729, 1747, 1759, 1789, 1815, 1821, 1825, 1849, &
                           1863, 1869, 1877, 1881, 1891, 1917, 1933, 1939, &
                           1969, 2011, 2035, 2041 /)

    polys = primdata(1:s)

  end subroutine get_primpolys_press_teukolsky 


  ! Return the primitive polynomials in the order as Bratley and Fox specify them.
  !
  ! Note:
  !   The order of Bratley and Fox should also be the order as from
  !   Sobol in 'USSR COMPUT. MATHS. MATH. PHYS. 16 (1977), pages 236-242'.
  !
  subroutine get_primpolys_bratley_fox(s, polys)
    integer(kind=i4b), intent(in)                :: s
    integer(kind=i4b), dimension(:), intent(out) :: polys

    integer(kind=i4b), dimension(39) :: primdata

    if (s>39) then
      stop "Bratley and Fox only specify the first 39 primitive polynomials!"
    end if
    primdata(1:39) = (/ 3,7,11,13,19,25,37,59,47,                 &
                       61,55,41,67,97,91,109,103,115,131,         &
                       193,137,145,143,241,157,185,167,229,171,   &
                       213,191,253,203,211,239,247,285,369,299 /)
    polys = primdata(1:s)
    
  end subroutine get_primpolys_bratley_fox


  ! Return the primitive polynomials in the order as Joe and Kuo specify them.
  !
  ! Note: 
  !   Joe and Kuo don't specify a primitive polynomial for the first dimension.
  !   All values from 1 to MAXCOL in the direction matrix are set to 1!  The
  !   polynomials as given here, are the polynomials that are used for
  !   constructing the direction matrix, starting from dimension 2!!!
  !
  subroutine get_primpolys_joe_kuo(s, polys)
    integer(kind=i4b), intent(in)                :: s
    integer(kind=i4b), dimension(:), intent(out) :: polys

    integer(kind=i4b), dimension(1110) :: primdata

   if (s>1110) then
     stop "Joe and Kuo only specify the first 1110 primitive polynomials!"
   end if

   primdata(1:210) = (/ 3,7,11,13,19,25,37,59,47,61,55,41,67,97,91, &
          109,103,115,131,193,137,145,143,241,157,185,167,229,171,213, &
          191,253,203,211,239,247,285,369,299,301,333,351,355,357,361, &
          391,397,425,451,463,487,501,529,539,545,557,563,601,607,617, &
          623,631,637,647,661,675,677,687,695,701,719,721,731,757,761, &
          787,789,799,803,817,827,847,859,865,875,877,883,895,901,911, &
          949,953,967,971,973,981,985,995,1001,1019,1033,1051,1063, &
          1069,1125,1135,1153,1163,1221,1239,1255,1267,1279,1293,1305, &
          1315,1329,1341,1347,1367,1387,1413,1423,1431,1441,1479,1509, &
          1527,1531,1555,1557,1573,1591,1603,1615,1627,1657,1663,1673, &
          1717,1729,1747,1759,1789,1815,1821,1825,1849,1863,1869,1877, &
          1881,1891,1917,1933,1939,1969,2011,2035,2041,2053,2071,2091, &
          2093,2119,2147,2149,2161,2171,2189,2197,2207,2217,2225,2255, &
          2257,2273,2279,2283,2293,2317,2323,2341,2345,2363,2365,2373, &
          2377,2385,2395,2419,2421,2431,2435,2447,2475,2477,2489,2503, &
          2521,2533,2551,2561,2567,2579,2581,2601,2633,2657,2669 /)
   primdata(211:400) = (/ 2681,2687,2693,2705,2717,2727,2731,2739, &
          2741,2773,2783,2793,2799,2801,2811,2819,2825,2833,2867,2879, &
          2881,2891,2905,2911,2917,2927,2941,2951,2955,2963,2965,2991, &
          2999,3005,3017,3035,3037,3047,3053,3083,3085,3097,3103,3159, &
          3169,3179,3187,3205,3209,3223,3227,3229,3251,3263,3271,3277, &
          3283,3285,3299,3305,3319,3331,3343,3357,3367,3373,3393,3399, &
          3413,3417,3427,3439,3441,3475,3487,3497,3515,3517,3529,3543, &
          3547,3553,3559,3573,3589,3613,3617,3623,3627,3635,3641,3655, &
          3659,3669,3679,3697,3707,3709,3713,3731,3743,3747,3771,3791, &
          3805,3827,3833,3851,3865,3889,3895,3933,3947,3949,3957,3971, &
          3985,3991,3995,4007,4013,4021,4045,4051,4069,4073,4179,4201, &
          4219,4221,4249,4305,4331,4359,4383,4387,4411,4431,4439,4449, &
          4459,4485,4531,4569,4575,4621,4663,4669,4711,4723,4735,4793, &
          4801,4811,4879,4893,4897,4921,4927,4941,4977,5017,5027,5033, &
          5127,5169,5175,5199,5213,5223,5237,5287,5293,5331,5391,5405, &
          5453,5523,5573,5591,5597,5611,5641,5703,5717,5721,5797,5821, &
          5909,5913 /)
   primdata(401:590) = (/ 5955,5957,6005,6025,6061,6067,6079,6081, &
          6231,6237,6289,6295,6329,6383,6427,6453,6465,6501,6523,6539, &
          6577,6589,6601,6607,6631,6683,6699,6707,6761,6795,6865,6881, &
          6901,6923,6931,6943,6999,7057,7079,7103,7105,7123,7173,7185, &
          7191,7207,7245,7303,7327,7333,7355,7365,7369,7375,7411,7431, &
          7459,7491,7505,7515,7541,7557,7561,7701,7705,7727,7749,7761, &
          7783,7795,7823,7907,7953,7963,7975,8049,8089,8123,8125,8137, &
          8219,8231,8245,8275,8293,8303,8331,8333,8351,8357,8367,8379, &
          8381,8387,8393,8417,8435,8461,8469,8489,8495,8507,8515,8551, &
          8555,8569,8585,8599,8605,8639,8641,8647,8653,8671,8675,8689, &
          8699,8729,8741,8759,8765,8771,8795,8797,8825,8831,8841,8855, &
          8859,8883,8895,8909,8943,8951,8955,8965,8999,9003,9031,9045, &
          9049,9071,9073,9085,9095,9101,9109,9123,9129,9137,9143,9147, &
          9185,9197,9209,9227,9235,9247,9253,9257,9277,9297,9303,9313, &
          9325,9343,9347,9371,9373,9397,9407,9409,9415,9419,9443,9481, &
          9495,9501,9505,9517,9529,9555,9557,9571,9585,9591,9607,9611, &
          9621,9625 /)
   primdata(591:764) = (/ 9631,9647,9661,9669,9679,9687,9707,9731, &
          9733,9745,9773,9791,9803,9811,9817,9833,9847,9851,9863,9875, &
          9881,9905,9911,9917,9923,9963,9973,10003,10025,10043,10063, &
          10071,10077,10091,10099,10105,10115,10129,10145,10169,10183, &
          10187,10207,10223,10225,10247,10265,10271,10275,10289,10299, &
          10301,10309,10343,10357,10373,10411,10413,10431,10445,10453, &
          10463,10467,10473,10491,10505,10511,10513,10523,10539,10549, &
          10559,10561,10571,10581,10615,10621,10625,10643,10655,10671, &
          10679,10685,10691,10711,10739,10741,10755,10767,10781,10785, &
          10803,10805,10829,10857,10863,10865,10875,10877,10917,10921, &
          10929,10949,10967,10971,10987,10995,11009,11029,11043,11045, &
          11055,11063,11075,11081,11117,11135,11141,11159,11163,11181, &
          11187,11225,11237,11261,11279,11297,11307,11309,11327,11329, &
          11341,11377,11403,11405,11413,11427,11439,11453,11461,11473, &
          11479,11489,11495,11499,11533,11545,11561,11567,11575,11579, &
          11589,11611,11623,11637,11657,11663,11687,11691,11701,11747, &
          11761,11773,11783,11795,11797,11817,11849,11855,11867,11869, &
          11873,11883,11919 /)
   primdata(765:935) = (/ 11921,11927,11933,11947,11955,11961, &
          11999,12027,12029,12037,12041,12049,12055,12095,12097,12107, &
          12109,12121,12127,12133,12137,12181,12197,12207,12209,12239, &
          12253,12263,12269,12277,12287,12295,12309,12313,12335,12361, &
          12367,12391,12409,12415,12433,12449,12469,12479,12481,12499, &
          12505,12517,12527,12549,12559,12597,12615,12621,12639,12643, &
          12657,12667,12707,12713,12727,12741,12745,12763,12769,12779, &
          12781,12787,12799,12809,12815,12829,12839,12857,12875,12883, &
          12889,12901,12929,12947,12953,12959,12969,12983,12987,12995, &
          13015,13019,13031,13063,13077,13103,13137,13149,13173,13207, &
          13211,13227,13241,13249,13255,13269,13283,13285,13303,13307, &
          13321,13339,13351,13377,13389,13407,13417,13431,13435,13447, &
          13459,13465,13477,13501,13513,13531,13543,13561,13581,13599, &
          13605,13617,13623,13637,13647,13661,13677,13683,13695,13725, &
          13729,13753,13773,13781,13785,13795,13801,13807,13825,13835, &
          13855,13861,13871,13883,13897,13905,13915,13939,13941,13969, &
          13979,13981,13997,14027,14035,14037,14051,14063,14085,14095, &
          14107,14113,14125,14137,14145 /)
   primdata(936:1106) = (/ 14151,14163,14193,14199,14219,14229, &
          14233,14243,14277,14287,14289,14295,14301,14305,14323,14339, &
          14341,14359,14365,14375,14387,14411,14425,14441,14449,14499, &
          14513,14523,14537,14543,14561,14579,14585,14593,14599,14603, &
          14611,14641,14671,14695,14701,14723,14725,14743,14753,14759, &
          14765,14795,14797,14803,14831,14839,14845,14855,14889,14895, &
          14909,14929,14941,14945,14951,14963,14965,14985,15033,15039, &
          15053,15059,15061,15071,15077,15081,15099,15121,15147,15149, &
          15157,15167,15187,15193,15203,15205,15215,15217,15223,15243, &
          15257,15269,15273,15287,15291,15313,15335,15347,15359,15373, &
          15379,15381,15391,15395,15397,15419,15439,15453,15469,15491, &
          15503,15517,15527,15531,15545,15559,15593,15611,15613,15619, &
          15639,15643,15649,15661,15667,15669,15681,15693,15717,15721, &
          15741,15745,15765,15793,15799,15811,15825,15835,15847,15851, &
          15865,15877,15881,15887,15899,15915,15935,15937,15955,15973, &
          15977,16011,16035,16061,16069,16087,16093,16097,16121,16141, &
          16153,16159,16165,16183,16189,16195,16197,16201,16209,16215, &
          16225,16259,16265,16273,16299 /)
   primdata(1107:1110) = (/ 16309,16355,16375,16381 /)

   ! Only return the part we need
   polys = primdata(1:s)
    
  end subroutine get_primpolys_joe_kuo


  ! References:
  !  [1] http://www.nr.com/contrib/jaeckel/ppm2.txt and
  !  [2] http://www.nr.com/contrib/
  !  [3] The CD accompanying Peter Jaeckel's book.
  !  [4] 'Monte Carlo methods in finance', Peter Jaeckel
  !
  subroutine get_primpolys_jaeckel(s, polys)
    integer(kind=i4b), intent(in)                :: s
    integer(kind=i4b), dimension(:), intent(out) :: polys

    integer(kind=i4b), dimension(31) :: primdata

    primdata(1:1) = (/ 3 /)
    primdata(2:2) = (/ 7 /)
    primdata(3:4) = (/ 11, 13 /)
    primdata(5:6) = (/ 19, 25 /)
    primdata(7:12) = (/ 37, 41, 47, 55, 59, 61 /)
    primdata(13:18) = (/ 67, 91, 97, 103, 109, 115 /)
    primdata(19:31) = (/ 131, 137, 143, 145, 157, 167, 171, 185, 191, 193, &
                         203, 211, 213 /)
    ! Only return the part we need
    polys = primdata(1:s)

  end subroutine get_primpolys_jaeckel


  ! The primitive polynomials from Sobol's 1976 paper
  !
  subroutine get_primpolys_sobol76(s, polys)
    integer(kind=i4b), intent(in)                :: s
    integer(kind=i4b), dimension(:), intent(out) :: polys

    integer(kind=i4b), dimension(15) :: primdata

   if (s>15) then
     stop "Sobol only specifies the first 15 primitive polynomials in his paper!"
   end if

    primdata(1:15) = (/ 3, 7, 11, 13, 19, 25, 37, 59, 47, 61, 55, 41, 67, 97, &
                        91 /)
    polys = primdata(1:s)

  end subroutine get_primpolys_sobol76


  ! Return the primitive polynomials in the order as specified by Watson.
  !
  ! Note that Watson only gives one polynomial for each degree, so he's not
  ! listing ALL primitive polynomials of a specific degree!!!
  ! ATTENTION: watch out, as these encodings can get very large...
  ! The last polynomials might not be represented well on your architecture... 
  subroutine get_primpolys_watson(s, polys)
    integer(kind=i4b), intent(in)                :: s
    integer(kind=i4b), dimension(:), intent(out) :: polys

    integer(kind=i4b), dimension(30) :: primdata

   if (s>30) then
     stop "You can only use the first 30 primitive polynomials specified by Watson!"
   end if

   primdata(1) = sum(2**(/1, 0/))
   primdata(2) = sum(2**(/2, 1, 0/))
   primdata(3) = sum(2**(/3, 1, 0/))
   primdata(4) = sum(2**(/4, 1, 0/))
   primdata(5) = sum(2**(/5, 2, 0/))
   primdata(6) = sum(2**(/6, 1, 0/))
   primdata(7) = sum(2**(/7, 1, 0/))
   primdata(8) = sum(2**(/8, 4, 3, 2, 0/))
   primdata(9) = sum(2**(/9, 4, 0/))
   primdata(10) = sum(2**(/10, 3, 0/))
   primdata(11) = sum(2**(/11, 2, 0/))
   primdata(12) = sum(2**(/12, 6, 4, 1, 0/))
   primdata(13) = sum(2**(/13, 4, 3, 1, 0/))
   primdata(14) = sum(2**(/14, 5, 3, 1, 0/))
   primdata(15) = sum(2**(/15, 1, 0/))
   primdata(16) = sum(2**(/16, 5, 3, 2, 0/))
   primdata(17) = sum(2**(/17, 3, 0/))
   primdata(18) = sum(2**(/18, 5, 2, 1, 0/))
   primdata(19) = sum(2**(/19, 5, 2, 1, 0/))
   primdata(20) = sum(2**(/20, 3, 0/))
   primdata(21) = sum(2**(/21, 2, 0/))
   primdata(22) = sum(2**(/22, 1, 0/))
   primdata(23) = sum(2**(/23, 5, 0/))
   primdata(24) = sum(2**(/24, 4, 3, 1, 0/))
   primdata(25) = sum(2**(/25, 3, 0/))
   primdata(26) = sum(2**(/26, 6, 2, 1, 0/))
   primdata(27) = sum(2**(/27, 5, 2, 1, 0/))
   primdata(28) = sum(2**(/28, 3, 0/))
   primdata(29) = sum(2**(/29, 2, 0/))
   primdata(30) = sum(2**(/30, 6, 4, 1, 0/))
   ! The rest overflows for the integer kind we use...
   !primdata(31) = sum(2**(/31, 3, 0/))
   !primdata(32) = sum(2**(/32, 7, 5, 3, 2, 1, 0/))
   !primdata(33) = sum(2**(/33, 6, 4, 1, 0/))
   !primdata(34) = sum(2**(/34, 7, 6, 5, 2, 1, 0/))
   !primdata(35) = sum(2**(/35, 2, 0/))
   !primdata(36) = sum(2**(/36, 6, 5, 4, 2, 1, 0/))
   !primdata(37) = sum(2**(/37, 5, 4, 3, 2, 1, 0/))
   !primdata(38) = sum(2**(/38, 6, 5, 1, 0/))
   !primdata(39) = sum(2**(/39, 4, 0/))
   !primdata(40) = sum(2**(/40, 5, 4, 3, 0/))
   !primdata(41) = sum(2**(/41, 3, 0/))
   !primdata(42) = sum(2**(/42, 5, 4, 3, 2, 1, 0/))
   !primdata(43) = sum(2**(/43, 6, 4, 3, 0/))
   !primdata(44) = sum(2**(/44, 6, 5, 2, 0/))
   !primdata(45) = sum(2**(/45, 4, 3, 1, 0/))
   !primdata(46) = sum(2**(/46, 8, 5, 3, 2, 1, 0/))
   !primdata(47) = sum(2**(/47, 5, 0/))
   !primdata(48) = sum(2**(/48, 7, 5, 4, 2, 1, 0/))
   !primdata(49) = sum(2**(/49, 6, 5, 4, 0/))
   !primdata(50) = sum(2**(/50, 4, 3, 2, 0/))
   !primdata(51) = sum(2**(/6, 3, 1, 0/))
   !primdata(52) = sum(2**(/52, 3, 0/))
   !primdata(53) = sum(2**(/53, 6, 2, 1, 0/))
   !primdata(54) = sum(2**(/54, 6, 5, 4, 3, 2, 0/))
   !primdata(55) = sum(2**(/55, 6, 2, 1, 0/))
   !primdata(56) = sum(2**(/56, 7, 4, 2, 0/))
   !primdata(57) = sum(2**(/57, 5, 3, 2, 0/))
   !primdata(58) = sum(2**(/58, 6, 5, 1, 0/))
   !primdata(59) = sum(2**(/59, 6, 5, 4, 3, 1, 0/))
   !primdata(60) = sum(2**(/60, 1, 0/))
   !primdata(61) = sum(2**(/61, 5, 2, 1, 0/))
   !primdata(62) = sum(2**(/62, 6, 5, 3, 0/))
   !primdata(63) = sum(2**(/63, 1, 0/))
   !primdata(64) = sum(2**(/64, 4, 3, 1, 0/))
   !primdata(65) = sum(2**(/65, 4, 3, 1, 0/))
   !primdata(66) = sum(2**(/66, 8, 6, 5, 3, 2, 0/))
   !primdata(67) = sum(2**(/67, 5, 2, 1, 0/))
   !primdata(68) = sum(2**(/68, 7, 5, 1, 0/))
   !primdata(69) = sum(2**(/69, 6, 5, 2, 0/))
   !primdata(70) = sum(2**(/70, 5, 3, 1, 0/))
   !primdata(71) = sum(2**(/71, 5, 3, 1, 0/))
   !primdata(72) = sum(2**(/72, 6, 4, 3, 2, 1, 0/))
   !primdata(73) = sum(2**(/73, 4, 3, 2, 0/))
   !primdata(74) = sum(2**(/74, 7, 4, 3, 0/))
   !primdata(75) = sum(2**(/75, 6, 3, 1, 0/))
   !primdata(76) = sum(2**(/76, 5, 4, 2, 0/))
   !primdata(77) = sum(2**(/77, 6, 5, 2, 0/))
   !primdata(78) = sum(2**(/78, 7, 2, 1, 0/))
   !primdata(79) = sum(2**(/79, 4, 3, 2, 0/))
   !primdata(80) = sum(2**(/80, 7, 5, 3, 2, 1, 0/))
   !primdata(81) = sum(2**(/81, 4, 0/))
   !primdata(82) = sum(2**(/82, 8, 7, 6, 4, 1, 0/))
   !primdata(83) = sum(2**(/83, 7, 4, 2, 0/))
   !primdata(84) = sum(2**(/84, 8, 7, 5, 3, 1, 0/))
   !primdata(85) = sum(2**(/85, 8, 2, 1, 0/))
   !primdata(86) = sum(2**(/86, 6, 5, 2, 0/))
   !primdata(87) = sum(2**(/87, 7, 5, 1, 0/))
   !primdata(88) = sum(2**(/88, 8, 5, 4, 3, 1, 0/))
   !primdata(89) = sum(2**(/89, 6, 5, 3, 0/))
   !primdata(90) = sum(2**(/90, 5, 3, 2, 0/))
   !primdata(91) = sum(2**(/91, 7, 6, 5, 3, 2, 0/))
   !primdata(92) = sum(2**(/92, 6, 5, 2, 0/))
   !primdata(93) = sum(2**(/93, 2, 0/))
   !primdata(94) = sum(2**(/94, 6, 5, 1, 0/))
   !primdata(95) = sum(2**(/95, 6, 5, 4, 2, 1, 0/))
   !primdata(96) = sum(2**(/96, 7, 6, 4, 3, 2, 0/))
   !primdata(97) = sum(2**(/97, 6, 0/))
   !primdata(98) = sum(2**(/98, 7, 4, 3, 2, 1, 0/))
   !primdata(99) = sum(2**(/99, 7, 5, 4, 0/))
   !primdata(100) = sum(2**(/100, 8, 7, 2, 0/))
   ! The last two are of degree 107 and 127!!!
   !primdata(101) = sum(2**(/107, 7, 5, 3, 2, 1, 0/))
   !primdata(102) = sum(2**(/127, 1, 0/))

   polys = primdata(1:s)

  end subroutine get_primpolys_watson

end module mod_primpoly
