! Module that implements extensible lattice sequences.
!
! References:
!
!   [1] `Extensible lattice sequences for quasi-Monte Carlo quadrature',
!       Hickernell, Fred J., Hong, Hee Sun, L'Ecuyer, Pierre, Lemieux
!       Christiane, SIAM Journal on Scientific Computing, vol. 22, 2002,
!       pages 1117--1138.
!
!   [2] `Constructing embedded lattice rules for multivariate integration',
!       Cools R., Kuo, Frances Y., Nuyens D.

module mod_extensible_lattice

  use numeric_kinds
  use mod_utilities
  use mod_radical_inverse

  private

  public :: init_extensible_lattice
  public :: next_extensible_lattice
  public :: next_extensible_lattice_gray
  public :: free_extensible_lattice
  public :: get_vector_nuyens_kuo

  ! The dimension of the extensible lattice sequence
  integer(kind=i4b), private                            :: ext_lat_s

  ! The generating vector
  integer(kind=i4b), dimension(:), allocatable, private :: ext_lat_z

  ! The step by which we increment the index for the radical inverse function
  integer(kind=i4b), dimension(:), allocatable, private :: ext_lat_step

  ! The index for the radical inverse function
  integer(kind=i4b), dimension(:), allocatable, private :: ext_lat_n

  ! The bases to be used for the radical inverse function
  integer(kind=i4b), dimension(:), allocatable, private :: ext_lat_b


  contains
  

    ! Do some initialization for the extensible lattice code.
    !
    ! Defaults values:
    !   b          = 2
    !   startindex = 0
    !   step       = 1
    !
    subroutine init_extensible_lattice(s, z, b, startindex, step)

      integer(kind=i4b), intent(in)                         :: s
      integer(kind=i4b), dimension(:), intent(in)           :: z
      integer(kind=i4b), dimension(:), intent(in), optional :: b
      integer(kind=i4b), dimension(:), intent(in), optional :: startindex
      integer(kind=i4b), dimension(:), intent(in), optional :: step


      ext_lat_s = s

      if (allocated(ext_lat_z)) then
        deallocate(ext_lat_z)
      end if
      allocate(ext_lat_z(s))
      ext_lat_z = z


      if (allocated(ext_lat_b)) then
        deallocate(ext_lat_b)
      end if
      allocate(ext_lat_b(s))
      if (present(b)) then
        ext_lat_b = b
      else
        ext_lat_b = 2
      end if


      if (allocated(ext_lat_n)) then
        deallocate(ext_lat_n)
      end if
      allocate(ext_lat_n(s))
      if (present(startindex)) then
        ext_lat_n = startindex
      else
        ext_lat_n = 0
      end if


      if (allocated(ext_lat_step)) then
        deallocate(ext_lat_step)
      end if
      allocate(ext_lat_step(s))
      if (present(step)) then
        ext_lat_step = step
      else
        ext_lat_step = 1
      end if

    end subroutine init_extensible_lattice


    ! Return the next point of the extensible lattice sequence.
    !
    subroutine next_extensible_lattice(x)

      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b) :: i

      do i=1,ext_lat_s
        x(i) = frac_part(radical_inverse(ext_lat_n(i), ext_lat_b(i))*ext_lat_z(i))
      end do

      ext_lat_n = ext_lat_n + ext_lat_step

    end subroutine next_extensible_lattice


    ! Return the next point of the extensible lattice sequence, using Graycode
    ! ordering for generating the points.
    !
    ! Note:
    !   Currently, we only use base 2 graycode ordering, but it shouldn't
    !   be much of a problem to generate general base b orderings.  This will
    !   be less efficient though...
    !
    subroutine next_extensible_lattice_gray(x)

      real(kind=qp), dimension(:), intent(out) :: x

      integer(kind=i4b) :: i
      !integer(kind=i4b) :: c

      do i=1,ext_lat_s


        ! The clever way of implementing it... look for the index
        ! of the bit to change and change the bit in order to update
        ! the value of the radical inverse.

        ! Step 1: find the rightmost zero bit in the previous value for k
        !c = find_rightmost_zerobitpos(ext_lat_n(i)-1)

        ! Step 2: change that bit in the value for the radical inverse

        ! Step 3: using the updated value of the radical inverse, calculate
        !         the new value of the extensible lattice coordinate
        !x(i) = frac_part(ext_lat_phi(i)*ext_lat_z(i))



        ! The easiest way to implement it... a little less efficient...
        x(i) = frac_part(radical_inverse(graycode(ext_lat_n(i)), &
                                         ext_lat_b(i)            )*ext_lat_z(i))


      end do

      ext_lat_n = ext_lat_n + ext_lat_step

    end subroutine next_extensible_lattice_gray


    subroutine free_extensible_lattice()

      if (allocated(ext_lat_z)) then
        deallocate(ext_lat_z)
      end if
      if (allocated(ext_lat_step)) then
        deallocate(ext_lat_step)
      end if
      if (allocated(ext_lat_n)) then
        deallocate(ext_lat_n)
      end if
      if (allocated(ext_lat_b)) then
        deallocate(ext_lat_b)
      end if

    end subroutine free_extensible_lattice


    ! Return the first s generating vectors given by Nuyens and Kuo's in [2].
    !
    subroutine get_vector_nuyens_kuo(s, z)
      integer(kind=i4b), intent(in)                :: s
      integer(kind=i4b), dimension(:), intent(out) :: z

      integer(kind=i4b), dimension(360) :: genvec

      if (s>360) then
        stop "ERROR: Nuyens and Kuo only specify the first 360 dimensions!"
      end if

      genvec(1:180) = &
        (/      1, 182667, 302247, 433461, 160317,  94461, 481331, 252345, &
           358305, 221771,  48157, 489023, 438503, 399693, 200585, 169833, &
           308325, 247437, 281713, 424209, 244841, 205461, 336811, 359375, &
            86263, 370621, 422443, 284811, 231547, 360239, 505287, 355195, &
            52937, 344561, 286935, 312429, 513879, 171905,  50603, 441451, &
           164379, 139609, 371213, 152351, 138607, 441127, 157037, 510073, &
           281681, 380297, 208143, 497641, 482925, 233389, 238553, 121499, &
           137783, 463115, 168681,  70699, 390955, 202915, 230815, 304301, &
           452315, 383101,  70519, 346393, 291605, 397773,  92713, 391775, &
           131613,  74351, 382127, 219343, 297125,  88545,  89837, 191863, &
           506647, 441649, 240063, 239067, 310875, 211625, 147791, 370849, &
           149445, 340329, 493031, 336897, 202595, 518247, 369599,  50453, &
           133655, 395941,  13871, 248123, 380725, 314879, 106903,  94505, &
           499197, 479671,  76553, 351609,  13815, 342069, 446599, 374429, &
           342011,   8365, 424079, 184381, 203853,   7937, 302301, 444747, &
           151423, 492123, 166279, 335461, 295093,  98821, 118331, 371123, &
           166909, 369553, 310959, 130595, 417025, 264103, 453203, 497985, &
           200509, 269743, 461919,  45927, 394663, 299155,  81299, 112403, &
           447473, 480325, 105053, 328455, 513239, 322199,  47537, 485183, &
           490687, 214311,  61871, 359761, 509981, 192829,  17075, 486463, & 
           463461, 438237, 444353, 381279, 425405, 179723, 203897, 300411, &
            35087, 501519, 172275, 278507, 473739, 141429, 117025, 516465, &
           204743, 505099, 359937, 359503 /)
       genvec(181:360) =                (/ 455605, 227317, 251285, 295275, &
           476155, 291155, 294443,  31913, 518445, 481917, 326981, 488711, &
           447541,  39629, 394681, 379411, 335309,  93541, 188491, 152371, &
           408829, 217957,  65145, 462013, 220627, 368989, 273585, 373297, &
           285793, 405807,  63693, 382573,  84291, 444801, 226471, 166195, &
           170689, 368423, 509169, 243221, 191447, 200867, 194633, 226469, &
           413665, 148007, 284505, 459795, 324557, 422149, 368087, 166133, &
           401441,  44077, 457535, 122611,  27489, 392389, 327231,  87841, &
             8601, 118219, 503827, 368993, 439439, 153211, 282303, 252999, &
            28617, 514327, 422201, 465011, 160907, 191687, 226205, 254941, &
           480375, 322857, 456621, 446175, 109049, 285305, 436731, 290763, &
           290439,  25363, 480371, 153337, 406899,  90863,  78537, 358757, &
            69087, 431749, 384083, 472419, 298375, 291499, 187231,   6967, &
           338357, 411965, 447817, 135463, 156061, 356943, 224431, 452203, &
           284195, 489071, 299537, 313173, 436265, 316447,   8353, 512723, &
           522699, 453557, 512829, 214315,  83785,  11217, 163315, 293397, &
           386597, 338761,  88653, 337581, 230703, 140519, 421913, 416151, &
           197025, 350607, 262579, 510879, 451713,   8619,  50451, 193793, &
           155739,  23611, 437689, 100267, 439671, 211341,  22951,  73405, &
           220037, 169733, 408633,  64171,  44141, 122149, 258491, 215021, &
            11507, 295341,  41981, 378867, 261545,  56423, 445605, 362017, &
           504897, 335763, 417407,  58033,  80239, 372353, 116163, 287617, &
           483795, 310473,  75779, 141943, 418435, 305609, 152901, 129869 /)

      z = genvec(1:s)

    end subroutine get_vector_nuyens_kuo

end module mod_extensible_lattice
