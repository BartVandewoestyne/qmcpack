! Module containing the main routine for the adaptive cube based QMC algorithm.
!
! Parameters we can play with
!
!    * See if we need to multiply by the hollow or the full volume
!      in the calculation of the variation estimator.
!
! References:
!
!   [1] `Accuracy and Stability of Numerical Algorithms', Higham, N. J.,
!       Society for Industrial and Applied Mathematics, Philadelphia, 2002,
!       ISBN 0-89871-355-2.
!
!   [2] `An adaptive approach to cube based quasi-Monte Carlo integration
!       on R^d', Pillards, T., Vandewoestyne, B. and Cools, R., submitted
!       for the proceedings of the MCM2007 conference.
!
module mod_adaptive_cube_qmc_rd

  use numeric_kinds
  use mod_function
  use mod_halton
  use mod_sobol

  private

  public  :: init_adaptive_cube_qmc
  public  :: adaptive_cube_qmc
  public  :: free_adaptive_cube_qmc

  private :: initialize_cubelist
  private :: initialize_cube
  private :: expand
  private :: set_points
  private :: estimate_integral
  private :: add_cube
  private :: print_cube_list
  private :: get_variation_estimator
  private :: get_total_var
  private :: add_point
  private :: get_value
  private :: get_volume
  private :: update_parameters
  private :: init_pointset
  private :: next_point
  private :: free_pointset

  ! The initial width of a hypercube (default = 2)
  real(kind=qp), private :: V = 2.0_qp

  ! The growth of the widths of the cubes (default = 2)
  real(kind=qp), private :: w = 2.0_qp

  ! If the algorithm notices that the number of points needed in a cube is lower
  ! than eps, we will not place points in this cube.
  real(kind=qp), private :: eps = 1.0_qp

  ! When above the treshold eps, we will use at least min_points points.  Using
  ! less will make the estimator too inaccurate.
  integer(kind=i4b), private :: min_points

  ! 
  integer(kind=i4b), dimension(:), allocatable, private :: N_step

  ! The dimension of the problem
  integer(kind=i4b), private :: d

  type, public :: hypercube
  
    ! the index of the `box'e
    integer(kind=i4b) :: j

    ! upper border of the hollow hypercube
    real(kind=qp)     :: u_border

    ! number of effectively used points to calculate the variation estimator
    ! (due to the characteristic function over I_j)
    integer(kind=i4b) :: nb_feval

    ! number of points assigned to this `hypercube'
    integer(kind=i4b) :: nb_points

    ! the maximum function value, calculated over all the points in this hollow
    ! hypercube
    real(kind=qp)     :: f_max

    ! the minimum function value, calculated over all the points in this hollow
    ! hypercube
    real(kind=qp)     :: f_min

    ! the sum of the absolute functionvalues at the points in this cube
    real(kind=qp)     :: abssum

    ! the sum of the function values of the (transformed) points that fall
    ! within this hollow hypercube
    real(kind=qp)     :: f_sum

    ! the sum of the squares of the function values of the (transformed) points
    ! that fall within this hollow hypercube
    real(kind=qp)     :: f2_sum

    ! needed for stable calculation of the sample variance (See [1], page 13)
    real(kind=qp)     :: m_n

    ! needed for stable calculation of the sample variance (See [1], page 13)
    real(kind=qp)     :: q_n

  end type hypercube


contains


  ! Initialize the algorithm with a certain initial width V of a
  ! hypercube and a certain initial growth w of the width of the cubes.
  ! If this initialization is forgotten or not done, then by
  ! default V=2 and w=2.
  !
  subroutine init_adaptive_cube_qmc(init_d, init_V, init_w, init_eps, &
                                    init_min_points, init_N_step)
    integer(kind=i4b), intent(in)               :: init_d
    real(kind=qp), intent(in)                   :: init_V
    real(kind=qp), intent(in)                   :: init_w
    real(kind=qp), intent(in)                   :: init_eps
    integer(kind=i4b), intent(in)               :: init_min_points
    integer(kind=i4b), dimension(:), intent(in) :: init_N_step

    d = init_d
    V = init_V
    w = init_w
    eps = init_eps
    min_points = init_min_points

    if (allocated(N_step)) then
      deallocate(N_step)
    end if
    allocate(N_step(size(init_N_step)))
    N_step = init_N_step

    !if (associated(N_step)) then
    !  N_step => null()
    !end if
    !N_step => init_N_step

  end subroutine init_adaptive_cube_qmc


  ! The main subroutine.
  !
  subroutine adaptive_cube_qmc(f,              &
                               params,         &
                               variation_type, &
                               point_family,   &
                               res)

    interface
      function f(x, fparams) result (y)
        use mod_function
        use numeric_kinds
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function f
    end interface
    type(functionparams), intent(in)            :: params
    character(len=*), intent(in)                :: variation_type, point_family
    real(kind=qp), dimension(:), intent(out)    :: res 

    integer(kind=i4b)                      :: i
    integer(kind=i4b)                      :: n_current
    type(hypercube), dimension(:), pointer :: cubelist

    call initialize_cubelist(f, params, cubelist, point_family)
    !call print_cube_list(cubelist)

    do i=1,size(N_step)

      n_current = N_step(i)
      write(unit=*, fmt="(A)") ""
      write(unit=*, fmt="(A, I0, A, I0)") &
                        "==> Running test for N_step(", i, ") = ", n_current

      call expand(f, params, cubelist, variation_type, n_current, point_family)
      call set_points(f, params, cubelist, variation_type, n_current, point_family)
      !call print_cube_list(cubelist)
      call estimate_integral(cubelist, res(i))

    end do

    !do i=1, size(res)
    !  write(unit=*, fmt="(A, F25.15)") "I = ", res(i)
    !end do

    deallocate(cubelist)

  end subroutine adaptive_cube_qmc


  ! Initialize the list of hypercubes with two cubes containing min_points
  ! points per hypercube.
  !
  subroutine initialize_cubelist(f, params, list, point_family)

    interface
      function f(x, fparams) result (y)
        !use qmcpack
        use mod_function
        use numeric_kinds
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function f
    end interface
    type(functionparams), intent(in)       :: params
    type(hypercube), dimension(:), pointer :: list
    character(len=*), intent(in)           :: point_family

    type(hypercube) :: cube1, cube2

    write(unit=*, fmt="(A)", advance="no") "Initializing cube-list... "

    allocate(list(2))

    call initialize_cube(cube1, 0)
    call update_parameters(f, params, cube1, min_points, point_family)

    call initialize_cube(cube2, 1)
    call update_parameters(f, params, cube2, min_points, point_family)

    list(1) = cube1
    list(2) = cube2

    write(unit=*, fmt="(A)") "done."

  end subroutine initialize_cubelist


  ! Add more cubes to the list of hypercubes if necessary.  When the algorithm
  ! calculates that the number of points needed in a cube is lower than
  ! eps, we will not place points in this cube.
  !
  subroutine expand(f,              &
                    params,         &
                    cubelist,       &
                    variation_type, &
                    n_current,      &
                    point_family)

    interface
      function f(x, fparams) result (y)
        use mod_function
        use numeric_kinds
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function f
    end interface
    type(functionparams), intent(in)       :: params
    type(hypercube), dimension(:), pointer :: cubelist
    real(kind=qp)                          :: total_var
    character(len=*), intent(in)           :: variation_type, point_family
    integer(kind=i4b), intent(in)          :: n_current

    type(hypercube) :: newcube

      
    write(unit=*, fmt="(A)") "Expanding cube-list if necessary... "

    total_var = get_total_var(cubelist, variation_type)
    do

      ! Keep on adding `cubes' until the number of points needed in a cube is
      ! lower than eps.
      !print *, "DEBUG: variation estimator for current last cube = ", &
      !  get_variation_estimator((cubelist(size(cubelist))), variation_type)
      !print *, "DEBUG: relative, this is = ", &
      !  get_variation_estimator((cubelist(size(cubelist))), variation_type)/total_var
      if (get_variation_estimator((cubelist(size(cubelist))), variation_type)/total_var*n_current > eps) then

        write(unit=*, fmt="(A)", advance="no") "  Adding extra cube to cube-list... " 

        call initialize_cube(newcube, size(cubelist))
        call update_parameters(f, params, newcube, min_points, point_family)

        call add_cube(cubelist, newcube)
        total_var = total_var &
            + get_variation_estimator(cubelist(size(cubelist)), variation_type)

        write(unit=*, fmt="(A)") "done."

      else
        exit
      end if

    end do

    write(unit=*, fmt="(A)") "done."

  end subroutine expand


  ! Return the sum of all variations over all the hypercubes.
  !
  function get_total_var(cubelist, variation_type) result (total_var)
    type(hypercube), dimension(:), pointer :: cubelist
    real(kind=qp)                          :: total_var
    character(len=*), intent(in)           :: variation_type
     
    integer(kind=i4b) :: i

    total_var = 0.0_qp
    do i=1,size(cubelist)
      total_var = total_var &
                        + get_variation_estimator(cubelist(i), variation_type)
    end do

  end function get_total_var


  ! Determine and set the amount of points for each hypercube, based on the
  ! variation information.
  !
  subroutine set_points(f, params, cube_list, variation_type, n_current, point_family)

    interface
      function f(x, fparams) result (y)
        use mod_function
        use numeric_kinds
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function f
    end interface
    type(functionparams), intent(in)       :: params
    type(hypercube), dimension(:), pointer :: cube_list
    character(len=*), intent(in)           :: variation_type, point_family
    integer(kind=i4b), intent(in)          :: n_current

    integer(kind=i4b) :: i
    real(kind=qp)     :: this_var
    real(kind=qp)     :: total_var
    integer(kind=i4b) :: total_nb_points

    write(unit=*, fmt="(A)") "Calculating and setting the "// &
        "amount of points for each hypercube from the cube-list... "

    total_nb_points = 0
    total_var = get_total_var(cube_list, variation_type)

    do i=1, size(cube_list)

      this_var = get_variation_estimator(cube_list(i), variation_type)

      call update_parameters(f, params, cube_list(i), &
                            ceiling(n_current*this_var/total_var), point_family)
      write(unit=*, fmt="(A, I3, A, I10, A)") &
                        "  Cube ", i, ": ", cube_list(i)%nb_points, " points."
      total_nb_points = total_nb_points + cube_list(i)%nb_points

    end do

    write(unit=*, fmt="(A)") "  ----------------------------"
    write(unit=*, fmt="(A, I10, A, I0, A)") &
      "  Total   : ", total_nb_points, " points, so ", &
      total_nb_points - n_current, " extra points were added."
    write(unit=*, fmt="(A)") "done."

  end subroutine set_points


  ! After determining the number of points for each hypercube, we can
  ! estimate the integral by iterating over the cube list.
  !
  subroutine estimate_integral(cubelist, res)
    type(hypercube), dimension(:), pointer :: cubelist
    real(kind=qp), intent(out)             :: res

    integer(kind=i4b) :: j

    write(unit=*, fmt="(A)") "Calculating integral..."

    res = 0.0_qp
    do j=1,size(cubelist)
      res = res + get_value(cubelist(j))

      !print *, "j = ", j-1, "    value = ", get_value(cubelist(j))
      !print *, "cube%f_sum = ", cubelist(j)%f_sum
      !print *, "cube%nb_feval = ", cubelist(j)%nb_feval
      !print *, "volume(cube) = ", get_volume(cubelist(j), "full")

    end do

    write(unit=*, fmt="(A, F25.15)") "  INTEGRAL APPROXIMATION = ", res
    write(unit=*, fmt="(A)") "done."

  end subroutine estimate_integral


  ! Add a hypercube to the list of hypercubes.
  !
  subroutine add_cube(list, cube)
    type(hypercube), dimension(:), pointer :: list
    type(hypercube), intent(in)            :: cube

    integer(kind=i4b)                        :: k
    type(hypercube), dimension(size(list)+1) :: x

    k = size(list)
    x(1:k) = list
    x(k+1) = cube
    deallocate(list)
    allocate(list(k+1))
    list = x

  end subroutine add_cube


  ! Initialize the j-th d-dimensional cube.
  !
  subroutine initialize_cube(cube, j)
    type(hypercube), intent(inout) :: cube
    integer(kind=i4b), intent(in)  :: j

    cube%j         = j
    cube%u_border  = 0.5_qp*w**(cube%j)*V
    cube%nb_feval    = 0
    cube%nb_points = 0
    cube%f_sum     = 0.0_qp
    cube%f2_sum    = 0.0_qp
    cube%m_n       = 0.0_qp
    cube%q_n       = 0.0_qp
    cube%abssum    = 0.0_qp
    cube%f_min     = huge(1.0_qp)
    cube%f_max     = -huge(1.0_qp)

  end subroutine initialize_cube


  ! Print the cube list.
  !
  subroutine print_cube_list(list)
    type(hypercube), dimension(:), pointer :: list

    integer(kind=i4b) :: i

    write(unit=*, fmt="(A)") "Cubelist now looks as follows:"
    do i=1,size(list)
      write(unit=*, fmt="(A, I2, A, I4, A, I8, A, ES25.15, A)") &
        "  (j=", list(i)%j, ", d = ", d, ", nb_points = ",   &
        list(i)%nb_points, ", volume = ", get_volume(list(i), "hollow"), ")"
    end do

  end subroutine print_cube_list


  ! Return the given variation of the given variation type for this
  ! hypercube.
  !
  function get_variation_estimator(cube, variation_type) result (res)
    type(hypercube), intent(in)  :: cube
    character(len=*), intent(in) :: variation_type
    real(kind=qp) :: res

    real(kind=qp) :: s

    select case (variation_type)
      
      case ("mathewei.5")

        s = d + 0.5_qp

        ! Wat Tim had
        res = (2*cube%u_border)**(-s)

        ! Wat Bart denkt dat het moet zijn (in geval V=2, w=2)
        !res = w**(-cube%j*(0.5-d))
        !res = 1.0_qp/(V**d*w**(cube%j*s)*(w**d-1.0_qp))

      case ("mathewei1")

        s = d + 1.0_qp

        ! Wat Tim had
        res = (2*cube%u_border)**(-s)

        ! Wat Bart denkt dat het moet zijn (in geval V=2, w=2)
        !res = w**(-cube%j*(1-d))
        !res = 1.0_qp/(V**d*w**(cube%j*s)*(w**d-1.0_qp))
        
      case ("mathewei2")

        s = d + 2.0_qp

        ! Wat Tim had
        res = (2*cube%u_border)**(-s)

        ! Wat Bart denkt dat het moet zijn (in geval V=2, w=2)
        !res = w**(-cube%j*(2-d))
        !res = 1.0_qp/(V**d*w**(cube%j*s)*(w**d-1.0_qp))
        
      case ("maxmin")
        res = cube%f_max-cube%f_min

      case ("sqrtmaxmin")
        res = sqrt(cube%f_max-cube%f_min)

      case ("avgabssum")
        res = cube%abssum/cube%nb_feval

      case ("avgabssum2")
        res = (cube%abssum/cube%nb_feval)**2.0_qp

      case ("variance")

        ! This way of calculating the variance suffers from severe roundoff
        ! errors...
        !res = (cube%f2_sum - (cube%f_sum)**2.0_qp/cube%nb_feval)/(cube%nb_feval-1.0_qp)

        ! See [1] on page 13 on how we compute the variance with a one-pass
        ! algorithm in a way with less roundoff errors.
        res = cube%q_n/(cube%nb_feval-1.0_qp)

      case ("stddev")

        ! This way of calculating the standard deviation suffers from severe
        ! roundoff errors...
        !res = sqrt(cube%f2_sum-cube%f_sum)

        ! See [1] on page 13 on how we compute the variance with a one-pass
        ! algorithm in a way with less roundoff errors.
        res = sqrt(cube%q_n/(cube%nb_feval-1.0_qp))

      case ("sqrtstddev")
        res = sqrt(sqrt(cube%q_n/(cube%nb_feval-1.0_qp)))

      case default

        print *, "WARNING: unknown variation type, using `standard deviation' as " &
               // "variation type!"
        res = sqrt(cube%q_n/(cube%nb_feval-1.0_qp))

    end select


    select case (variation_type)

      case ("mathewei1", "mathewei2", "mathewei.5")

        ! Wat Tim had
        ! TODO:
        !   * check of dit niet "full" moet zijn!
        !     => Blijkbaar maakt het niet al te veel uit voor de grafieken.
        !        Check waarom het niet uitmaakt!
        !     => Eventueel kan je het ook anders oplossen door hier
        !        gewoon res = res*1 te zetten en in de bovenstaande code dan
        !        de juiste fractie te berekenen.
        res = res*get_volume(cube, "hollow")

        ! Wat Bart denkt dat het moet zijn...
        !res = res*get_volume(cube, "full")
        !res = res*1.0_qp

      case default

        ! Wat Tim had
        ! TODO:
        !   * check of dit niet "full" moet zijn?  Zie midden paper
        !     Tim pag. 11.
        !     => Blijkbaar maakt het niet al te veel uit voor de grafieken.
        !        Check waarom het niet veel uitmaakt!
        res = sqrt(res*get_volume(cube, "hollow"))

        ! Wat Bart denkt dat het moet zijn...
        !res = sqrt(res*get_volume(cube, "full"))

    end select
        
  end function get_variation_estimator


  ! Return the volume of this `hypercube'.  The volume_type can be:
  !
  !  "full"   => the volume of a full hypercube.
  !  "hollow" => the volume of a hollow hypercube.
  !
  function get_volume(cube, volume_type) result (vol)
    type(hypercube), intent(in)  :: cube
    character(len=*), intent(in) :: volume_type
    real(kind=qp)                :: vol

    select case (volume_type)

      case ("full")

        vol = (V*w**cube%j)**d

      case ("hollow")

        if (cube%j == 0) then
          vol = (V*w**cube%j)**d
        else
          vol = ((V*w**cube%j)**d)*(w**d-1.0_qp)
        end if

      case default

        print *, "ERROR: wrong volume type!"

    end select

  end function get_volume


  ! Adds the following point from the discrepancy sequence to the box and
  ! updates all variables.
  ! (timp: this was not used, i added it now to update_parameters...)
  !
  subroutine add_point(f, params, cube, x)

    interface
      function f(x, fparams) result (y)
        !use qmcpack
        use mod_function
        use numeric_kinds
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function f
    end interface
    type(functionparams), intent(in)        :: params
    type(hypercube), intent(inout)          :: cube
    real(kind=qp), dimension(:), intent(in) :: x

    real(kind=qp) :: fx


    cube%nb_points = cube%nb_points + 1

    ! If the point falls within this hollow hypercube...
    ! TODO:
    !   * Check if we should also check that maxval(abs(x)) < cube%u_border
    if (maxval(abs(x)) > 0.5_qp*w**(cube%j-1.0_qp)*V .or. cube%j==0) then
      fx = f(x, params)

      ! ... update the characteristics that we have to keep track of for the
      ! variation estimation.
      ! (See also MWAVariation.m)
      cube%nb_feval = cube%nb_feval + 1
      cube%f_max = max(cube%f_max, fx)
      cube%f_min = min(cube%f_min, fx)
      cube%abssum = cube%abssum + abs(fx)
      cube%f_sum = cube%f_sum + fx
      cube%f2_sum = cube%f2_sum + fx**2

      ! See [1], page 13 on how the variance is calculated in a way without
      ! too much roundoff error.
      cube%q_n = cube%q_n + ((cube%nb_feval-1)*(fx-cube%m_n)**2)/cube%nb_feval
      cube%m_n = cube%m_n + (fx-cube%m_n)/cube%nb_feval

    end if
      
  end subroutine add_point


  ! Returns an estimate for f on this hypercube.
  !
  ! TODO:
  !   * Check if we need cube%nb_points and if we need 'full' for the volume.
  !
  function get_value(cube) result (res)
    type(hypercube), intent(in) :: cube
    real(kind=qp)               :: res

    res = cube%f_sum/cube%nb_points*get_volume(cube, "full")

  end function get_value


  ! Update the variation-parameters for the specified cube, using the
  ! specified function and evaluating for nb_points points.
  !
  subroutine update_parameters(f, params, cube, nb_points, point_family)

    interface
      function f(x, fparams) result (y)
        !use qmcpack
        use mod_function
        use numeric_kinds
        real(kind=qp), dimension(:), intent(in) :: x
        type(functionparams), intent(in)        :: fparams
        real(kind=qp)                           :: y
      end function f
    end interface
    type(functionparams), intent(in) :: params
    type(hypercube), intent(inout)   :: cube
    integer(kind=i4b), intent(in)    :: nb_points
    character(len=*), intent(in)     :: point_family

    integer(kind=i4b)                    :: i
    real(kind=qp), dimension(d)     :: x
    integer(kind=i4b), dimension(d) :: start_index

    if (nb_points > cube%nb_points) then

      start_index = cube%nb_points + 1  
      call init_pointset(start_index, nb_points, point_family)

      do i = cube%nb_points+1, nb_points

        ! Get next (quasi)-random point.
        call next_point(x, point_family)

        ! Scale x
        ! TODO: x werd door Tim verkeerd gescaleerd!
        ! Dit geeft de resultaten uit de paper.
        !x = (x-0.5_qp*2.0_qp)*cube%u_border
        ! Dit is de juiste scalering maar geeft rare resultaten???
        x = (x-0.5_qp)*V*(w**(cube%j))

        if (any(x>cube%u_border) .or. any(x<-cube%u_border)) then
          stop "ERROR: x is being scaled wrong in update_parameters(...)!"
        end if

        call add_point(f, params, cube, x)

      end do

      call free_pointset(point_family)

    end if

  end subroutine update_parameters
    
    
  subroutine init_pointset(start_index, max_points, point_family)

    ! initialises the point set
    ! adding new point set must be done here
    !   but also in next_point and free_pointset
    ! to add parameters, add special case e.g. sobol2
      
    integer(kind=i4b), intent(in)               :: max_points
    integer(kind=i4b), dimension(:), intent(in) :: start_index
    character(len=*), intent(in)                :: point_family
    type(soboltype)                             :: mysoboltype
    
    select case (point_family)

      case ("random", "pseudorandom", "mc")

      case ("halton", "quasirandom", "qmc")

        call init_halton(d, start_index)

      case ("sobol_BratleyFox")
    
        mysoboltype%poly_order = "BratleyFox"
        mysoboltype%dirnumbers = "BratleyFox"
        mysoboltype%use_ones = .true.
        mysoboltype%use_antonov_saleev = .true.
        call init_sobol(max_points, d, start_index,mysoboltype)

      case ("sobol_NumericalRecipes")
    
        mysoboltype%poly_order = "NumericalRecipes"
        mysoboltype%dirnumbers = "NumericalRecipes"
        mysoboltype%use_ones = .false.
        mysoboltype%use_antonov_saleev = .true.
        call init_sobol(max_points, d, start_index,mysoboltype)

      case ("sobol","sobol_JoeKuo")

        mysoboltype%poly_order = "JoeKuo"
        mysoboltype%dirnumbers = "JoeKuo"
        mysoboltype%use_ones = .true.
        mysoboltype%use_antonov_saleev = .true.
        call init_sobol(max_points, d, start_index,mysoboltype)

      case default

        call init_halton(d, start_index)

        print *, "WARNING: using default point set: Halton!"

    end select
    
  end subroutine init_pointset
    
    
  ! Returns the next point from the sequence.  Note that
  ! the sequence must have been initiated first!
  !
  subroutine next_point(x, point_family)
    real(kind=qp), dimension(:), intent(inout) :: x
    character(len=*), intent(in)               :: point_family
    
    select case (point_family)

        case ("random", "pseudorandom", "mc")

          call random_number(x)

        case ("halton", "quasirandom", "qmc")

          call next_halton(x)

        case ("sobol")
  
          call next_sobol(x)

        case default

          call next_halton(x)

      end select
    
    end subroutine next_point
    
        
    ! Un-initializes the pointset.
    !
    subroutine free_pointset(point_family)
      character(len=*), intent(in) :: point_family

      select case (point_family)

        case ("random", "pseudorandom", "mc")

        case ("halton", "quasirandom", "qmc")

          call free_halton()

        case ("sobol")
  
          !call free_sobol()

        case default

          call free_halton()

      end select
    
    end subroutine free_pointset


    subroutine free_adaptive_cube_qmc()

      if (allocated(N_step)) then
        deallocate(N_step)
      end if

      !N_step => null()

    end subroutine free_adaptive_cube_qmc
    
    
end module mod_adaptive_cube_qmc_rd
