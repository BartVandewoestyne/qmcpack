! Program to test and try to reproduce Tims adaptive qMC results (that were
! written in Matlab).
!
! Parameters we can play with:
!
!       * width of the initial cube V
!       * multiplier w
!       * the epsilon value
!       * the minimum number of points per box
!       * the startindex of the low discrepancy pointset
!
! Notes:
!   * Tim did runs up to N=2^28.
!
! References:
!
!   [1] `An adaptive approach to cube based quasi-Monte Carlo integration
!       on R^d', Pillards, T., Vandewoestyne, B. and Cools, R., submitted
!       for the proceedings of the MCM2007 conference.
!
! TODO:
!   * check what happens with w=2^(1/d) instead of w=2^(2/d)
!     => I think Tim's results were with w=2^(2/d), because with this value
!        for w the results are most similar... those with w=2^(1/d) slightly
!        differ...
!
program pub_adaptive_cube_qmc_rd1

  use qmcpack

  integer(kind=i4b)                           :: i, j
  integer(kind=i4b), dimension(18), parameter :: n_steps = (/ (2**i, i=1,18) /)
  real(kind=qp), dimension(size(n_steps))     :: res
  type(functionparams)                        :: myparams
  integer(kind=i4b)                           :: s, nb_columns
  integer(kind=i4b)                           :: min_points
  real(kind=qp)                               :: my_epsilon
  real(kind=qp)                               :: V, w
  real(kind=qp)                               :: exact
  real(kind=qp), dimension(:,:), allocatable  :: rel_err
  character(len=6)                            :: point_family = "halton"
  character(len=45)                           :: file_name
  character(len=*), parameter, dimension(1:6) :: var_estimator = &
                                                        (/ "mathewei1", &
                                                           "mathewei2", &
                                                           "variance ", &
                                                           "stddev   ", &
                                                           "maxmin   ", &
                                                           "avgabssum" /)

  ! PARAMETERS CHOSEN FOR THE ALGORITHM

  ! The epsilon.
  ! NOTE: i think this is 0.5 in Tim's Matlab code.  See MWABoxes.m???
  my_epsilon = 1.0_qp
  
  ! The initial width.
  V = 2.0_qp

  ! The minimum number of points for a cube.
  min_points = 100

  ! END OF THE PART WITH PARAMETERS FOR THE ALGORITHM

  
  ! Allocate the array that we're going to use to write the relative errors to
  ! a file.  The columns represent the different methods for variation
  ! estimation.
  nb_columns = size(var_estimator)
  allocate(rel_err(size(res), nb_columns))


  ! TIMS FIRST TESTFUNCTION IN 10 DIMENSIONS

  ! The dimension of the problem.
  s = 10
  w = 2.0_qp**(2.0_qp/s)
  !w = 2.0_qp
  call init_adaptive_cube_qmc(s, V, w, my_epsilon, min_points, n_steps)

  ! Calculate the exact value of the integral in s dimensions.
  exact = f_timps_cube01_exact_inf(s)

  call create_filename(s, "f_timps_cube01", point_family, file_name)

  do j=1,nb_columns

    print *, "Method: "//trim(var_estimator(j))
    call adaptive_cube_qmc(f_timps_cube01_1p,      &
                           myparams,               &
                           trim(var_estimator(j)), &
                           point_family,           &
                           res)
     rel_err(:,j) = abs(exact-res)/exact

  end do

  print *, "EXACT VALUE = ", exact
  call write_to_file(rel_err, file_name)


  ! TIMS FIRST TESTFUNCTION IN 15 DIMENSIONS

  ! The dimension of the problem.
  s = 15
  w = 2.0_qp**(2.0_qp/s)
  !w = 2.0_qp
  call init_adaptive_cube_qmc(s, V, w, my_epsilon, min_points, n_steps)

  ! Calculate the exact value of the integral in s dimensions.
  exact = f_timps_cube01_exact_inf(s)

  call create_filename(s, "f_timps_cube01", point_family, file_name)

  do j=1,nb_columns

    print *, "Method: "//trim(var_estimator(j))
    call adaptive_cube_qmc(f_timps_cube01_1p,      &
                           myparams,               &
                           trim(var_estimator(j)), &
                           point_family,           &
                           res)
     rel_err(:,j) = abs(exact-res)/exact

  end do

  print *, "EXACT VALUE = ", exact
  call write_to_file(rel_err, file_name)
  call free_adaptive_cube_qmc()


  ! TIMS SECOND TESTFUNCTION IN 10 DIMENSIONS

  ! The dimension of the problem.
  s = 10
  w = 2.0_qp**(2.0_qp/s)
  !w = 2.0_qp
  call init_adaptive_cube_qmc(s, V, w, my_epsilon, min_points, n_steps)

  ! Calculate the exact value of the integral in s dimensions.
  exact = f_timps_cube02_exact_inf(s)

  call create_filename(s, "f_timps_cube02", point_family, file_name)

  do j=1,nb_columns

    print *, "Method: "//trim(var_estimator(j))
    call adaptive_cube_qmc(f_timps_cube02_1p,      &
                           myparams,               &
                           trim(var_estimator(j)), &
                           point_family,           &
                           res)
     rel_err(:,j) = abs(exact-res)/exact

  end do

  print *, "EXACT VALUE = ", exact
  call write_to_file(rel_err, file_name)
  call free_adaptive_cube_qmc()


  ! TIMS SECOND TESTFUNCTION IN 15 DIMENSIONS

  ! The dimension of the problem.
  s = 15
  w = 2.0_qp**(2.0_qp/s)
  !w = 2.0_qp
  call init_adaptive_cube_qmc(s, V, w, my_epsilon, min_points, n_steps)

  ! Calculate the exact value of the integral in s dimensions.
  exact = f_timps_cube02_exact_inf(s)

  call create_filename(s, "f_timps_cube02", point_family, file_name)

  do j=1,nb_columns

    print *, "Method: "//trim(var_estimator(j))
    call adaptive_cube_qmc(f_timps_cube02_1p,      &
                           myparams,               &
                           trim(var_estimator(j)), &
                           point_family,           &
                           res)
     rel_err(:,j) = abs(exact-res)/exact

  end do

  print *, "EXACT VALUE = ", exact
  call write_to_file(rel_err, file_name)
  call free_adaptive_cube_qmc()

  deallocate(rel_err)

contains

  subroutine create_filename(s, function_name, point_family, file_name)
    integer(kind=i4b), intent(in) :: s
    character(len=*), intent(in)  :: function_name
    character(len=*), intent(in)  :: point_family

    character(len=*), intent(out) :: file_name

    write(unit=file_name, fmt="(A,I2.2,A)") &
      "adaptive_qmc_dim", s, "_"//trim(function_name)//"_"//trim(point_family)// ".dat"

  end subroutine create_filename


  subroutine write_to_file(rel_err, filename)
    real(kind=qp), dimension(:,:), intent(in) :: rel_err
    character(len=*), intent(in)              :: filename

    integer(kind=i4b) :: my_unit
    integer(kind=i4b) :: j
    character(len=20) :: fmt_string

    write(unit=fmt_string, fmt="(A, I0.0, A)") "(I9, ", nb_columns, "ES25.15e4)"
   
    call get_unit(my_unit)
    call checked_open(my_unit, filename, "write")
    write(unit=my_unit, fmt="(A1, A8, 6A25)") "#", "N", var_estimator
    do j = 1,size(rel_err, dim=1)
      write(unit=my_unit, fmt=fmt_string) n_steps(j), rel_err(j,:)
    end do
    call checked_close(my_unit)
    
  end subroutine write_to_file
  
end program pub_adaptive_cube_qmc_rd1
