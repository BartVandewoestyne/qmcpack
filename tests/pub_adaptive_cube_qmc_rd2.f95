! This program tries to reproduce Tim's Mathe and Wei results appearing in [1].
! Tim did his experiments in Matlab.  This code checks his results with our
! Fortran implementation of Mathe and Wei's algorithm.
!
! References:
!
!   [1] `An adaptive approach to cube based quasi-Monte Carlo integration
!       on R^d', Pillards, T., Vandewoestyne, B. and Cools, R., submitted
!       for the proceedings of the MCM2007 conference.
!
! For his Mathe-Wei tests, Tim used dimension+[0.5, 1, 2] as decay factors.
!
program pub_adaptive_cube_qmc_rd2

  use qmcpack

  integer(kind=i4b)                            :: d, s, max_size
  integer(kind=i4b)                            :: i, j, k
  real(kind=qp), dimension(:,:), allocatable   :: x
  integer(kind=i4b), dimension(2), parameter   :: mw_factors = (/ 1, 2 /)
  integer(kind=i4b), dimension(18), parameter  :: n_steps = (/ (2**i, i=1,18) /)
  !integer(kind=i4b), dimension(1), parameter  :: n_steps = (/ 10 /)
  real(kind=qp)                                :: res
  real(kind=qp), dimension(:,:), allocatable   :: rel_err
  real(kind=qp)                                :: V, w
  real(kind=qp)                                :: exact
  character(len=45)                            :: file_name
  integer(kind=i4b), dimension(:), allocatable :: startindex
  type(functionparams)                         :: dummy_params


  ! The initial width.
  !V = 1.0_qp
  V = 2.0_qp

  allocate(rel_err(size(n_steps), size(mw_factors)))


  ! TIMS FIRST TESTFUNCTION IN 10 DIMENSIONS
  s = 10
  exact = f_timps_cube01_exact_inf(s)

  do k = 1,size(mw_factors)

    d = s + mw_factors(k)
    w = 2.0_qp**(2.0_qp/s)
    !w = 2.0_qp
    call init_mathe_wei(V, w)

    call create_filename(s, "f_timps_cube01", "halton", file_name)

    do i=1,size(n_steps)

      call calculate_max_pointset_size(n_steps(i), s, "M_SR", d, max_size)
      print *, "maximum pointset size = ", max_size
      allocate(x(s, max_size), startindex(s))
      startindex = 0
      call init_halton(s, startindex)
      do j=1,max_size
        call next_halton(x(:,j))
      end do
      call mathe_wei(n_steps(i), x, f_timps_cube01_np, dummy_params, rho_one_np, dummy_params, "M_SR", d, res)
      write(unit=*, fmt="(A16, F25.15)") "approximation = ", res
      write(unit=*, fmt="(A16, F25.15)") "exact value   = ", exact
      rel_err(i,k) = abs((exact-res)/exact)
      deallocate(x, startindex)

    end do

  end do
  call write_to_file(rel_err, file_name, mw_factors)

  ! EINDE TIMS FIRST TESTFUNCTION IN 10 DIMENSIONS


  ! TIMS FIRST TESTFUNCTION IN 15 DIMENSIONS
  s = 15
  exact = f_timps_cube01_exact_inf(s)

  do k = 1,size(mw_factors)

    d = s + mw_factors(k)
    w = 2.0_qp**(2.0_qp/s)
    !w = 2.0_qp
    call init_mathe_wei(V, w)

    call create_filename(s, "f_timps_cube01", "halton", file_name)

    do i=1,size(n_steps)

      call calculate_max_pointset_size(n_steps(i), s, "M_SR", d, max_size)
      print *, "maximum pointset size = ", max_size
      allocate(x(s, max_size), startindex(s))
      startindex = 0
      call init_halton(s, startindex)
      do j=1,max_size
        call next_halton(x(:,j))
      end do
      call mathe_wei(n_steps(i), x, f_timps_cube01_np, dummy_params, rho_one_np, dummy_params, "M_SR", d, res)
      write(unit=*, fmt="(A16, F25.15)") "approximation = ", res
      write(unit=*, fmt="(A16, F25.15)") "exact value   = ", exact
      rel_err(i,k) = abs((exact-res)/exact)
      deallocate(x, startindex)

    end do

  end do
  call write_to_file(rel_err, file_name, mw_factors)

  ! EINDE TIMS FIRST TESTFUNCTION IN 15 DIMENSIONS


  ! TIMS SECOND TESTFUNCTION IN 10 DIMENSIONS
  s = 10
  exact = f_timps_cube02_exact_inf(s)

  do k = 1,size(mw_factors)

    d = s + mw_factors(k)
    w = 2.0_qp**(2.0_qp/s)
    !w = 2.0_qp
    call init_mathe_wei(V, w)

    call create_filename(s, "f_timps_cube02", "halton", file_name)

    do i=1,size(n_steps)

      call calculate_max_pointset_size(n_steps(i), s, "M_SR", d, max_size)
      print *, "maximum pointset size = ", max_size
      allocate(x(s, max_size), startindex(s))
      startindex = 0
      call init_halton(s, startindex)
      do j=1,max_size
        call next_halton(x(:,j))
      end do
      call mathe_wei(n_steps(i), x, f_timps_cube02_np, dummy_params, rho_one_np, dummy_params, "M_SR", d, res)
      write(unit=*, fmt="(A16, F25.15)") "approximation = ", res
      write(unit=*, fmt="(A16, F25.15)") "exact value   = ", exact
      rel_err(i,k) = abs((exact-res)/exact)
      deallocate(x, startindex)

    end do

  end do
  call write_to_file(rel_err, file_name, mw_factors)

  ! EINDE TIMS SECOND TESTFUNCTION IN 10 DIMENSIONS


  ! TIMS SECOND TESTFUNCTION IN 15 DIMENSIONS
  s = 15
  exact = f_timps_cube02_exact_inf(s)

  do k = 1,size(mw_factors)

    d = s + mw_factors(k)
    w = 2.0_qp**(2.0_qp/s)
    !w = 2.0_qp
    call init_mathe_wei(V, w)

    call create_filename(s, "f_timps_cube02", "halton", file_name)

    do i=1,size(n_steps)

      call calculate_max_pointset_size(n_steps(i), s, "M_SR", d, max_size)
      print *, "maximum pointset size = ", max_size
      allocate(x(s, max_size), startindex(s))
      startindex = 0
      call init_halton(s, startindex)
      do j=1,max_size
        call next_halton(x(:,j))
      end do
      call mathe_wei(n_steps(i), x, f_timps_cube02_np, dummy_params, rho_one_np, dummy_params, "M_SR", d, res)
      write(unit=*, fmt="(A16, F25.15)") "approximation = ", res
      write(unit=*, fmt="(A16, F25.15)") "exact value   = ", exact
      rel_err(i,k) = abs((exact-res)/exact)
      deallocate(x, startindex)

    end do

  end do
  call write_to_file(rel_err, file_name, mw_factors)

  ! EINDE TIMS SECOND TESTFUNCTION IN 15 DIMENSIONS

  deallocate(rel_err)


contains

  subroutine create_filename(s, function_name, point_family, file_name)
    integer(kind=i4b), intent(in) :: s
    character(len=*), intent(in)  :: function_name
    character(len=*), intent(in)  :: point_family
    character(len=*), intent(out) :: file_name

    write(unit=file_name, fmt="(A,I2.2,A)") &
      "mathe_wei_dim", s, "_"//trim(function_name)//"_"//trim(point_family)//".dat"

  end subroutine create_filename


  subroutine write_to_file(rel_err, filename, mw_factors)
    real(kind=qp), dimension(:,:), intent(in)   :: rel_err
    character(len=*), intent(in)                :: filename
    integer(kind=i4b), dimension(:), intent(in) :: mw_factors

    integer(kind=i4b) :: my_unit
    integer(kind=i4b) :: j
    character(len=20) :: fmt_string


    call get_unit(my_unit)
    call checked_open(my_unit, filename, "write")

    write(unit=fmt_string, fmt="(A, I0.0, A)") "(A1, A8, ", size(mw_factors), "I25)"
    write(unit=my_unit, fmt=fmt_string) "#", "N", mw_factors

    write(unit=fmt_string, fmt="(A, I0.0, A)") "(I9, ", size(mw_factors), "ES25.15e4)"
    do j = 1,size(rel_err, dim=1)
      write(unit=my_unit, fmt=fmt_string) n_steps(j), rel_err(j,:)
    end do
    call checked_close(my_unit)

  end subroutine write_to_file


end program pub_adaptive_cube_qmc_rd2
