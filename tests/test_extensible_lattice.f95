program test_extensible_lattice

  use qmcpack

  integer(kind=i4b)                            :: N, s
  integer(kind=i4b), dimension(:), allocatable :: z
  integer(kind=i4b)                            :: i
  real(kind=qp), dimension(:), allocatable     :: x
  real(kind=qp), dimension(:), allocatable     :: x_1p
  integer(kind=i4b), dimension(:), allocatable :: startindex, step, b
  character(len=10)                            :: fmt_string
  integer(kind=i4b)                            :: my_unit
  logical                                      :: equal
  integer                                      :: ios

  call show_test_header("MOD_EXTENSIBLE_LATTICE")
  call init_assert()

  ! Write out an extensible lattice to a file.
  !
  call start_test("Writing out extensible lattice with 128 points to extensible_lattice_128p_3d.dat...")
  N = 128
  s = 3
  allocate(x(s), z(s))
  call get_vector_nuyens_kuo(s, z)
  call get_unit(unit_number)
  open(unit=unit_number, file="extensible_lattice_128p_3d.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F20.15)"
  call init_extensible_lattice(s, z)
  do i=1,N
    call next_extensible_lattice(x)
    write(unit=unit_number, fmt=fmt_string) x(:)
  end do
  call free_extensible_lattice()
  deallocate(x, z)
  close(unit=unit_number)
  call stop_test()


  ! Write out an extensible lattice with graycode ordening to a file.
  !
  call start_test("Writing out extensible lattice with graycode ordening with 128 points to extensible_lattice_128p_3d.dat...")
  N = 128
  s = 3
  allocate(x(s), z(s))
  call get_vector_nuyens_kuo(s, z)
  call get_unit(unit_number)
  open(unit=unit_number, file="extensible_lattice_gray_128p_3d.dat", iostat=ios, &
       status="replace", access="sequential", action="write")
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F20.15)"
  call init_extensible_lattice(s, z)
  do i=1,N
    call next_extensible_lattice_gray(x)
    write(unit=unit_number, fmt=fmt_string) x(:)
  end do
  call free_extensible_lattice()
  deallocate(x, z)
  close(unit=unit_number)
  call stop_test()


  call start_test("Testing next_extensible_lattice_gray() by writing points to" // &
    " extensible_lattice_gray_30p_10d.dat and comparing with Dirk's implementation...")
  n = 30
  s = 10
  write(unit=fmt_string, fmt="(A, I0.0, A)") "(", s, "F9.6)"
  allocate(x_1p(s), startindex(s), step(s), z(s), b(s))
  startindex = 0
  step = 1
  b = 2
  call get_vector_nuyens_kuo(s, z)
  call init_extensible_lattice(s, z, b, startindex, step)
  call get_unit(my_unit)
  open(unit=my_unit, file="extensible_lattice_gray_30p_10d.dat", iostat=ios,        &
       status="replace", access="sequential", action="write")
  if (ios == 0) then
    do i=1,n
      call next_extensible_lattice_gray(x_1p)
      write(unit=my_unit, fmt=fmt_string) x_1p
    end do
    close(unit=my_unit)
  else
    print *, "ERROR: could not write extensible lattice (graycode version) dataset to file!"
    call increase_nb_assert_errors(1)
  end if
  call compare_pointsets("extensible_lattice_gray_30p_10d.dat", &
                         "test_data/extensible_lattice_gray_dirk_30p_10d.dat", &
                         0.00000000000001_qp, equal)
  if (.not. equal) then
    call increase_nb_assert_errors(1)
  end if
  call free_extensible_lattice()
  deallocate(x_1p, startindex, step, z, b)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_extensible_lattice
