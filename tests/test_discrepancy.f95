program test_discrepancy

  use qmcpack

  real(kind=qp), dimension(:), allocatable           :: x
  real(kind=qp), dimension(:,:), allocatable         :: x_np
  real(kind=qp), dimension(:), allocatable           :: tn2, d2
  real(kind=qp)                                      :: myresult
  integer(kind=i4b)                                  :: i
  integer(kind=i4b)                                  :: s, N

  call show_test_header("MOD_DISCREPANCY")
  call init_assert()


  call start_test("Performing simple test for T_N_star_squared with 2x3 pointset...")
  s = 2
  N = 3
  allocate(x_np(s,N), tn2(N))
  x_np(1,:) = (/1, 3, 5/)
  x_np(2,:) = (/2, 4, 6/)
  call T_N_star_squared(x_np, tn2)
  call assert(tn2(1),     1.0_qp/9)
  call assert(tn2(2),  -457.0_qp/18)
  call assert(tn2(3), -1321.0_qp/9)
  deallocate(x_np, tn2)
  call stop_test()


  call start_test("Performing simple test for T_N_star_squared with 1x2 pointset...")
  s = 1
  N = 2
  allocate(x_np(s,N), tn2(N))
  x_np(1,:) = (/0.6_qp, 0.2_qp/)
  call T_N_star_squared(x_np, tn2)
  call assert(tn2(N), 1.0_qp/30)
  deallocate(x_np, tn2)
  call stop_test()


  call start_test("Performing simple test for D2_hickernell...")
  s = 3
  N = 5
  allocate(x_np(s,N), tn2(N), d2(N))
  call random_number(x_np)
  call D2_hickernell_squared(x_np, d2)
  call D2_hickernell_straightforward(x_np, tn2)
  call assert(d2, tn2)
  deallocate(x_np, tn2, d2)
  call stop_test()


  call start_test("Performing simple test for discd2 with 2x3 pointset...")
  s = 2
  N = 3
  allocate(x_np(N,s), tn2(N))
  x_np(1,:) = (/1, 2/)
  x_np(2,:) = (/3, 4/)
  x_np(3,:) = (/5, 6/)
  call discd2(N, s, x_np, tn2)
  call assert(tn2(1),     1.0_qp/9)
  call assert(tn2(2),  -457.0_qp/18)
  call assert(tn2(3), -1321.0_qp/9)
  deallocate(x_np, tn2)
  call stop_test()


  call start_test("Performing simple test for discd2 with 1x2 pointset...")
  s = 1
  N = 2
  allocate(x_np(N,s), tn2(N))
  x_np(1,:) = (/0.6_qp/)
  x_np(2,:) = (/0.2_qp/)
  call discd2(N, s, x_np, tn2)
  call assert(tn2(N), 1.0_qp/30.0_qp)
  deallocate(x_np, tn2)
  call stop_test()


  call start_test("Performing simple test for T_N_extreme_squared with 2x2 pointset...")
  s = 2
  N = 2
  allocate(x_np(s,N), tn2(N))
  x_np(1,:) = (/1, 3/)
  x_np(2,:) = (/2, 4/)
  call T_N_extreme_squared(x_np, tn2)
  call assert(tn2(N), 24-18+1.0_qp/144)
  deallocate(x_np, tn2)
  call stop_test()


  call start_test("Testing d_n_star1d()...")
  allocate(x(9))
  x = (/ (0.1_qp*i, i=1,9) /)
  call d_n_star1d(x, myresult)
  call assert(myresult, 0.1_qp)
  deallocate(x)
  ! The discrepancy of (2*i-1)/(2*N) = 1/(2*N)
  allocate(x(7))
  x = (/ ((2.0_qp*i-1)/(2*7), i=1,7) /)
  call d_n_star1d(x, myresult)
  call assert(myresult, 0.5_qp/7)
  deallocate(x)
  call stop_test()


  call start_test("Testing d_n_extreme1d()...")
  allocate(x(9))
  x = (/ (0.1_qp*i, i=1,9) /)
  call d_n_extreme1d(x, myresult)
  call assert(myresult, 0.2_qp)
  deallocate(x)
  call stop_test()

  
  ! Calculate discrepancy for random sequence
  !s = 16 
  !N = 100000

  !allocate(x(N,s))
  !allocate(tn2(N))

  !call random_number(x)
  !call halton(N,s,1, allprimes(1:s), x)
  !call T_N_star_squared(transpose(x), tn2)
  !open(unit_number, file='halton_tn_star_squared.dat', IOSTAT=IERR, ERR=10, &
  !      status='replace', access='sequential', action='write')
  !write(unit_number, FMT='(F32.30)') sqrt(tn2)
  !close(unit_number)


  !call halton_bw(N,s,1, x)
  !call T_N_star_squared(transpose(x), tn2)
  !open(unit_number, file='halton_bw_tn_star_squared.dat', IOSTAT=IERR, ERR=10, &
  !      status='replace', access='sequential', action='write')
  !write(unit_number, FMT='(F32.30)') sqrt(tn2)
  !close(unit_number)

  !deallocate(x)
  !deallocate(tn2)


!10 if (IERR .GT. 0) then
!     print *, "Error opening file"
!   end if


  call show_test_summary(get_nb_assert_errors())


contains

  ! This is a straightforward, but slower version of D2_hickernell to
  ! compare the results of D2_hickernell with.
  !
  subroutine D2_hickernell_straightforward(x, D2)

    real(kind=qp), dimension(:,:), intent(in) :: x
    real(kind=qp), dimension(:), intent(out)  :: D2

    integer(kind=i4b) :: s, nbpoints, N, i, k
    real(kind=qp)     :: temp

    s = size(x, 1)
    nbpoints = size(x, 2)

    do N=1,nbpoints

      D2(N) = (4.0_qp/3)**s

      temp = 0
      do i=1,N
        temp = temp + product(1.5_qp-0.5*x(:,i)**2)
      end do
      D2(N) = D2(N) - 2.0_qp/N*temp

      temp = 0
      do i=1,N
        do k=1,N
          temp = temp + product(2-max(x(:,i),x(:,k)))
        end do
      end do
      D2(N) = D2(N) + temp/N**2

    end do

  end subroutine D2_hickernell_straightforward

end program test_discrepancy
