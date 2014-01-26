! Module to test the routines for the Haber sequence.
!
! Note:
!   * It seems like we lose precision quite fast... this is
!     probably due to the fact that we keep only the fractional part
!     of some large numbers.

program test_haber

  use qmcpack

  integer(kind=i4b)                            :: n, s
  real(kind=qp), dimension(:), allocatable     :: x_1p
  real(kind=qp), dimension(:,:), allocatable   :: x_np
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b)                            :: i
  integer(kind=i4b), dimension(:), allocatable :: myprimes

  call show_test_header("MOD_HABER")
  call init_assert()

  call start_test("Testing block assignment with haber()...")
  n = 1000
  s = 3
  allocate(x_np(n,s), startindex(s), myprimes(s))
  startindex = 1
  myprimes = primes(s)
  call haber(n, s, startindex, myprimes, x_np)
  call assert(x_np(500,:),                 &
              (/ 0.2486872301669791_qp,    &
                 0.3636480018612929_qp,    &
                 0.5141818486736156_qp /)) 
  call assert(x_np(1000,:),                &
              (/ 0.8879677340719252_qp,    & 
                 0.4291882230854104_qp,    &
                 0.0227386447430527_qp /)) ! Small precision!
  deallocate(x_np, startindex, myprimes)
  call stop_test()


  call start_test("Testing next_haber() against values calculated with Maple...")
  allocate(x_1p(3), startindex(3))
  startindex = 1
  call init_haber(3, startindex)
  do i=1,1000
    call next_haber(x_1p)
  end do
  call assert(x_1p, (/ 0.8879677340719252_qp,    &
                       0.4291882230854104_qp,    &
                       0.0227386447430527_qp /)) ! Small precision!
  call free_haber()
  deallocate(x_1p, startindex)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_haber
