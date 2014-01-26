program test_utilities

  use qmcpack

  integer(kind=i4b)                            :: i, j
  integer(kind=i4b)                            :: u, v, w
  integer(kind=i4b), dimension(:), allocatable :: myperm
  integer(kind=i4b), dimension(:), allocatable :: bincoef
  integer(kind=i4b), dimension(:), allocatable :: n_vector

  call show_test_header("MOD_UTILITIES")
  call init_assert()
  

  call start_test("Testing get_nb_digits()...")
  call assert(get_nb_digits(0, 2), 1)
  call assert(get_nb_digits(1, 2), 1)
  call assert(get_nb_digits(1, 2), 1)
  call assert(get_nb_digits(7, 2), 3)
  call assert(get_nb_digits(0, 7), 1)
  call assert(get_nb_digits(56, 7), 3)
  call assert(get_nb_digits(0, 10), 1)
  call assert(get_nb_digits(1, 10), 1)
  call assert(get_nb_digits(100000, 10), 6)
  call assert(get_nb_digits(1234567, 10), 7)
  call assert(get_nb_digits(2, 11), 1)
  call assert(get_nb_digits(121, 11), 3)
  call stop_test()


  call start_test("Testing find_rightmost_zerobitpos()...")
  ! 0 = 000 -> 1
  call assert(find_rightmost_zerobitpos(0), 1)
  ! 1 = 001 -> 2
  call assert(find_rightmost_zerobitpos(1), 2)
  ! 2 = 010 -> 1
  call assert(find_rightmost_zerobitpos(2), 1)
  ! 3 = 011 -> 3
  call assert(find_rightmost_zerobitpos(3), 3)
  ! 143423431 = 1000100011000111011111000111 -> 4
  call assert(find_rightmost_zerobitpos(143423431), 4)
  ! 2^30+2^28-1 = 1001111111111111111111111111111 -> 29
  call assert(find_rightmost_zerobitpos(1342177279), 29)
  call stop_test()


  call start_test("Testing btest_real(x, k)...")
  call assert(btest_real(0.0_qp, -1), .false.)
  call assert(btest_real(0.6875_qp, -1), .true.)
  call assert(btest_real(0.6875_qp, -2), .false.)
  call assert(btest_real(0.6875_qp, -3), .true.)
  call assert(btest_real(0.6875_qp, -4), .true.)
  call assert(btest_real(1.0_qp, -1), .false.)
  call assert(btest_real(0.00390625_qp, -8), .true.)
  call stop_test()


  call start_test("Testing graycode(n)...")
  call assert(graycode(2), 3)
  call assert(graycode(3), 2)
  call assert(graycode(4), 6)
  call assert(graycode(5), 7)
  call assert(graycode(6), 5)
  call assert(graycode(7), 4)
  call assert(graycode(8), 12)
  call assert(graycode(9), 13)
  call assert(graycode(10), 15)
  call assert(graycode(11), 14)
  call assert(graycode(12), 10)
  call assert(graycode(13), 11)
  call assert(graycode(14), 9)
  call assert(graycode(15), 8)
  call assert(graycode(16), 24)
  call assert(graycode(17), 25)
  call assert(graycode(18), 27)
  call assert(graycode(19), 26)
  call assert(graycode(20), 30)
  call assert(graycode(21), 31)
  call assert(graycode(22), 29)
  call assert(graycode(23), 28)
  call assert(graycode(999), 532)
  call assert(graycode(1000), 540)
  call stop_test()


  call start_test("Testing b_ary_graycode(n)...")
  call assert(b_ary_graycode(4321, 10), 4999)
  call assert(b_ary_graycode(5432, 10), 5999)
  call assert(b_ary_graycode(3715, 10), 3444)
  call assert(b_ary_graycode(9, 11), 9)
  call assert(b_ary_graycode(10, 11), 10)
  call assert(b_ary_graycode(11, 11), 21)
  call assert(b_ary_graycode( (/ 10, 11 /), 11), (/ 10, 21 /))
  call stop_test()


  call start_test("Testing frac_part()...")
  call assert(frac_part(0.0_qp), 0.0_qp)
  call assert(frac_part(1.0_qp), 0.0_qp)
  call assert(frac_part(0.1234_qp), 0.1234_qp)
  call assert(frac_part(11.2246204830937302_qp), 0.2246204830937302_qp)
  call assert(frac_part(3.789_qp), 0.789_qp)
  call assert(frac_part((/3.789_qp, 0.1234_qp/)), (/0.789_qp, 0.1234_qp/))
  call stop_test()


  call start_test("Testing create_n_vector()...")
  allocate(n_vector(get_size_n_vector(1000, 10)))
  call create_n_vector(1000, n_vector)
!  call assert(n_vector, (/ 10, 100, 1000/))
  deallocate(n_vector)
!  call assert(create_n_vector(2000, 10), (/ 10, 100, 1000/))
!  call assert(create_n_vector(32, 2), (/ 2, 4, 8, 16, 32/))
  call stop_test()


  call start_test("Testing randperm()...")
  do i=1,20
    allocate(myperm(i))
    call randperm(i, myperm)
    print *, " Random permutation for base ", i, " = ", myperm
    deallocate(myperm)
  end do
  call stop_test()


  call start_test("Testing random_integer() by generating integers 0<= x <5...")
    do i=1,1000
      call random_integer(5, j)
      write(unit=*, fmt="(I2)", advance="no") j
      if (j >= 5) then
        print *, " ERROR: found an integer which is equal to max_int !!!"
        exit
      end if
    end do
    write(unit=*, fmt="(A)") ""
  call stop_test()


  call start_test("Testing binomial_coefficients()...")
  allocate(bincoef(1))
  call binomial_coefficients(0, bincoef)
  call assert(bincoef, (/ 1 /))
  deallocate(bincoef)
  allocate(bincoef(2))
  call binomial_coefficients(1, bincoef)
  call assert(bincoef, (/ 1, 1 /))
  deallocate(bincoef)
  allocate(bincoef(3))
  call binomial_coefficients(2, bincoef)
  call assert(bincoef, (/ 1, 2, 1 /))
  deallocate(bincoef)
  allocate(bincoef(4))
  call binomial_coefficients(3, bincoef)
  call assert(bincoef, (/ 1, 3, 3, 1/))
  deallocate(bincoef)
  allocate(bincoef(5))
  call binomial_coefficients(4, bincoef)
  call assert(bincoef, (/ 1, 4, 6, 4, 1/))
  deallocate(bincoef)
  allocate(bincoef(6))
  call binomial_coefficients(5, bincoef)
  call assert(bincoef, (/ 1, 5, 10, 10, 5, 1 /))
  deallocate(bincoef)
  allocate(bincoef(7))
  call binomial_coefficients(6, bincoef)
  call assert(bincoef, (/ 1, 6, 15, 20, 15, 6, 1/))
  deallocate(bincoef)
  call stop_test()


  call start_test("Testing nchoosek(n,k)...")
  call assert(nchoosek(0,0), 1)
  call assert(nchoosek(1,0), 1)
  call assert(nchoosek(1,1), 1)
  call assert(nchoosek(2,0), 1)
  call assert(nchoosek(2,1), 2)
  call assert(nchoosek(2,2), 1)
  call assert(nchoosek(3,0), 1)
  call assert(nchoosek(3,1), 3)
  call assert(nchoosek(3,2), 3)
  call assert(nchoosek(3,3), 1)
  call assert(nchoosek(4,0), 1)
  call assert(nchoosek(4,1), 4)
  call assert(nchoosek(4,2), 6)
  call assert(nchoosek(4,3), 4)
  call assert(nchoosek(4,4), 1)
  call assert(nchoosek(5,0), 1)
  call assert(nchoosek(5,1), 5)
  call assert(nchoosek(5,2), 10)
  call assert(nchoosek(5,3), 10)
  call assert(nchoosek(5,4), 5)
  call assert(nchoosek(5,5), 1)
  call stop_test()


  call start_test("Testing search_permutation_polynomial(p, u, v, w)...")
  call search_permutation_polynomial(3, u, v, w)
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_utilities
