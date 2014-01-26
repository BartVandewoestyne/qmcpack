! References:
!
!   [1] Shaw, J. E., `A quasirandom approach to integration in Bayesian
!       statistics', Annals of Statistics, vol. 16, 1988, pages 895-914.
!
! Note: watch out when testing against sequences for which the pointset
!       depends on the number of points!  For each test against a
!       discrepancy value, the pointset should be re-generated with the
!       specified number of points!!!
!
! TODO:
!   * Add tests for Irrational sequences and Sobol' sequence.

program pub_shaw88

  use qmcpack

  real(kind=qp), dimension(:,:), allocatable   :: x2
  real(kind=qp)                                :: myresult
  integer(kind=i4b), dimension(:), allocatable :: startindex

  write(unit=*, fmt="(A)") "######### REPRODUCING RESULTS FROM shaw88 ##########"
  call init_assert()

  ! See page 904 of [1]
  call start_test("Testing d_n_star2d() with Haber points...")
  allocate(x2(256, 2), startindex(2))
  startindex = 1
  call haber(256, 2, startindex, x=x2)
  call d_n_star2d(transpose(x2(1:32,:)), myresult)
  call assert(myresult, 0.175382_qp, 0.000001_qp)
  call d_n_star2d(transpose(x2(1:64,:)), myresult)
  call assert(myresult, 0.109077_qp, 0.000001_qp)
  call d_n_star2d(transpose(x2(1:128,:)), myresult)
  call assert(myresult, 0.116184_qp, 0.000001_qp)
  call d_n_star2d(transpose(x2(1:256,:)), myresult)
  call assert(myresult, 0.064437_qp, 0.000001_qp)
  deallocate(x2, startindex)
  call stop_test()


  ! See page 904 of [1]
  call start_test("Testing d_n_star2d() with Halton points...")
  allocate(x2(256, 2), startindex(2))
  startindex = 0
  call halton(256, 2, startindex, x=x2)
  call d_n_star2d(transpose(x2(1:32,:)), myresult)
  call assert(myresult, 0.104167_qp, 0.000001_qp)
  call d_n_star2d(transpose(x2(1:64,:)), myresult)
  call assert(myresult, 0.052083_qp, 0.000001_qp)
  call d_n_star2d(transpose(x2(1:128,:)), myresult)
  call assert(myresult, 0.036651_qp, 0.000001_qp)
  call d_n_star2d(transpose(x2(1:256,:)), myresult)
  call assert(myresult, 0.018760_qp, 0.000001_qp)
  deallocate(x2, startindex)
  call stop_test()
  

  ! See page 904 of [1]
  call start_test("Testing d_n_star2d() with Hammersley point set...")
  allocate(x2(256, 2), startindex(1))
  startindex = 0

  call hammersley(32, 2, primes(1), startindex, def="Shaw", x=x2(1:32,:))
  call d_n_star2d(transpose(x2(1:32,:)), myresult)
  call assert(myresult, 0.097656_qp, 0.000001_qp)

  call hammersley(64, 2, primes(1), startindex, def="Shaw", x=x2(1:64,:))
  call d_n_star2d(transpose(x2(1:64,:)), myresult)
  call assert(myresult, 0.053711_qp, 0.000001_qp)

  call hammersley(128, 2, primes(1), startindex, def="Shaw", x=x2(1:128,:))
  call d_n_star2d(transpose(x2(1:128,:)), myresult)
  call assert(myresult, 0.029541_qp, 0.000001_qp)

  call hammersley(256, 2, primes(1), startindex, def="Shaw", x=x2(1:256,:))
  call d_n_star2d(transpose(x2(1:256,:)), myresult)
  call assert(myresult, 0.016052_qp, 0.000001_qp)
  deallocate(x2, startindex)
  call stop_test()


  call show_test_summary(get_nb_assert_errors())

end program pub_shaw88
