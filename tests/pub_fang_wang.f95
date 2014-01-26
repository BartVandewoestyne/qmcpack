! References:
!
!  [1] `Number-theoretic Methods in Statistics', Wang, Y. and Fang, K.-T.,
!      1994, page 27.
!
! TODO:
!   * Check why we cannot reproduce the results from table 32 in [1].  Possible
!     explanations could be...
!              - error in the routine D_N_STAR2D from the module MOD_DISCREPANCY
!              - rounding errors
!              - typo's in the book
!              - slightly other sequences were used in the book

program pub_fang_wang

  use qmcpack

  real(kind=qp), dimension(:), allocatable     :: x
  real(kind=qp), dimension(:,:), allocatable   :: x2
  real(kind=qp)                                :: myresult
  integer(kind=i4b), dimension(:), allocatable :: startindex

  print *, "######### REPRODUCING RESULTS FROM FANG AND WANG'S BOOK ##########"

  call init_assert()

  ! See page 27 of [1]
  call start_test("Testing next_gp_set_c() function...")
  call init_gp_set_c(2, 7)
  allocate(x(2))
  call next_gp_set_c(x)
  call assert(x, (/ 0.247_qp, 0.555_qp /), 0.001_qp)
  call next_gp_set_c(x)
  call assert(x, (/ 0.494_qp, 0.110_qp /), 0.001_qp)
  call next_gp_set_c(x)
  call assert(x, (/ 0.741_qp, 0.665_qp /), 0.001_qp)
  call next_gp_set_c(x)
  call assert(x, (/ 0.988_qp, 0.220_qp /), 0.001_qp)
  call next_gp_set_c(x)
  call assert(x, (/ 0.235_qp, 0.775_qp /), 0.001_qp)
  call next_gp_set_c(x)
  call assert(x, (/ 0.482_qp, 0.330_qp /), 0.001_qp)
  call next_gp_set_c(x)
  call assert(x, (/ 0.729_qp, 0.885_qp /), 0.001_qp)
  call next_gp_set_c(x)
  call assert(x, (/ 0.976_qp, 0.440_qp /), 0.001_qp)
  deallocate(x)
  call stop_test()


  ! See page 27 of [1]
  call start_test("Testing gp_set_c()...")
  allocate(x2(8, 2))
  call gp_set_c(n=8, s=2, p=7, x=x2)
  call assert(x2(1,:), (/ 0.247_qp, 0.555_qp /), 0.001_qp)
  call assert(x2(2,:), (/ 0.494_qp, 0.110_qp /), 0.001_qp)
  call assert(x2(3,:), (/ 0.741_qp, 0.665_qp /), 0.001_qp)
  call assert(x2(4,:), (/ 0.988_qp, 0.220_qp /), 0.001_qp)
  call assert(x2(5,:), (/ 0.235_qp, 0.775_qp /), 0.001_qp)
  call assert(x2(6,:), (/ 0.482_qp, 0.330_qp /), 0.001_qp)
  call assert(x2(7,:), (/ 0.729_qp, 0.885_qp /), 0.001_qp)
  call assert(x2(8,:), (/ 0.976_qp, 0.440_qp /), 0.001_qp)
  deallocate(x2)
  call stop_test()


  ! NOT REPRODUCIBLE!
  call start_test("Testing d_n_star2d() with square root points...")
  allocate(x2(2584, 2), startindex(2))
  startindex = 1
  call square_root_seq(2584, 2, startindex, x=x2)
  call d_n_star2d(transpose(x2(1:34,:)), myresult)
  call assert(myresult, 0.0770_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:89,:)), myresult)
  call assert(myresult, 0.0469_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:144,:)), myresult)
  call assert(myresult, 0.0249_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:233,:)), myresult)
  call assert(myresult, 0.0239_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:610,:)), myresult)
  call assert(myresult, 0.0152_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:1597,:)), myresult)
  call assert(myresult, 0.0048_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:2584,:)), myresult)
  call assert(myresult, 0.0028_qp, 0.0001_qp)
  deallocate(x2, startindex)
  call stop_test()


  ! NOT REPRODUCIBLE!
  call start_test("Testing d_n_star2d() with gp_set_b()...")
  allocate(x2(2584, 2), startindex(2))
  startindex = 1
  call gp_set_b(n=2584, s=2, p=2, n_start=startindex, x=x2)
  call d_n_star2d(transpose(x2(1:34,:)), myresult)
  call assert(myresult, 0.0693_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:89,:)), myresult)
  call assert(myresult, 0.0484_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:144,:)), myresult)
  call assert(myresult, 0.0253_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:233,:)), myresult)
  call assert(myresult, 0.0180_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:610,:)), myresult)
  call assert(myresult, 0.0086_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:1597,:)), myresult)
  call assert(myresult, 0.0050_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:2584,:)), myresult)
  call assert(myresult, 0.0037_qp, 0.0001_qp)
  deallocate(x2, startindex)
  call stop_test()


  ! NOT REPRODUCIBLE!
  call start_test("Testing d_n_star2d() with gp_set_c()...")
  allocate(x2(2584, 2), startindex(2))
  startindex = 1
  call gp_set_c(n=2584, s=2, p=7, n_start=startindex, x=x2)
  call d_n_star2d(transpose(x2(1:34,:)), myresult)
  call assert(myresult, 0.1853_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:89,:)), myresult)
  call assert(myresult, 0.0786_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:144,:)), myresult)
  call assert(myresult, 0.0481_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:233,:)), myresult)
  call assert(myresult, 0.0350_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:610,:)), myresult)
  call assert(myresult, 0.0205_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:1597,:)), myresult)
  call assert(myresult, 0.0072_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:2584,:)), myresult)
  call assert(myresult, 0.0053_qp, 0.0001_qp)
  deallocate(x2, startindex)
  call stop_test()


  write(unit=*, fmt="(A)") "Testing d_n_star2d() with Hammersley point set..."
  write(unit=*, fmt="(A)") "  Note: the values published in the book are obtained"
  write(unit=*, fmt="(A)") "        using another definition of the Hammersley point"
  write(unit=*, fmt="(A)") "        set, namely the one from Shaw88!!!" 
  allocate(x2(2584, 2), startindex(1))
  startindex = 1
  call hammersley(34, 2, primes(1), startindex, def="Shaw", x=x2(1:34,:))
  call d_n_star2d(transpose(x2(1:34,:)), myresult)
  call assert(myresult, 0.0864_qp, 0.0001_qp)

  call hammersley(89, 2, primes(1), startindex, def="Shaw", x=x2(1:89,:))
  call d_n_star2d(transpose(x2(1:89,:)), myresult)
  call assert(myresult, 0.0406_qp, 0.0001_qp)

  call hammersley(144, 2, primes(1), startindex, def="Shaw", x=x2(1:144,:))
  call d_n_star2d(transpose(x2(1:144,:)), myresult)
  call assert(myresult, 0.0251_qp, 0.0001_qp)

  call hammersley(233, 2, primes(1), startindex, def="Shaw", x=x2(1:233,:))
  call d_n_star2d(transpose(x2(1:233,:)), myresult)
  call assert(myresult, 0.0169_qp, 0.0001_qp)

  call hammersley(610, 2, primes(1), startindex, def="Shaw", x=x2(1:610,:))
  call d_n_star2d(transpose(x2(1:610,:)), myresult)
  call assert(myresult, 0.0073_qp, 0.0001_qp)

  call hammersley(1597, 2, primes(1), startindex, def="Shaw", x=x2(1:1597,:))
  call d_n_star2d(transpose(x2(1:1597,:)), myresult)
  call assert(myresult, 0.0031_qp, 0.0001_qp)

  call hammersley(2584, 2, primes(1), startindex, def="Shaw", x=x2(1:2584,:))
  call d_n_star2d(transpose(x2(1:2584,:)), myresult)
  call assert(myresult, 0.0020_qp, 0.0001_qp)
  deallocate(x2, startindex)
  call stop_test()


  ! NOT REPRODUCIBLE!
  call start_test("Testing d_n_star2d() with Haber points...")
  allocate(x2(2584, 2), startindex(2))
  startindex = 1
  call haber(2584, 2, startindex, x=x2)
  call d_n_star2d(transpose(x2(1:34,:)), myresult)
  call assert(myresult, 0.1644_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:89,:)), myresult)
  call assert(myresult, 0.0947_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:144,:)), myresult)
  call assert(myresult, 0.0838_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:233,:)), myresult)
  call assert(myresult, 0.0604_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:610,:)), myresult)
  call assert(myresult, 0.0327_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:1597,:)), myresult)
  call assert(myresult, 0.0222_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:2584,:)), myresult)
  call assert(myresult, 0.0217_qp, 0.0001_qp)
  deallocate(x2, startindex)
  call stop_test()
  

  call start_test("Testing d_n_star2d() with Halton points...")
  allocate(x2(2584, 2), startindex(2))
  startindex = 1
  call halton(2584, 2, startindex, x=x2)
  call d_n_star2d(transpose(x2(1:34,:)), myresult)
  call assert(myresult, 0.1106_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:89,:)), myresult)
  call assert(myresult, 0.0432_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:144,:)), myresult)
  call assert(myresult, 0.0328_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:233,:)), myresult)
  call assert(myresult, 0.0207_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:610,:)), myresult)
  call assert(myresult, 0.0109_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:1597,:)), myresult) ! @home, this is where ifort 9.0 crashes...
  call assert(myresult, 0.0052_qp, 0.0001_qp)
  call d_n_star2d(transpose(x2(1:2584,:)), myresult)
  call assert(myresult, 0.0030_qp, 0.0001_qp)
  deallocate(x2, startindex)
  call stop_test()
  
  
  call show_test_summary(get_nb_assert_errors())

end program pub_fang_wang
