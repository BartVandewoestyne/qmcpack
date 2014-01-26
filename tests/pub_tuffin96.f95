program pub_tuffin96

  use qmcpack

  real(kind=qp), dimension(:,:), allocatable    :: seq_halton
  real(kind=qp), dimension(:), allocatable      :: tn2
  integer(kind=i4b), dimension(16)              :: myprimes
  integer(kind=i4b)                             :: mydim
  real(kind=qp), dimension(15)                  :: tolerances
  real(kind=qp), dimension(15)                  :: halton_exact = &
        (/ 4.34e-2_qp, 1.67e-2_qp, 6.97e-3_qp, 2.49e-3_qp, 1.91e-3_qp, &
           1.18e-3_qp, 9.86e-4_qp, 6.76e-4_qp, 4.18e-4_qp, 3.63e-4_qp, &
           2.46e-4_qp, 1.94e-4_qp, 1.71e-4_qp, 1.38e-4_qp, 1.04e-4_qp /)
  real(kind=qp), dimension(15)                  :: bw_exact = &
        (/ 4.86e-2_qp, 9.72e-3_qp, 6.15e-3_qp, 1.42e-3_qp, 6.79e-4_qp, &
           2.78e-4_qp, 1.81e-4_qp, 1.03e-4_qp, 1.76e-5_qp, 1.13e-5_qp, &
           5.34e-6_qp, 2.34e-6_qp, 1.60e-6_qp, 4.07e-7_qp, 2.11e-7_qp /)

  print *, "################# REPRODUCING ARTICLE RESULTS #################"
  print *, "Article: A new permutation choice in Halton sequences"
 
  call init_assert()

  myprimes = primes(16)


  call start_test("Checking numbers for 'Halton'...")
  do mydim = 2,16

    allocate(seq_halton(myprimes(mydim)-1, mydim), tn2(myprimes(mydim)-1))
    call halton(myprimes(mydim)-1, mydim, x=seq_halton)
    call T_N_star_squared(transpose(seq_halton), tn2)
    deallocate(seq_halton)

    tolerances = (/ 0.01e-2_qp, 0.01e-2_qp, 0.01e-3_qp, 0.01e-3_qp, &
                    0.01e-3_qp, 0.01e-3_qp, 0.01e-4_qp, 0.01e-4_qp, &
                    0.01e-4_qp, 0.01e-4_qp, 0.01e-4_qp, 0.01e-4_qp, &
                    0.01e-4_qp, 0.01e-4_qp, 0.01e-4_qp /)

    call assert(tn2(myprimes(mydim)-1), &
                                    halton_exact(mydim-1), tolerances(mydim-1))
    deallocate(tn2)

  end do
  call stop_test()


  call start_test("Checking numbers for 'B W'...")
  do mydim = 2, 16

    allocate(seq_halton(myprimes(mydim)-1, mydim), tn2(myprimes(mydim)-1))
    call halton_bw(myprimes(mydim)-1, mydim, x=seq_halton)
    call T_N_star_squared(transpose(seq_halton), tn2)
    deallocate(seq_halton)

    tolerances = (/ 0.01e-2_qp, 0.01e-3_qp, 0.01e-3_qp, 0.01e-3_qp, &
                    0.01e-4_qp, 0.01e-4_qp, 0.01e-4_qp, 0.01e-4_qp, &
                    0.01e-5_qp, 0.01e-5_qp, 0.01e-6_qp, 0.01e-6_qp, &
                    0.01e-6_qp, 0.01e-7_qp, 0.01e-7_qp /)

    call assert(tn2(myprimes(mydim)-1), &
                                    bw_exact(mydim-1), tolerances(mydim-1))
    deallocate(tn2)

  end do
  call stop_test()


  call show_test_summary(get_nb_assert_errors())

end program pub_tuffin96
