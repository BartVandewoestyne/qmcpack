! References:
!
!   [1] 'Good permutations for extreme discrepancy', Faure, Henri, Journal of
!       Number Theory, vol. 42, pages 47-56, 1992.

program test_permutation

  use qmcpack

  integer(kind=i4b)                            :: i
  integer, parameter                           :: nbdims=20
  type(permutation_pointer), dimension(nbdims) :: x
  character(len=10)                            :: formatstring

  call show_test_header("MOD_PERMUTATION")
  call init_assert()

  call start_test("Testing init_permutations()...")
  x = init_permutations(allprimes(1:nbdims))
  do i=1,nbdims
    write(unit=formatstring, fmt="(A,I0,A,I0,A)") "(", size(x(i)%p), "I", get_nb_digits(maxval(x(nbdims)%p), 10)+1, ")"
    write(unit=*, fmt=formatstring) x(i)%p
    deallocate(x(i)%p)
  end do
  call stop_test()


  call start_test("Testing permutations by Braaten and Weller:")
  x(1:16) = init_perm_braaten_weller(16)
  do i=1,16
    write(unit=formatstring, fmt="(A,I0,A,I0,A)") "(", size(x(i)%p), "I", get_nb_digits(maxval(x(16)%p), 10)+1, ")"
    write(unit=*, fmt=formatstring) x(i)%p
    deallocate(x(i)%p)
  end do
  call stop_test()


  call start_test("Testing permutation MCL by Tuffin:")
  x(1:16) = init_perm_mcl(16)
  do i=1,16
    write(unit=formatstring, fmt="(A,I0,A,I0,A)") "(", size(x(i)%p), "I", get_nb_digits(maxval(x(16)%p), 10)+1, ")"
    write(unit=*, fmt=formatstring) x(i)%p
    deallocate(x(i)%p)
  end do
  call stop_test()


  call start_test("Testing reverse permutations:")
  x(1:nbdims) = init_perm_reverse(20)
  do i=1,20
    write(unit=formatstring, fmt="(A,I0,A,I0,A)") "(", size(x(i)%p), "I", get_nb_digits(maxval(x(16)%p), 10)+1, ")"
    write(unit=*, fmt=formatstring) x(i)%p
    deallocate(x(i)%p)
  end do
  call stop_test()


  ! See [1]
  call start_test("Testing permutations by Faure...")
  write(unit=*, fmt="(A)", advance="no") "  Testing generation of all bases... "
  x(1:7) = init_perm_faure(7, .false.)
  call assert(x(1)%p, (/ 0, 1 /))
  call assert(x(2)%p, (/ 0, 1, 2 /))
  call assert(x(3)%p, (/ 0, 2, 1, 3 /))
  call assert(x(4)%p, (/ 0, 3, 2, 1, 4 /))
  call assert(x(5)%p, (/ 0, 2, 4, 1, 3, 5 /))
  call assert(x(6)%p, (/ 0, 2, 5, 3, 1, 4, 6 /))
  call assert(x(7)%p, (/ 0, 4, 2, 6, 1, 5, 3, 7 /))
  do i=1,7
    deallocate(x(i)%p)
  end do
  call stop_test()


  call start_test("  Testing generation of only prime bases... ")
  x(1:4) = init_perm_faure(4, .true.)
  call assert(x(1)%p, (/ 0, 1 /))
  call assert(x(2)%p, (/ 0, 1, 2 /))
  call assert(x(3)%p, (/ 0, 3, 2, 1, 4 /))
  call assert(x(4)%p, (/ 0, 2, 5, 3, 1, 4, 6 /))
  do i=1,4
    deallocate(x(i)%p)
  end do
  call stop_test()

  call show_test_summary(get_nb_assert_errors())

end program test_permutation
