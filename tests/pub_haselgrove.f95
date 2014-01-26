! References:
!
!   [1] "A Method for Numerical Integration", C. B. Haselgrove, Mathematics
!       of Computation, Vol. 15, 1961, pages 323-337.

program pub_haselgrove

  use qmcpack

  integer(kind=i4b)    :: s
  real(kind=qp)        :: myresult
  type(functionparams) :: params

  call show_test_header("MOD_HASELGROVE")
  call init_assert()

  s = 5
 
  allocate(params%a(s), params%u(s))
  params%a = 1
  params%u = 0

  call init_haselgrove(s)

  ! TODO: check if there is a problem with very large N
  call haselgrove(f_haselgrove61_1p, params, 120000, 2, myresult)
  print *, "myresult = ", myresult
  
!  call haselgrove(f_haselgrove61_1p, params, 20000, 1, myresult)
!  print *, "myresult = ", myresult

  call show_test_summary(get_nb_assert_errors())

end program pub_haselgrove
