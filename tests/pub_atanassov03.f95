program pub_atanassov03

  use qmcpack

  integer(kind=i4b)                            :: n, s
  real(kind=qp), dimension(:), allocatable     :: x
  type(functionparams)                         :: params
  type(integrationbounds)                      :: bounds
  real(kind=qp), dimension(5)                  :: res
  real(kind=qp), dimension(5)                  :: abs_err
  real(kind=qp), dimension(1)                  :: mysum
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b)                            :: i

  print *, "#################### REPRODUCING ARTICLE RESULTS #################"
  print *, "Article: Generating and Testing the Modified Halton Sequences"
  print *, "         Atanassov, E. and Durchove, M."

  call init_assert()

  n = 100000

  print *, "Reproducing Halton results..."
  s = 5
  print *, "Dimension = ", s
  allocate(x(s), startindex(s), params%a(s), params%u(s))
  call allocate_integration_bounds(bounds, s)
  params%a = 4.0_qp
  params%u = 2.0_qp
  call set_unit_hypercube(bounds)
  startindex = 1
  call init_halton(s, startindex, primes(s))
  mysum = 0.0_qp
  res = 0.0_qp
  abs_err = 0.0_qp
  do i=1,n
    call next_halton(x)
    mysum = mysum + f_fox86_1p(x, params)
    call update(res, abs_err, mysum, 1.0_qp, i)
  end do
  print *, "res = "
  write(unit=*, fmt="(ES8.2)") res
  print *, "abs_err = "
  write(unit=*, fmt="(ES8.2)") abs_err
  deallocate(x, startindex, params%a, params%u, bounds%lb, bounds%ub)


  call show_test_summary(get_nb_assert_errors())


contains


    subroutine update(res, abs_err, mysum, exact, n)
      real(kind=qp), dimension(:), intent(inout)       :: abs_err
      real(kind=qp), dimension(:), intent(inout)       :: res
      real(kind=qp), dimension(:), intent(in)          :: mysum
      real(kind=qp), intent(in)                        :: exact
      integer(kind=i4b), intent(in)                    :: n

      if (n == 1000) then
        res(1) = mysum(1)/n
        abs_err(1) = abs(exact-mysum(1)/n)
      else if (n == 10000) then
        res(2) = mysum(1)/n
        abs_err(2) = abs(exact-mysum(1)/n)
      else if (n == 100000) then
        res(3) = mysum(1)/n
        abs_err(3) = abs(exact-mysum(1)/n)
      else if (n == 1000000) then
        res(4) = mysum(1)/n
        abs_err(4) = abs(exact-mysum(1)/n)
      else if (n == 10000000) then
        res(5) = mysum(1)/n
        abs_err(5) = abs(exact-mysum(1)/n)
      end if
      
    end subroutine update

end program pub_atanassov03
