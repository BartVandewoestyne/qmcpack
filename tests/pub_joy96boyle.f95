! References:
!
!   [1] 'Quasi-Monte Carlo Methods in Numerical Finance', Corwin Joy,
!       Phelim P. Boyle, Ken Seng Tan, Management Science, Vol. 42, No. 6,
!       June 1996.
!
program pub_joy96boyle

  use qmcpack

  integer(kind=i4b)                            :: n, s
  real(kind=qp), dimension(:,:), allocatable   :: x
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b), dimension(:), allocatable :: step
  integer(kind=i4b)                            :: i
  real(kind=qp), dimension(:), allocatable     :: x_1p
  integer                                      :: my_unit

  print *, "#################### REPRODUCING ARTICLE RESULTS #################"
  print *, "Article: Joy, Boyle, Tan, 1996"
  call init_assert()


  call start_test("Reproducing results from Table 1 on page 930")
  n = 11
  s = 3
  allocate(x(n,s), step(s), startindex(s))
  step = 1
  startindex = 1
  call init_faure(n, s, "None", init_startindex=startindex, init_step=step)
  do i=1,n
    call next_faure(x(i,:))
  end do
  call assert(x(1,:),  (/  1.0_qp/3,   1.0_qp/3,   1.0_qp/3 /))
  call assert(x(2,:),  (/  2.0_qp/3,   2.0_qp/3,   2.0_qp/3 /))
  call assert(x(3,:),  (/  1.0_qp/9,   4.0_qp/9,   7.0_qp/9 /))
  call assert(x(4,:),  (/  4.0_qp/9,   7.0_qp/9,   1.0_qp/9 /))
  call assert(x(5,:),  (/  7.0_qp/9,   1.0_qp/9,   4.0_qp/9 /))
  call assert(x(6,:),  (/  2.0_qp/9,   8.0_qp/9,   5.0_qp/9 /))
  call assert(x(7,:),  (/  5.0_qp/9,   2.0_qp/9,   8.0_qp/9 /))
  call assert(x(8,:),  (/  8.0_qp/9,   5.0_qp/9,   2.0_qp/9 /))
  call assert(x(9,:),  (/  1.0_qp/27, 16.0_qp/27, 13.0_qp/27 /))
  call assert(x(10,:), (/ 10.0_qp/27, 25.0_qp/27, 22.0_qp/27 /))
  call assert(x(11,:), (/ 19.0_qp/27,  7.0_qp/27,  4.0_qp/27 /))
  deallocate(x, step , startindex)
  call free_faure()
  call stop_test()


  ! Note: the startindex is not mentioned, but it is probably 2**4 or 2**4-1.
  ! (for 2**4-1 we don't have the point in the upper right corner but we do
  ! have the point in the lower right corner.  For 2**4 we don't have the
  ! point in the lower right corner but we do have the point in the upper right
  ! corner.)
  !
  call start_test("Reproducing points from Figure 1 and 2 on page 931")
  n = 128
  s = 2
  allocate(x_1p(s), step(s), startindex(s))
  step = 1
  startindex = 2**4-1
  call get_unit(my_unit)
  call checked_open(my_unit, "pub_joy96boyle.dat", "write")
  call init_faure(n, s, "None", init_startindex=startindex, init_step=step)
  do i=1,n
    call next_faure(x_1p)
    write(unit=my_unit, fmt="(2F18.15)") x_1p
  end do
  deallocate(x_1p, step , startindex)
  call free_faure()
  call stop_test()

end program pub_joy96boyle
