program pub_fox86

  use qmcpack

  integer(kind=i4b)                            :: n, s
  integer(kind=i4b)                            :: i, j, k
  integer(kind=i4b), dimension(:), allocatable :: startindex
  integer(kind=i4b), dimension(:), allocatable :: step
  real(kind=qp), dimension(:), allocatable     :: x
  real(kind=qp)                                :: myres
  real(kind=qp)                                :: mysum
  type(functionparams)                         :: params
  integer(kind=i4b), dimension(6)              :: Nvec = (/ 500, 1000, 7000, 20000, 40000, 100000 /)
  integer(kind=i4b), dimension(6)              :: svec = (/ 4, 7, 13, 20, 25, 40 /)

  do i = 1,size(Nvec)

    n = Nvec(i)

    do j = 1,size(svec)

      s = svec(j)

      allocate(x(s), startindex(s), step(s), params%a(s))
      step = 1
      startindex = prime_ge(s)**4-1
      params%a = 0
      mysum = 0.0_qp

      call init_faure(n, s, init_scrambletype="None", init_startindex=startindex, init_step=step)
      do k = 1,n

        call next_faure(x)

        mysum = mysum + f_wang03fang_1p(x, params)

        ! Compute value for integral
        myres = mysum/k

      end do
      deallocate(x, startindex, step, params%a)
      call free_faure()

      write(unit=*, fmt="(F6.3, A)", advance="no") myres, "   "

    end do

    write(unit=*, fmt="(A)") ""

  end do


end program pub_fox86
