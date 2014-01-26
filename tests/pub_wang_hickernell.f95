program pub_wang_hickernell

  use qmcpack

  ! The number of replications for the randomizations
  integer(kind=i4b), parameter :: M = 32

  ! The parameters for the testfunction
  type(functionparams) :: params

  ! The values at which we want to see our results
  integer(kind=i4b), dimension(3), parameter :: checkvalues = (/ 1024, 4096, 16384 /)
  integer(kind=i4b), dimension(4), parameter :: dimensions = (/ 5, 10, 20, 50 /)
  !integer(kind=i4b), parameter        :: dimensions(1) = (/ 50 /)

  integer(kind=i4b) :: i, j
  integer(kind=i4b) :: my_unit


  print *, "#################### REPRODUCING ARTICLE RESULTS #################"
  print *, "Article: 'Randomized Halton Sequences', Wang and Hickernell."

  call random_seed()

  call get_unit(my_unit)
  open(unit=my_unit, file="tables.tex", status="replace", action="write")
 
  do i=1,size(dimensions)

    call write_table_header(dimensions(i))

    allocate(params%a(dimensions(i)))
    params%a = 0.01_qp
    call run_test(params, dimensions(i), M)
    params%a = 1.0_qp
    call run_test(params, dimensions(i), M)
    params%a = real((/ (j, j=1,dimensions(i)) /), kind=dp)
    call run_test(params, dimensions(i), M)
    params%a = real((/ (j*j, j=1,dimensions(i)) /), kind=dp)
    call run_test(params, dimensions(i), M)
    deallocate(params%a)

    call write_table_footer()

  end do

  close(unit=my_unit)


contains

  subroutine run_test(params, s, M)

    type(functionparams), intent(in) :: params
    integer(kind=i4b), intent(in)    :: s
    integer(kind=i4b), intent(in)    :: M

    integer(kind=i4b)                              :: i, j, c
    integer(kind=i4b), dimension(s)                :: startindex
    real(kind=qp), dimension(1)                    :: mysum
    real(kind=qp), dimension(s,1)                  :: x
    real(kind=qp), dimension(s)                    :: x1

    real(kind=qp), dimension(size(checkvalues))    :: MC_res
    real(kind=qp), dimension(M, size(checkvalues)) :: MC_res_k
    real(kind=qp), dimension(size(checkvalues))    :: MC_sigma

    real(kind=qp), dimension(size(checkvalues))    :: Shift_res
    real(kind=qp), dimension(M, size(checkvalues)) :: Shift_res_k
    real(kind=qp), dimension(size(checkvalues))    :: Shift_sigma

    real(kind=qp), dimension(size(checkvalues))    :: RanStart_res
    real(kind=qp), dimension(M, size(checkvalues)) :: RanStart_res_k
    real(kind=qp), dimension(size(checkvalues))    :: RanStart_sigma

    real(kind=qp), dimension(size(checkvalues))    :: Sqrt_res
    real(kind=qp), dimension(M, size(checkvalues)) :: Sqrt_res_k
    real(kind=qp), dimension(size(checkvalues))    :: Sqrt_sigma


    write(unit=*, fmt="(A20, I2)") "Dimension s = ", s
    write(unit=*, fmt="(A10, F10.2)") "a = ", params%a(s)

    ! Calculate value of integral for plain Monte Carlo
    do i=1,M
      mysum = 0.0_qp 
      do j=1,checkvalues(size(checkvalues))
        call random_number(x)
        mysum = mysum + f_wang_hickernell(x, params)
        do c=1,size(checkvalues)
          if (j==checkvalues(c)) then
            MC_res_k(i,c) = mysum(1)/j
          end if 
        end do
      end do
    end do
    ! Average these values
    MC_res = sum(MC_res_k, dim=1)/M
    ! Calculate the sample variance
    MC_sigma = sum((MC_res_k - spread(MC_res, 1, M))**2, dim=1)/(M*(M-1))
    write(unit=*, fmt="(A20, 3(ES10.2))") "MC_res = ", MC_res
    write(unit=*, fmt="(A20, 3(ES10.2))") "MC_sigma = ", sqrt(MC_sigma)



    ! Calculate value of integral for Random Shift Halton
    startindex = 5000 ! Wang and Hickernell skipped first 5000 points
    do i=1,M
      call init_halton(s, startindex)
      call init_random_shift(s)
      mysum = 0.0_qp 
      do j=1,checkvalues(size(checkvalues))
        call next_halton(x1)
        call random_shift(x1)
        x = reshape(x1, (/ s, 1 /))
        mysum = mysum + f_wang_hickernell(x, params)
        do c=1,size(checkvalues)
          if (j==checkvalues(c)) then
            Shift_res_k(i,c) = mysum(1)/j
          end if 
        end do
      end do
    end do
    ! Average these values
    Shift_res = sum(Shift_res_k, dim=1)/M
    ! Calculate the sample variance
    Shift_sigma = sum((Shift_res_k - spread(Shift_res, 1, M))**2, dim=1)/(M*(M-1))
    write(unit=*, fmt="(A20, 3(ES10.2))") "Shift_res = ", Shift_res
    write(unit=*, fmt="(A20, 3(ES10.2))") "Shift_sigma = ", sqrt(Shift_sigma)



    ! Calculate value of integral for Random Start Halton
    do i=1,M
      call random_startindex(startindex, huge(i4b)-16385)
      call init_halton(s, startindex)
      mysum = 0.0_qp 
      do j=1,checkvalues(size(checkvalues))
        call next_halton(x1)
        x = reshape(x1, (/ s, 1 /))
        mysum = mysum + f_wang_hickernell(x, params)
        do c=1,size(checkvalues)
          if (j==checkvalues(c)) then
            RanStart_res_k(i,c) = mysum(1)/j
          end if 
        end do
      end do
    end do
    ! Average these values
    RanStart_res = sum(RanStart_res_k, dim=1)/M
    ! Calculate the sample variance
    RanStart_sigma = sum((RanStart_res_k - spread(RanStart_res, 1, M))**2, dim=1)/(M*(M-1))
    write(unit=*, fmt="(A20, 3(ES10.2))") "RanStart_res = ", RanStart_res
    write(unit=*, fmt="(A20, 3(ES10.2))") "RanStart_sigma = ", sqrt(RanStart_sigma)



    ! Calculate value of integral for Random Start Sqrt
    do i=1,M
      call random_startindex(startindex, 100000000)
      call init_square_root_seq(s, primes(s), startindex)
      mysum = 0.0_qp 
      do j=1,checkvalues(size(checkvalues))
        call next_square_root(x1)
        x = reshape(x1, (/ s, 1 /))
        mysum = mysum + f_wang_hickernell(x, params)
        do c=1,size(checkvalues)
          if (j==checkvalues(c)) then
            Sqrt_res_k(i,c) = mysum(1)/j
          end if 
        end do
      end do
    end do
    ! Average these values
    Sqrt_res = sum(Sqrt_res_k, dim=1)/M
    write(unit=*, fmt="(A20, 3(ES7.2))") "Sqrt_res = ", Sqrt_res
    ! Calculate the sample variance
    Sqrt_sigma = sum((Sqrt_res_k - spread(Sqrt_res, 1, M))**2, dim=1)/(M*(M-1))
    write(unit=*, fmt="(A20, 3(ES10.2))") "Sqrt_sigma = ", sqrt(Sqrt_sigma)

    call write_dataline(checkvalues, MC_sigma, Shift_sigma, RanStart_sigma, Sqrt_sigma)

  end subroutine run_test


  subroutine write_table_header(s)
    integer(kind=i4b), intent(in)        :: s
    
    write(unit=my_unit, fmt="(A)") backslash//"begin{tabular}{|c|c|cccc|ccc|} "//backslash//"hline"
    write(unit=my_unit, fmt="(A)") "Dimension &  N & "//backslash//"multicolumn{4}{|c|}" // &
        "{Estimated Standard Deviation} & "//backslash//"multicolumn{3}{|c|} " // &
        " {Relative Efficiency}"//backslash//backslash
    write(unit=my_unit, fmt="(A,I2,A)") "$s=", s, "$ &       & MC & Shift & " // &
        "Ran.St. & Sqrt & Shift & Ran.St. & sqrt "//backslash//backslash//backslash//"hline"

  end subroutine write_table_header


  subroutine write_table_footer()

    write(unit=my_unit, fmt="(A)") backslash//"end{tabular}"

  end subroutine write_table_footer


  ! Write out one line of data for a certain alhpa_k
  subroutine write_dataline(N, MC_sigma, Shift_sigma, RanStart_sigma, Sqrt_sigma)
    integer(kind=i4b), dimension(:), intent(in) :: N
    real(kind=qp), dimension(:), intent(in)     :: MC_sigma
    real(kind=qp), dimension(:), intent(in)     :: Shift_sigma
    real(kind=qp), dimension(:), intent(in)     :: RanStart_sigma
    real(kind=qp), dimension(:), intent(in)     :: Sqrt_sigma

    integer(kind=i4b) :: i

    do i=1,size(N)
      write(unit=my_unit, &
        fmt="(A, I5, 4(A, ES10.2), 3(A, I7), A)") &
        " & ", N(i), " & ", &
        sqrt(MC_sigma(i)),       " & ", &
        sqrt(Shift_sigma(i)),    " & ", &
        sqrt(RanStart_sigma(i)), " & ", &
        sqrt(Sqrt_sigma(i)),     " & ", &
        int(MC_sigma(i)/Shift_sigma(i),    kind=i4b), " & ", &
        int(MC_sigma(i)/RanStart_sigma(i), kind=i4b), " & ", &
        int(MC_sigma(i)/Sqrt_sigma(i),     kind=i4b), " "//backslash//backslash
    end do
    write(unit=my_unit, fmt="(A)") backslash//"hline"

  end subroutine write_dataline

end program pub_wang_hickernell
