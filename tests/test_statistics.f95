program test_statistics

  use qmcpack

  real(kind=qp)     :: x
  integer(kind=i4b) :: i
  type(statistic)   :: homer_stat, marge_stat

  call show_test_header("MOD_STATISTICS")
  call init_assert()

  call create_statistic("Homer")
  do i=1,10
    x = i
    call update_statistic("Homer", x)
  end do

  call create_statistic("Marge")
  do i=1,10
    x = 2*i
    call update_statistic("Marge", x)
  end do

  call get_statistic("Homer", homer_stat)
  call get_statistic("Marge", marge_stat)

  call assert(homer_stat%average, 5.5_qp)    
  call assert(homer_stat%nb_values, 10)    

  call assert(marge_stat%average, 11.0_qp)    
  call assert(marge_stat%nb_values, 10)    
  
  call free_statistics()

  call show_test_summary(get_nb_assert_errors())

end program test_statistics
