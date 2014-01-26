! Module implementing several kinds of Probability Density Functions.
!
! References:
!
!  [1] Press, William H. and Teukolsky, Saul A. and Vetterlin, William T. and
!      Flannery, Brian P., `Numerical Recipes in Fortran 77', second edition,
!      Cambridge Univeristy Press, 1992, ISBN 0-521-43064-X.

module mod_pdf

  use numeric_kinds
  use mod_constants

!  implicit none
  private

  public :: randn_box_muller
  public :: norm_pdf


contains


  ! Successive calls to this function return independent, normally distributed
  ! pseudo-random numbers with zero mean and unit standard deviation.
  !
  ! Note:
  !   * To generate the random numbers from [0,1) that are necessary for the
  !     Box-Muller method, we simply use the RANDOM_NUMBER intrinsic here.
  !
  ! See [1], page 280 for details on the algorithm.
  !
  subroutine randn_box_muller(x, mu, sigma)

    real(kind=qp), intent(out) :: x
    real(kind=qp), intent(in)  :: mu, sigma

    logical, save       :: have_second_sample = .false.
    real(kind=qp), save :: saved_sample
    real(kind=qp)       :: myrand1, myrand2, v1, v2, rsq, fac

    if (.not. have_second_sample) then

      do

        ! We don't have an extra deviate handy, so pick two uniform numbers in
        ! the square extending from -1 to +1 in each direction.
        call random_number(myrand1)
        call random_number(myrand2)
        v1 = 2.0_qp*myrand1 - 1.0_qp
        v2 = 2.0_qp*myrand2 - 1.0_qp

        ! See if they are in the unit circle, and if they are not, try again.
        rsq = v1*v1 + v2*v2
        if (rsq < 1.0_qp .and. rsq /= 0.0_qp) then
          exit
        end if

      end do

      ! Now make the Box-Muller transformation to get two normal deviates.
      ! Return one and save the other for next time.
      fac = sqrt(-2.0_qp*log(rsq)/rsq)
      saved_sample = v1*fac*sigma + mu
      x = v2*fac*sigma + mu

      ! Set the flag indicating that we still have an extra deviate handy.
      have_second_sample = .true.

    else ! We have an extra deviate handy,...

      ! ... so return it,...
      x = saved_sample

      ! ... and unset the flag.
      have_second_sample = .false.

    end if
    
  end subroutine randn_box_muller


  ! Computes the normal probability density function with the given mean
  ! and variance for each component of x.
  !
  elemental function norm_pdf(x, mean, variance) result (y)
    real(kind=qp), intent(in) :: x
    real(kind=qp), intent(in) :: mean, variance
    real(kind=qp)             :: y

    y = 1/(sqrt(variance*2*pi))*exp( -(x-mean)**2/(2*variance) )

  end function norm_pdf

end module mod_pdf
