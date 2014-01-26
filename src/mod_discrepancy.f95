! Module that implements methods to calculate several types of discrepancy.
!
! References:
!
!   [1] `A Method for Exact Calculation of the Discrepancy of
!       Low-dimensional Finite Point Sets I'.
!
!   [2] Niederreiter, Harald, `Random Number Generation and Quasi-Monte Carlo
!       methods', 1992, ISBN 0-89871-295-5.
!
!   [3] James, F., Hoogland, J. and Kleiss, R., `Multidimensional sampling
!       for simulation and integration: measures, discrepancies, and
!       and quasi-random numbers', Computer Physics Communications, vol. 99,
!       pages 180-220, 1997.
!
!   [4] Warnock, Tony T., `Computational Investigations of Low-Discrepancy
!       Point Sets' in `Applications of Number Theory to Numerical Analysis',
!       pages 319-343, 1972.
!
!   [5] Hickernell, Fred. J., `A generalized discrepancy and quadrature error
!       bound', Mathematics of Computation, Vol. 67, Nb. 221, January 1998,
!       pages 299-322.
!
!   [6] Joe, Stephen, `An Average L_2 Discrepancy for Number-Theoretic Rules',
!       SIAM Journal on Numerical Analysis, Vol. 36, No. 6 (Sep.-Oct., 1999),
!       1949-1961.
!
! TODO:
!   * implement the discrepancy from formula (5.1c) in [5], which is the same
!     as the discrepancy from formula (1.4) in [6].

module mod_discrepancy

  use numeric_kinds
  use mod_sort

  private

  public :: T_N_star_squared
  public :: T_N_extreme_squared
  public :: D2_hickernell_squared
  public :: d_n_star1d
  public :: d_n_star2d
  public :: d_n_extreme1d
  public :: discd2

contains

  ! Calculate the squared L_2 star-discrepancy of an [sxN] pointset where s
  ! is the dimension of the pointset, and N is the number of points.
  !
  subroutine T_N_star_squared(x, tn2)

    real(kind=qp), dimension(:,:), intent(in) :: x
    real(kind=qp), dimension(:), intent(out)  :: tn2

    integer(kind=i4b)                        :: s, nbpoints, N, i
    real(kind=qp), dimension(:), allocatable :: P, Q, J

    s = size(x, 1)
    nbpoints = size(x, 2)

    allocate(P(nbpoints))
    allocate(Q(nbpoints))
    allocate(J(nbpoints))
      
    P(1) = product(1.0_qp-x(:,1))
    Q(1) = product(1.0_qp-x(:,1)**2)
    J(1) = P(1) - 2.0_qp**(-s+1)*Q(1) + 3.0_qp**(-s)
    
    do N = 2, nbpoints

      P(N) = P(N-1) &
         + 2*sum(product(1-max(x(:,1:N-1),spread(x(:,N),2,N-1)) ,1)) &
         + product(1-x(:,N))

      Q(N) = Q(N-1) + product(1.0_qp-x(:,N)**2)
      J(N) = P(N) - 2.0_qp**(-s+1)*Q(N)*N + 3.0_qp**(-s)*N*N

    end do

    do i=1,nbpoints
      tn2(i) = (J(i)/i)/i
    end do

    deallocate(P)
    deallocate(Q)
    deallocate(J) 

  end subroutine T_N_star_squared


  ! Calculate the squared L_2 extreme-discrepancy of an [sxN] pointset,
  ! where s is the dimension of the pointset and N is the number of points.
  !
  subroutine T_N_extreme_squared(x, tn2)

    real(kind=qp), dimension(:,:), intent(in) :: x
    real(kind=qp), dimension(:), intent(out)  :: tn2

    integer(kind=i4b)                        :: s, nbpoints, N, i
    real(kind=qp), dimension(:), allocatable :: P, Q, J

    s = size(x, 1)
    nbpoints = size(x, 2)

    allocate(P(nbpoints))
    allocate(Q(nbpoints))
    allocate(J(nbpoints))

    P(1) = product( (1-x(:,1))*x(:,1) )
    Q(1) = product( (1-x(:,1))*x(:,1) )
    J(1) = P(1) - 2.0_qp**(-s+1)*Q(1) + 12.0_qp**(-s)

    do N = 2, nbpoints

      P(N) = P(N-1) &
         + 2*sum(product( (1-max(x(:,1:N-1),spread(x(:,N),2,N-1)))      &
                            *min(x(:,1:N-1),spread(x(:,N),2,N-1)) , 1)) &
         + product( (1-x(:,N))*x(:,N) )

      Q(N) = Q(N-1) + product( (1.0_qp-x(:,N))*x(:,N) )
      J(N) = P(N) - 2.0_qp**(-s+1)*Q(N)*N + 12.0_qp**(-s)*N*N

    end do

    do i=1,nbpoints
      tn2(i) = (J(i)/i)/i
    end do

    deallocate(P)
    deallocate(Q)
    deallocate(J)

  end subroutine T_N_extreme_squared


  ! Compute the D^2 discrepancy of an [sxN] pointset as given in [5] (formula ?)
  ! and [6] (formula 1.4).
  !
  ! TODO:
  !   * check for correctness!
  !
  subroutine D2_hickernell_squared(x, D2)

    real(kind=qp), dimension(:,:), intent(in) :: x
    real(kind=qp), dimension(:), intent(out)  :: D2

    integer(kind=i4b)                        :: s, nbpoints, N, i
    real(kind=qp), dimension(:), allocatable :: P, Q, J

    s = size(x, 1)
    nbpoints = size(x, 2)

    allocate(P(nbpoints))
    allocate(Q(nbpoints))
    allocate(J(nbpoints))

    P(1) = product(2.0_qp-x(:,1))
    Q(1) = product(3.0_qp-x(:,1)**2)
    J(1) = P(1) - 2.0_qp**(-s+1)*Q(1) + (4.0_qp/3)**(s)

    do N = 2, nbpoints

      P(N) = P(N-1) &
         + 2*sum( product(2-max(x(:,1:N-1),spread(x(:,N),2,N-1)) ,1) ) &
         + product(2-x(:,N))

      Q(N) = Q(N-1) + product(3.0_qp-x(:,N)**2)
      J(N) = P(N) - 2.0_qp**(-s+1)*Q(N)*N + ((4.0_qp/3)**s)*N*N

    end do

    do i=1,nbpoints
      D2(i) = (J(i)/i)/i
    end do

    deallocate(P)
    deallocate(Q)
    deallocate(J)

  end subroutine D2_hickernell_squared


  ! Calculate the star discrepancy as defined with the infinity
  ! norm (supremum) of a 1-dimensional sequence (See [2], Theorem 2.6 p. 15)
  !
  ! TODO:
  !   * can we change this method so that instead of a single real it returns
  !     an array D_N_STAR of the same size of X and containing at position i the
  !     1D star discrepancy of the points x(1),...,x(i)?
  !
  subroutine d_n_star1d(x, d_n_star)
    real(kind=qp), intent(in), dimension(:) :: x
    real(kind=qp), intent(out)              :: d_n_star 

    integer(kind=i4b)                     :: N
    real(kind=qp), dimension(size(x))     :: x_sorted
    integer(kind=i4b), dimension(size(x)) :: n_small
    integer(kind=i4b)                     :: i

    N = size(x)
    n_small = (/ (i, i=1,N) /)

    ! The elements must be sorted to calculate this type of discrepancy!
    x_sorted = x
    call quicksort(x_sorted)

    d_n_star = 0.5_qp/N &
                 + maxval(abs(x_sorted-(2.0_qp*n_small-1.0_qp)/(2.0_qp*N)))

  end subroutine d_n_star1d


  ! Calculate the L_{\infty} extreme discrepancy of a 1-dimensional sequence.
  ! (See [2], Theorem 2.7 page 16)
  !
  subroutine d_n_extreme1d(x, d_n_extreme)
    real(kind=qp), intent(in), dimension(:) :: x
    real(kind=qp), intent(out)              :: d_n_extreme

    integer(kind=i4b)                     :: N
    real(kind=qp), dimension(size(x))     :: x_sorted
    integer(kind=i4b), dimension(size(x)) :: n_small
    integer(kind=i4b)                     :: i

    N = size(x)
    n_small = (/ (i, i=1,N) /)

    ! The elements must be sorted to calculate this type of discrepancy!
    x_sorted = x
    call quicksort(x_sorted)

    d_n_extreme = 1.0_qp/N + maxval(real(n_small, kind=dp)/N-x_sorted) &
                           - minval(real(n_small, kind=dp)/N-x_sorted)

  end subroutine d_n_extreme1d


  ! See [1], page 117 for the details of the algorithm for this routine.
  !
  subroutine d_n_star2d(x, discr)
    real(kind=qp), dimension(:,:), intent(in) :: x
    real(kind=qp), intent(out)                :: discr

    real(kind=qp), dimension(size(x, 2))                   :: x_sorted
    real(kind=qp), dimension(2,0:size(x, 2)+1)             :: sorted_by_x
    real(kind=qp), dimension(0:size(x, 2), 0:size(x, 2)+1) :: xi
    real(kind=qp), dimension(0:size(x, 2), 0:size(x, 2))   :: mymax
    integer(kind=i4b), dimension(size(x, 2))               :: indices
    integer(kind=i4b)                                      :: i, l, k, N, s

    s = size(x, 1)
    N = size(x, 2)

    if (s == 2) then

      indices = (/ (i, i=1,size(x,2)) /)

      ! First make sure the pointset is sorted according to
      ! its first dimension, and also permute the second dimension
      ! accordingly.
      x_sorted = x(1,:)
      call quicksort(x_sorted, indices)
      sorted_by_x(1,:) = (/ 0.0_qp, x_sorted, 1.0_qp /)
      sorted_by_x(2,:) = (/ 0.0_qp, x(2,indices), 1.0_qp /)

      xi = 0.0_qp
      do l=0,N
        xi(l, 0:l+1) = (/ sorted_by_x(2, 0:l), 1.0_qp /)
        call quicksort(xi(l, 0:l+1))
      end do

      mymax = 0.0_qp
      do l=0,N
        do k=0,l
          mymax(l, k) = & 
             max( abs(real(k, kind=dp)/N-sorted_by_x(1,l)*xi(l,k)),     &
                  abs(real(k, kind=dp)/N-sorted_by_x(1,l+1)*xi(l,k+1)) )
        end do
      end do

      ! Note: it is ok to take the maximum over the whole matrix my_max
      ! instead of just taking the maximum over just the lower triangular
      ! part.  The upper triangular part is zero anyway.
      discr = maxval(mymax)

    else
      print *, "ERROR: dimension of pointset must be 2!!!"
    end if
    
  end subroutine d_n_star2d


  ! Compute the quadratic discrepancy for the point set
  !
  !     x(i,mu),  i=1,2,...,ntot  mu=1,2,...,ns
  !
  ! using the code from [3].
  !
  ! The output is the array d2(n), n=1,2,...,ntot and gives the
  ! discrepancy^2 for the first n points.
  !
  ! Notes:
  !
  !   * The code from [3] was syntactically slightly adapted
  !     to make it F-compliant.  We also return discrepancy^2
  !     instead of discrepancy^2*ntot^2.
  !
  !   * Timing results (see timings directory) indicate that
  !     this version is a bit slower than our T_N_star_squared()
  !     subroutine.
  !
  subroutine discd2(ntot, ns, x, d2)
    integer(kind=i4b), intent(in)             :: ntot
    integer(kind=i4b), intent(in)             :: ns
    real(kind=qp), dimension(:,:), intent(in) :: x
    real(kind=qp), dimension(:), intent(out)  :: d2

    real(kind=qp)     :: c2, c3
    real(kind=qp)     :: an, bn
    real(kind=qp)     :: a, b
    real(kind=qp)     :: temp
    integer(kind=i4b) :: i
    integer(kind=i4b) :: n, mu

    ! Initialize a few constants
    c2 = 0.5_qp**ns
    c3 = (1.0_qp/3)**ns
    bn = 0.0_qp

    ! Loop over the number of points
    do n=1,ntot

      ! Compute b(n) and a(n,n) for the new point
      a = 1.0_qp
      b = 1.0_qp
      do mu=1,ns
        a = a*(1-x(n,mu))
        b = b*(1-x(n,mu)*x(n,mu))
      end do
      b = c2*b

      ! Update running value of sum_b
      bn = bn + b

      if (n==1) then ! Case n=1

        d2(n) = a - 2.0_qp*bn + c3

      else           ! Case n>1

        ! Sum the a(i,n) for i=1 to n-1
        an = 0.0_qp
        do i=1,n-1
          temp = 1.0_qp
          do mu=1,ns
            temp = temp*(1-max(x(i,mu),x(n,mu)))
          end do
          an = an + temp
        end do

        ! Give d2(n) for n>1 by relating it to d2(n-1)
        d2(n) = d2(n-1) + 2*an + a - 2*bn - 2.0_qp*(n-1)*b + (2.0_qp*n-1)*c3

      end if

    end do

    d2 = d2/(/ (i, i=1,ntot) /)
    d2 = d2/(/ (i, i=1,ntot) /)

  end subroutine discd2


end module mod_discrepancy
