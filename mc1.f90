program mc1
  implicit none

  integer, parameter :: nmc = 100       ! number of trials
  integer, parameter :: n = 100000      ! sample size
  real, dimension(nmc) :: mu            ! result of each trial
  real :: mean, stdev                   ! mean and standard deviation
  integer :: j

  !$OMP PARALLEL DO
  do j = 1, nmc
     print *, 'Experiment', j
     mu(j) = monte_carlo(n)
  end do
  !$OMP END PARALLEL DO

  mean = sum(mu) / dble(nmc)
  stdev = sqrt( sum( (mu - mean)**2 ) ) / dble(nmc)

  print *, 'mean', mean
  print *, 'stdev', stdev

contains

  ! Draws n Uniform(0,1) random numbers and returns the sample mean
  function monte_carlo(n) result(y)
    integer, intent(in) :: n
    integer :: i
    real :: y, u

    y = 0.d0
    do i = 1, n
       call random_number(u)            ! draw from Uniform(0,1)
       y = y + u                        ! sum the draws
    end do
    y = y / dble(n)                     ! return the sample mean
  end function monte_carlo

end program mc1
