module chardatamod

use iso_fortran_env

implicit none

public  :: mean
public  :: stdev
public  :: lambda
public  :: boxcox
private :: llik

integer, parameter :: i4 = int32
integer, parameter :: sp = real32
integer, parameter :: dp = real64

! ------------------------------
! interfaces for overloaded functions for the statistical routines for creating z-scores
! NB the lambda function is not overloaded because it only functions properly with values in real64

interface mean
  module procedure mean_sp,mean_dp
end interface mean

interface stdev
  module procedure stdev_sp,stdev_dp
end interface stdev

interface boxcox
  module procedure boxcox_sp,boxcox_dp
end interface boxcox

contains

! ------------------------------

real(sp) function mean_sp(vals)

implicit none

real(sp), dimension(:), intent(in) :: vals

integer :: n

! ---

n = size(vals)

mean_sp = sum(vals) / real(n)

end function mean_sp

! ------------------------------

real(dp) function mean_dp(vals)

implicit none

real(dp), dimension(:), intent(in) :: vals

integer(i4) :: n

! ---

n = size(vals)

mean_dp = sum(vals) / real(n)

end function mean_dp

! ------------------------------

real(sp) function stdev_sp(vals)

implicit none

real(sp), dimension(:), intent(in) :: vals

real(sp) :: vbar
real(sp) :: n

! -----

n = 1. / real(size(vals))

vbar = sum(vals) * n

stdev_sp = sqrt(n * sum((vals - vbar)**2))

end function stdev_sp

! ------------------------------

real(dp) function stdev_dp(vals)

implicit none

real(dp), dimension(:), intent(in) :: vals

real(dp) :: vbar
real(dp) :: n

! -----

n = 1. / real(size(vals))

vbar = sum(vals) * n

stdev_dp = sqrt(n * sum((vals - vbar)**2))

end function stdev_dp

! ------------------------------

real(dp) function lambda(vals,rng)

! estimate the lambda parameter used in the Box-Cox transformation
! using a maximum likelihood with simple bisection in the range of (-rng,rng)

implicit none

real(dp), dimension(:), intent(in) :: vals
real(dp)              , intent(in) :: rng

! rng is the range over which to estimate lambda. In edge cases, lambda will converge on one of the boundaries

! tolerance parameter, exit the loop when the improvement in likelihood is smaller than this value
real(dp), parameter :: eps = 1.e-8  

! maximum number of iterations
integer, parameter :: imax = 100

real(dp) :: lam0
real(dp) :: lam1
real(dp) :: lam2

real(dp) :: lik0
real(dp) :: lik1
real(dp) :: lik2

real(dp) :: d0
real(dp) :: d1

integer :: i

! -----

lam0 = -rng
lam1 =  0._dp
lam2 =  rng

lik0 = llik(vals,lam0)
lik1 = llik(vals,lam1)
lik2 = llik(vals,lam2)

do i = 1,imax

  d0 = abs(lik1 - lik0)
  d1 = abs(lik2 - lik1)
  
  if (d0 < d1) then  ! go down
    lik2 = lik1
    lam2 = lam1
    lam1 = lam0 + (lam1 - lam0) / 2._dp
  else               ! go up
    lik0 = lik1
    lam0 = lam1
    lam1 = lam1 + (lam2 - lam1) / 2._dp
  end if

  lik1 = llik(vals,lam1)
  
  if (abs(lik1-lik0) < eps .and. abs(lik2-lik1) < eps) then
    lambda = lam1
    return
  end if
  
end do

write(0,*)'estlambda warning: maximum number of iterations reached'

end function lambda

! ------------------------------

real(dp) function llik(vals,lam)

! calculate the likelihood metric for lambda given values and an input lambda

implicit none

real(dp), dimension(:), intent(in) :: vals
real(dp), intent(in) :: lam

real(dp) :: n

real(dp), dimension(size(vals)) :: bcvals
real(dp), dimension(size(vals)) :: zeta
real(dp), dimension(size(vals)) :: vlog

real(dp) :: vbar
real(dp) :: vexp
real(dp) :: zbar
real(dp) :: rss

! -----

n = real(size(vals))

vlog = log(vals)

vbar = mean(vlog)
vexp = exp(vbar)

bcvals = boxcox(lam,vals)

zeta = bcvals / (vexp**(lam - 1._dp))

zbar = mean(zeta)

rss = sum((zeta - zbar)**2)

llik = -0.5_dp * n * log(rss)

end function llik

! ------------------------------

function boxcox_sp(lambda,vals)

implicit none

real(sp), intent(in) :: lambda
real(sp), dimension(:), intent(in) :: vals

real(sp), dimension(size(vals)) :: boxcox_sp

! -----

if (lambda /= 0._sp) then

  boxcox_sp = (vals**lambda - 1._sp) / lambda

else

  boxcox_sp = log(vals)

end if

end function boxcox_sp

! ------------------------------

function boxcox_dp(lambda,vals)

implicit none

real(dp), intent(in) :: lambda
real(dp), dimension(:), intent(in) :: vals

real(dp), dimension(size(vals)) :: boxcox_dp

! -----

if (lambda /= 0._dp) then

  boxcox_dp = (vals**lambda - 1._dp) / lambda

else

  boxcox_dp = log(vals)

end if

end function boxcox_dp

! ------------------------------

end module chardatamod