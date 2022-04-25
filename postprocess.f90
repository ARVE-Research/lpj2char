program postprocess

! compile with Makefile

! calculate Z-score transformed burned area fraction and fire carbon emissions from LPJ-LMfire output
! and write output to a pre-generated netCDF file

! base period for transformation is selected by setting the parameters y0 and y1 [0-6000 BP] below
! for comparison with GCD records, these need to be the same as specified in the paleofire pfTransform() call

use iso_fortran_env
use netcdf
use chardatamod, only : mean,stdev,boxcox,estlambda => lambda

implicit none

integer, parameter :: sp = real32
integer, parameter :: dp = real64

integer, parameter :: y0 = 0
integer, parameter :: y1 = 6000

! Following GCD convention, the range of lambda is limited to (-2,2) and alpha is 0.01,
! but both of these values could be customized here
real(dp), parameter :: lrng  = 2.     ! search range for lambda
real(dp), parameter :: alpha = 0.01   ! offset value for raw data to eliminate zeroes

integer :: ncid
integer :: dimid
integer :: varid
integer :: status

integer :: xlen
integer :: ylen
integer :: tlen
integer :: olen

integer :: x
integer :: y

integer :: t0
integer :: t1

real(sp) :: bfmin
real(sp) :: bfmax

real(dp) :: lambda

real(sp) :: missing

character(200) :: lpjfile
character(200) :: outfile

integer, dimension(1) :: pos

real(sp), allocatable, dimension(:) :: lon
real(sp), allocatable, dimension(:) :: lat
integer,  allocatable, dimension(:) :: time

real(sp), allocatable, target, dimension(:,:,:) :: burnedf
real(sp), allocatable, target, dimension(:,:,:) :: acflux_fire

real(sp), pointer, dimension(:) :: vals

real(dp), allocatable, dimension(:) :: mmtransform
real(sp), allocatable, dimension(:) :: bctransform
real(sp), allocatable, dimension(:) :: bctransform_minmax

real(sp), allocatable, dimension(:,:,:) :: bctransform_zscore

real(sp) :: bctransform_minmax_mean

real(sp), dimension(2) :: rng

! -------------------------

call getarg(1,lpjfile)

status = nf90_open(lpjfile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inq_dimid(ncid,'lon',dimid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=xlen)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inq_dimid(ncid,'lat',dimid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=ylen)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inq_dimid(ncid,'time',dimid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (status /= nf90_noerr) call netcdf_err(status)

! -------------------------

allocate(lon(xlen))

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_get_var(ncid,varid,lon)
if (status /= nf90_noerr) call netcdf_err(status)

allocate(lat(ylen))

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_get_var(ncid,varid,lat)
if (status /= nf90_noerr) call netcdf_err(status)

allocate(time(tlen))

status = nf90_inq_varid(ncid,'time',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_get_var(ncid,varid,time)
if (status /= nf90_noerr) call netcdf_err(status)

pos = minloc(abs(y0-time))

t0 = pos(1)

pos = minloc(abs(y1-time))

t1 = pos(1)

! write(0,'(a,i6,a,i6)')'averaging period:',y0,' to',y1

! -------------------------

allocate(burnedf(xlen,ylen,tlen))
allocate(acflux_fire(xlen,ylen,tlen))

status = nf90_inq_varid(ncid,'burnedf',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_get_var(ncid,varid,burnedf)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inq_varid(ncid,'acflux_fire',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_get_var(ncid,varid,acflux_fire)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call netcdf_err(status)

! -------------------------

call getarg(2,outfile)

status = nf90_open(outfile,nf90_write,ncid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_put_var(ncid,varid,lon)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_put_var(ncid,varid,lat)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_inq_varid(ncid,'time',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_put_var(ncid,varid,time(t0:t1))
if (status /= nf90_noerr) call netcdf_err(status)

! -------------------------

olen = 1+t1-t0

allocate(mmtransform(olen))

allocate(bctransform(olen))
allocate(bctransform_minmax(olen))

allocate(bctransform_zscore(xlen,ylen,olen))

! -------------------------
! calculations for burned area fraction

status = nf90_inq_varid(ncid,'zscore_burnedf',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_get_att(ncid,varid,'missing_value',missing)
if (status /= nf90_noerr) call netcdf_err(status)

bctransform_zscore = missing

do y = 1,ylen
  do x = 1,xlen

    vals => burnedf(x,y,t0:t1)
    
    if (count(vals > 0.) == 0) cycle
    
    ! (1) min-max
  
    bfmin = minval(vals)
    bfmax = maxval(vals)

    mmtransform = (vals - bfmin) / (bfmax - bfmin)
    
    ! Following Bart's (2022) advice, lambda is estimated after adding alpha.

    mmtransform = mmtransform + alpha

    ! NB In cases with relatively few fire observations relative to zero background,
    ! lambda ends up being on the boundary of the search domain. 
    
    lambda = estlambda(mmtransform,lrng)

    ! (2) box-cox

    bctransform = sngl(boxcox(lambda,mmtransform))
    
    ! (3) min-max (again)
    
    bfmin = minval(bctransform)
    bfmax = maxval(bctransform)
    
    bctransform_minmax = (bctransform - bfmin) / (bfmax - bfmin)
    
    ! (4) Z-score
    
    bctransform_minmax_mean = mean(bctransform_minmax)
    
    bctransform_zscore(x,y,:) = (bctransform_minmax - bctransform_minmax_mean) / stdev(bctransform_minmax)

  end do
end do

status = nf90_put_var(ncid,varid,bctransform_zscore)
if (status /= nf90_noerr) call netcdf_err(status)

rng(1) = minval(bctransform_zscore,mask=bctransform_zscore/=missing)
rng(2) = maxval(bctransform_zscore,mask=bctransform_zscore/=missing)

status = nf90_put_att(ncid,varid,'actual_range',rng)
if (status /= nf90_noerr) call netcdf_err(status)

! -------------------------
! calculations for fire carbon emissions

status = nf90_inq_varid(ncid,'zscore_acflux',varid)
if (status /= nf90_noerr) call netcdf_err(status)

status = nf90_get_att(ncid,varid,'missing_value',missing)
if (status /= nf90_noerr) call netcdf_err(status)

bctransform_zscore = missing

do y = 1,ylen
  do x = 1,xlen

    vals => acflux_fire(x,y,t0:t1)
    
    if (count(vals > 0.) == 0) cycle
    
    ! (1) min-max
  
    bfmin = minval(vals)
    bfmax = maxval(vals)

    mmtransform = (vals - bfmin) / (bfmax - bfmin)
    
    ! Following Bart's (2022) advice, lambda is estimated after adding alpha.

    mmtransform = mmtransform + alpha

    ! NB In cases with relatively few fire observations relative to zero background,
    ! lambda ends up being on the boundary of the search domain. 
    
    lambda = estlambda(mmtransform,lrng)

    ! (2) box-cox

    bctransform = sngl(boxcox(lambda,mmtransform))
    
    ! (3) min-max (again)
    
    bfmin = minval(bctransform)
    bfmax = maxval(bctransform)
    
    bctransform_minmax = (bctransform - bfmin) / (bfmax - bfmin)
    
    ! (4) Z-score
    
    bctransform_minmax_mean = mean(bctransform_minmax)
    
    bctransform_zscore(x,y,:) = (bctransform_minmax - bctransform_minmax_mean) / stdev(bctransform_minmax)

  end do
end do

status = nf90_put_var(ncid,varid,bctransform_zscore)
if (status /= nf90_noerr) call netcdf_err(status)

rng(1) = minval(bctransform_zscore,mask=bctransform_zscore/=missing)
rng(2) = maxval(bctransform_zscore,mask=bctransform_zscore/=missing)

status = nf90_put_att(ncid,varid,'actual_range',rng)
if (status /= nf90_noerr) call netcdf_err(status)

! -------------------------

status = nf90_close(ncid)
if (status /= nf90_noerr) call netcdf_err(status)

! -------------------------

contains

subroutine netcdf_err(ncstat)

! handles netCDF errors; writes the error message if an error is encountered.

use netcdf, only : nf90_strerror

implicit none

integer, intent(in) :: ncstat

! ------

write(0,'(a,i5,a,a)')' NetCDF error ',ncstat,' encountered: ',trim(nf90_strerror(ncstat))
stop

end subroutine netcdf_err

! -------------------------

end program postprocess
