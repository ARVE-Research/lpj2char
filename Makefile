# Makefile for the lpj2char program

FC=gfortran
FCFLAGS=-ffree-form -ffree-line-length-none

# type nf-config --all to find out the path to your netcdf library
netcdf=/home/public/easybuild/software/netCDF-Fortran/4.5.4-gompi-2021b

NC_LIB=$(netcdf)/lib
NC_INC=$(netcdf)/include

CPPFLAGS = -I$(NC_INC)
LDFLAGS  = -L$(NC_LIB)
LIBS     = -lnetcdff
FCFLAGS  = # -ffree-line-length-none -finit-local-zero -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal -g -fbacktrace -Wall -pedantic

#---------------------------------------------

OBJS = chardatamod.f90 \
			 postprocess.f90

#---------------------------------------------

.SUFFIXES: .o .f90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all::	postprocess

postprocess: $(OBJS)
	$(FC) $(FCFLAGS) -o lpj2char $(OBJS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

clean::
	-rm lpj2char *.mod *.o
