USEPGPLOT = -DPGPLOT
PGPLIBS   =  -lpgplot -lcpgplot

CC=icc
F77=ifort
MPICC=mpicc
MPIF77=mpif90

CFLAGS=-std=c99 -O2 ${USEPGPLOT}
FFLAGS=-nofor_main -O2 -g


