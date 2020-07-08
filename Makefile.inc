USEPGPLOT = -DPGPLOT
PGPLIBS   =  -lpgplot -lcpgplot

CC=gcc
F77=gfortran
MPICC=mpicc
MPIF77=mpif90

CFLAGS=-std=c99 -O2 ${USEPGPLOT}
FFLAGS=-O2 -g


