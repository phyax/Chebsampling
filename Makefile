MPIFC = mpif90

# Intel compiler
FC90 = ifort
CC = icc

OPTS90 = -O3 -r8
CCOPTS = -O3 -std=c99

# # GNU compiler
# FC90 = gfortran
# CC = gcc
# 
# OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
# CCOPTS = -O3 -std=c99


# linkage rules
chebsampling : main.o distr.o invsampling.o libmfft1.o \
	ppush2.o pplib2.o dtimer.o
	$(MPIFC) $(OPTS90) -o chebsampling main.o invsampling.o distr.o \
		libmfft1.o libmfft1_h.o modmfft1.o ppush2.o pplib2.o dtimer.o	

# compilation rules
dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

pplib2.o : pplib2.f90
	$(MPIFC) $(OPTS90) -o pplib2.o -c pplib2.f90

ppush2.o : ppush2.f90 invsampling.o
	$(FC90) $(OPTS90) -o ppush2.o -c ppush2.f90

libmfft1.o : libmfft1.f
	$(FC90) $(OPTS90) -o libmfft1.o -c libmfft1.f

libmfft1_h.o : libmfft1_h.f90
	$(FC90) $(OPTS90) -o libmfft1_h.o -c libmfft1_h.f90

modmfft1.o : modmfft1.f90 libmfft1_h.o
	$(FC90) $(OPTS90) -o modmfft1.o -c modmfft1.f90

distr.o : distr.f90
	$(FC90) $(OPTS90) -o distr.o -c distr.f90

invsampling.o : invsampling.f90 modmfft1.o
	$(FC90) $(OPTS90) -o invsampling.o -c invsampling.f90

main.o : main.f90 invsampling.o pplib2.o ppush2.o distr.o
	$(FC90) $(OPTS90) -o main.o -c main.f90

clean :
	rm -f *.o *.mod

clobber : clean
	rm -f chebsampling
