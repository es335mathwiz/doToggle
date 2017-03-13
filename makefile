#identify operating system
UNAME= $(shell uname)
NUWEBFLAGS = -t
SPAMADIR=../sparseAMA

ifeq ($(UNAME),Linux)
#compilers
CC = gcc
FCFLAGS = -c -O2 -I $(SPAMADIR)/src/main/include   -I /msu/res5/software/myUsr/include/  -I../CStochSims/
FCFLAGS = -c -g -Wall -I $(SPAMADIR)/src/main/include   -I /msu/res5/software/myUsr/include  -I../CStochSims/
#lapack
LAPACKLIBS=   -L /msu/res5/software/ARPACK96forCluster -larpack_linux -L/msu/res5/software/lapackGithubForCluster -llapack -lrefblas
CUNITLIBS= -L /msu/res5/software/myUsr/lib/ -l cunit
endif

ifeq ($(UNAME),Darwin)
#compilers
CC = gcc-6
FCFLAGS = -c -O2 -I$(SPAMADIR)/src/main/include   -I /Users/garyanderson/myUsr/include/  -I../CStochSims/
FCFLAGS = -c -Wall -g -I $(SPAMADIR)/src/main/include   -I /Users/garyanderson/myUsr/include/  -I../CStochSims/
#lapack
LAPACKLIBS=  -L /Users/garyanderson/ARPACK96/  -larpack_MACOS -L /Users/garyanderson/lapack-release/ -llapack -lrefblas
CUNITLIBS= -L /Users/garyanderson/myUsr/lib -l cunit
endif

echo:
	@echo $(FCFLAGS)

#compilers
FC = gfortran
SPARSEAMALIB= -L../sparseAMA -lsparseAMA
STOCHSIMSLIB= -L../CStochSims -lstochSims

.c.o:
	$(CC) $(FCFLAGS) -c $*.c

.f.o:
	$(FC) $(FCFLAGS) -c $*.f


rbcTryC.o:	rbcTryC.c
	$(CC) $(FCFLAGS)  -o rbcTryC.o rbcTryC.c
rbcTryCDrv.o:	rbcTryCDrv.c
	$(CC) $(FCFLAGS) -o rbcTryCDrv.o rbcTryCDrv.c
rbcTryCSupport.o:	rbcTryCSupport.c
	$(CC) $(FCFLAGS) -o rbcTryCSupport.o rbcTryCSupport.c
rbcTryCData.o:	rbcTryCData.c
	$(CC) $(FCFLAGS) -o rbcTryCData.o rbcTryCData.c
rbcTryCShocks.o:	rbcTryCShocks.c
	$(CC) $(FCFLAGS) -o rbcTryCShocks.o rbcTryCShocks.c


runrbcTryC.o:	runrbcTryC.c runrbcTryCLocalDefs.h
	$(CC) $(FCFLAGS) runrbcTryC.c 



runrbcTryC: runrbcTryC.o	rbcTryC.o \
	rbcTryCDrv.o rbcTryCSupport.o \
	rbcTryCData.o rbcTryCShocks.o
	$(FC) -o runrbcTryC -g runrbcTryC.o  rbcTryC.o \
	rbcTryCDrv.o rbcTryCSupport.o \
	rbcTryCData.o rbcTryCShocks.o  $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS)

clean: 
	rm -f *.o rbcTryC libstochSims
