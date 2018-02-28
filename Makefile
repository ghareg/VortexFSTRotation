CC = icpc #g++
CFLAGS = -Wall -O3 -std=c++11 -fopenmp #-gdwarf-2
MKLINC = -I${MKL_ROOT}/include
CPLUSINC = -I/share/usr/libs/gsl/2.3/intel/include
LFGSL = -L/share/usr/libs/gsl/2.3/intel/lib -lgsl -lgslcblas
LFMKL = -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -liomp5 -L/global/u/areg.ghazaryan/usr/lib -lsuperlu_5.1
LFARP = -L/share/usr/arpack/3.1.5_intel/lib -larpack
LFLAGS = -lpthread -lgsl -lgslcblas

SOURCES = main.cpp \
		  basis.cpp \
		  matrix.cpp \
		  operator.cpp \
		  diag.cpp

OBJS = ${SOURCES:.cpp=.o}


all: $(OBJS)
	$(CC) $(OBJS) -o vortexUR $(LFGSL) $(LFMKL) $(LFARP)
	rm -f *.o


$(OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) $(MKLINC) $(CPLUSINC) $< -o $@

clean:
	rm -f *.o vortexUR*

