include ../../make.inc

ifndef MPI
MPI=1
endif

ifeq ($(MPI),1)
CC=$(MPICC)
CXX=$(MPICXX)
else
CC=$(TCC)
CXX=$(TCXX)
endif

OS=$(shell uname)

CFLAGS=-Wextra -Wall -pedantic-errors -O3 -fPIC -fno-common $(ARCH_FLAGS)
LDFLAGS=$(ARCH_FLAGS)
DEFINES=-D __BLAS_LEGACY__ $(D2_DEFINES)
INCLUDES=-I../../include/ -I$(MOSEK)/h $(CBLAS_INC)
MOSEKLIB=-L$(MOSEK)/bin -Wl,-rpath,$(MOSEK)/bin $(MOSEK_BIN)
LIBRARIES=-Wl,-rpath,../../ -Wl,-rpath,$(MOSEK)/bin $(BLAS_LIB) $(OTHER_LIB)



ifeq ($(OS), Darwin)
LIB=\
	../../libad2c.dylib\
	../../libmosek64_wrapper.dylib
else
LIB=\
	../../libad2c.so\
	../../libmosek64_wrapper.so
endif

all: protein

protein: d2_protein_ngram.cc $(LIB)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $@ $^ $(LIBRARIES)


clean:
	@rm -f protein
