include make.inc

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

ARCH_FLAGS=-m64
CFLAGS=-Wextra -Wall -pedantic-errors -O3 $(ARCH_FLAGS)
LDFLAGS=$(ARCH_FLAGS)
DEFINES=-D __BLAS_LEGACY__
INCLUDES=-Iinclude/ -I$(MOSEK)/h $(CBLAS_INC)
LIBRARIES=-L$(MOSEK)/bin -Wl,-rpath,$(MOSEK)/bin -lmosek64 -lpthread $(BLAS_LIB) $(OTHER_LIB)


C_SOURCE_FILES=\
	src/d2_clustering.c\
	src/d2_clustering_io.c\
	src/d2_clustering_util.c\
	src/d2_math.c\
	src/blas_like.c\
	src/d2_centroid_util.c\
	src/d2_centroid_rand.c\
	src/d2_centroid_Bregman.c\
	src/d2_centroid_GradDecent.c\
	src/d2_centroid_ADMM.c\

CPP_SOURCE_FILES=\
	src/d2_solver_mosek.cc\
	src/util.cc

SOURCE_FILES_WITH_MAIN=\
	src/main.cc

C_SOURCE_OBJECTS=\
	$(patsubst %.c, %.o, $(C_SOURCE_FILES))\

CPP_SOURCE_OBJECTS=\
	$(patsubst %.cc, %.o, $(CPP_SOURCE_FILES))

LIB_SOURCE_OBJECTS=\
	$(C_SOURCE_OBJECTS)\
	src/d2_solver_mosek.o # use mosek solver as the driver

ALL_OBJECTS=\
	$(C_SOURCE_OBJECTS)\
	$(CPP_SOURCE_OBJECTS)\
	$(patsubst %.cc, %.o, $(SOURCE_FILES_WITH_MAIN))

DEPENDENCY_FILES=\
	$(patsubst %.o, %.d, $(ALL_OBJECTS))

all: d2 protein

%.o: %.c Makefile
	@# Make dependecy file
	$(CC) -MM -MT $@ -MF $(patsubst %.c,%.d,$<) $(CFLAGS) $(DEFINES) $(INCLUDES) $<
	@# Compile
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

%.o: %.cc Makefile
	@# Make dependecy file
	$(CXX) -MM -MT $@ -MF $(patsubst %.cc,%.d,$<) $(CFLAGS) $(DEFINES) $(INCLUDES) $<
	@# Compile
	$(CXX) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

-include $(DEPENDENCY_FILES)

ifeq ($(OS), Darwin)
d2: $(ALL_OBJECTS)
	$(CXX) $(LDFLAGS) $(DEFINES) -o $@ $(ALL_OBJECTS) $(LIBRARIES)
	install_name_tool -change  @loader_path/libmosek64.$(MOSEK_VERSION).dylib  $(MOSEK)/bin/libmosek64.$(MOSEK_VERSION).dylib $@

data/protein_seq/protein: data/protein_seq/d2_protein_ngram.cc $(LIB_SOURCE_OBJECTS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $@ $^ $(LIBRARIES)
	install_name_tool -change  @loader_path/libmosek64.$(MOSEK_VERSION).dylib  $(MOSEK)/bin/libmosek64.$(MOSEK_VERSION).dylib $@

else 
d2: $(ALL_OBJECTS)
	$(CXX) $(LDFLAGS) $(DEFINES) -o $@ $(ALL_OBJECTS) $(LIBRARIES)
data/protein_seq/protein: data/protein_seq/d2_protein_ngram.cc $(LIB_SOURCE_OBJECTS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $@ $^ $(LIBRARIES)
endif


protein: data/protein_seq/protein

.PHONY: clean
clean:
	@rm -f *_test d2 data/protein_seq/protein
	@for pattern in '*.o' '*.d'; do \
		find . -name "$$pattern" | xargs rm; \
	done
