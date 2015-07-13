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
VERSION=0.1

CFLAGS=-Wextra -Wall -pedantic-errors -O3 -fPIC -fno-common $(ARCH_FLAGS)
LDFLAGS=$(ARCH_FLAGS) 
DEFINES=-D __BLAS_LEGACY__ $(D2_DEFINES)
INCLUDES=-Iinclude/ -I$(MOSEK)/h $(CBLAS_INC)
MOSEKLIB=-L$(MOSEK)/bin -Wl,-rpath,$(MOSEK)/bin $(MOSEK_BIN)
LIBRARIES=-Wl,-rpath,. -Wl,-rpath,$(MOSEK)/bin $(BLAS_LIB) $(OTHER_LIB)

C_SOURCE_FILES=\
	src/d2/clustering.c\
	src/d2/clustering_io.c\
	src/d2/clustering_util.c\
	src/d2/math.c\
	src/utils/blas_like32.c\
	src/utils/blas_like64.c\
	src/d2/centroid_util.c\
	src/d2/centroid_rand.c\
	src/d2/centroid_Bregman.c\
	src/d2/centroid_GradDecent.c\
	src/d2/centroid_ADMM.c\

CPP_SOURCE_FILES=\
	src/d2/solver_mosek.cc\

SOURCE_FILES_WITH_MAIN=\
	src/app/util.cc\
	src/app/main.cc

C_SOURCE_OBJECTS=\
	$(patsubst %.c, %.o, $(C_SOURCE_FILES))\

CPP_SOURCE_OBJECTS=\
	$(patsubst %.cc, %.o, $(CPP_SOURCE_FILES))

LIB_SOURCE_OBJECTS=\
	$(C_SOURCE_OBJECTS)\
	src/d2/solver_mosek.o # use mosek solver as the driver

ALL_OBJECTS=\
	$(C_SOURCE_OBJECTS)\
	$(CPP_SOURCE_OBJECTS)\
	$(patsubst %.cc, %.o, $(SOURCE_FILES_WITH_MAIN))

DEPENDENCY_FILES=\
	$(patsubst %.o, %.d, $(ALL_OBJECTS))

ifeq ($(OS), Darwin)
LIB=\
	libad2c.dylib\
	libmosek64_wrapper.dylib
else
LIB=\
	libad2c.so\
	libmosek64_wrapper.so
endif


all: d2 protein

lib: $(LIB)

%.o: %.c Makefile
	@# Make dependecy file
	$(CC) -MM -MT $@ -MF $(patsubst %.c,%.d,$<) $(CFLAGS) $(DEFINES) $(INCLUDES) $<
	@# Compile
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -g -c -o $@ $<

%.o: %.cc Makefile
	@# Make dependecy file
	$(TCXX) -MM -MT $@ -MF $(patsubst %.cc,%.d,$<) $(CFLAGS) $(DEFINES) $(INCLUDES) $<
	@# Compile
	$(TCXX) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

d2: src/app/util.cc src/app/main.cc $(LIB)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $@ $^ $(LIBRARIES)

data/protein_seq/protein: data/protein_seq/d2_protein_ngram.cc $(LIB)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $@ $^ $(LIBRARIES)


ifeq ($(OS), Darwin)

libmosek64_wrapper.dylib: src/d2/solver_mosek.o
	$(TCXX) -dynamiclib $(DEFINES) $(INCLUDES) -Wl,-install_name,$@ -compatibility_version 7.0 -current_version $(MOSEK_VERSION) -o libmosek64_wrapper.$(MOSEK_VERSION).dylib $< $(MOSEKLIB) $(LIBRARIES) 
	install_name_tool -change  @loader_path/libmosek64.$(MOSEK_VERSION).dylib  $(MOSEK)/bin/libmosek64.$(MOSEK_VERSION).dylib libmosek64_wrapper.$(MOSEK_VERSION).dylib
	ln -sf libmosek64_wrapper.$(MOSEK_VERSION).dylib $@ 

libad2c.dylib: $(C_SOURCE_OBJECTS) libmosek64_wrapper.dylib
	$(CC) -dynamiclib $(DEFINES) $(INCLUDES) -Wl,-install_name,$@ -current_version $(VERSION) -o libad2c.$(VERSION).dylib $^ $(LIBRARIES)
	ln -sf libad2c.$(VERSION).dylib $@

else
libmosek64_wrapper.so: src/d2/solver_mosek.o
	$(TCXX) -shared $(LDFLAGS) $(DEFINES) $(INCLUDES) -Wl,-soname,$@ -o libmosek64_wrapper.$(MOSEK_VERSION).so $< $(MOSEKLIB)
	ln -sf libmosek64_wrapper.$(MOSEK_VERSION).so $@ 

libad2c.so: $(C_SOURCE_OBJECTS) libmosek64_wrapper.so
	$(CC) -shared $(LDFLAGS) $(DEFINES) $(INCLUDES) -Wl,-soname,$@ -o libad2c.$(VERSION).so $^ -Wl,-rpath,. $(LIBRARIES)
	ln -sf libad2c.$(VERSION).so $@ 
endif


-include $(DEPENDENCY_FILES)


protein: data/protein_seq/protein

.PHONY: clean test
clean:
	@rm -f *.so
	@rm -f d2 data/protein_seq/protein
	@for pattern in '*.o' '*.d'; do \
		find . -name "$$pattern" | xargs rm; \
	done

ifeq ($(MPI),1)
test:
	cd test && ./test_mpi.sh
else
test:
	cd test && ./test_single.sh
endif



