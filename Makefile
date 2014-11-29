CC=gcc -std=c99
CXX=g++

MOSEK=$(HOME)/mosek/7/tools/platform/osx64x86

ARCH_FLAGS=
CFLAGS=-Wextra -Wall -pedantic-errors -O3 $(ARCH_FLAGS)
LDFLAGS=$(ARCH_FLAGS) -L$(MOSEK)/bin -lmosek64 -pthread
DEFINES=
INCLUDES=-Iinclude/ -I$(MOSEK)/h
LIBRARIES=-framework accelerate


C_SOURCE_FILES=\
	src/d2_clustering.c\
	src/d2_math.c\
	src/blas_like.c\
	src/d2_centroid_rand.c\
	src/d2_centroid_Bregman.c\
	src/d2_solver_mosek.c

CPP_SOURCE_FILES=\
	src/util.cc

SOURCE_FILES_WITH_MAIN=\
	src/main.cc

SOURCE_OBJECTS=\
	$(patsubst %.c, %.o, $(C_SOURCE_FILES))\
	$(patsubst %.cc, %.o, $(CPP_SOURCE_FILES))

ALL_OBJECTS=\
	$(SOURCE_OBJECTS)\
	$(patsubst %.cc, %.o, $(SOURCE_FILES_WITH_MAIN))

DEPENDENCY_FILES=\
	$(patsubst %.o, %.d, $(ALL_OBJECTS))

all: d2

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

d2: $(ALL_OBJECTS)
	$(CXX) $(LDFLAGS) $(DEFINES) -o $@ $(ALL_OBJECTS) $(LIBRARIES)
	install_name_tool -change  @loader_path/libmosek64.7.0.dylib  @loader_path/../../../mosek/7/tools/platform/osx64x86/bin/libmosek64.7.0.dylib $@

transportation_test: src/transportation_test.c src/d2_solver_mosek.o src/blas_like.o
	$(CC) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $@ $^ $(LIBRARIES)
	install_name_tool -change  @loader_path/libmosek64.7.0.dylib  @loader_path/../../../mosek/7/tools/platform/osx64x86/bin/libmosek64.7.0.dylib $@

.PHONY: clean
clean:
	@rm -f *_test d2
	@for pattern in '*.o' '*.d'; do \
		find . -name "$$pattern" | xargs rm; \
	done
