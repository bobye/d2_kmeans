CC=gcc -std=c99
CXX=g++

ARCH_FLAGS=
CFLAGS=-Wextra -Wall -pedantic-errors $(ARCH_FLAGS) -O3
LDFLAGS=$(ARCH_FLAGS)
DEFINES=
INCLUDES=-Iinclude/
LIBRARIES=


C_SOURCE_FILES=\
	src/d2_clustering.c

CPP_SOURCE_FILES=


SOURCE_FILES_WITH_MAIN=src/main.cc

SOURCE_OBJECTS=\
	$(patsubst %.c, %.o, $(C_SOURCE_FILES))\
	$(patsubst %.c, %.o, $(CPP_SOURCE_FILES))

ALL_OBJECTS=\
	$(SOURCE_OBJECTS) \
	$(patsubst %.cc, %.o, $(SOURCE_FILES_WITH_MAIN))

DEPENDENCY_FILES=\
	$(patsubst %.o, %.d, $(ALL_OBJECTS))

all: test

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

test: $(ALL_OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(ALL_OBJECTS) $(LIBRARIES)

.PHONY: clean
clean:
	@rm test
	@for pattern in '*.o' '*.d'; do \
		find . -name "$$pattern" | xargs rm; \
	done
