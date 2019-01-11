# compilers and flags

use_torc?=0
use_omp?=0

CC  := gcc
LD  := gcc
CXX := g++
LDX := g++

ifeq ($(use_torc),1)
	CC     := mpicc -fopenmp
	CFLAGS += -D_USE_TORC_=1 `torc_cflags`
	LDLIBS += `torc_libs`
endif

ifeq ($(use_omp),1)
    CFLAGS += -D_USE_OPENMP_=1 -fopenmp
    LDLIBS += -fopenmp
endif


CFLAGS += -O3 -std=c99
CFLAGS += -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += -Wall -Wno-unused-function
CFLAGS += `gsl-config --cflags`

CXXFLAGS += -O3 -std=c++11
CXXFLAGS += -Wall -Wno-unused-function

LDLIBS += `gsl-config --libs`  -lm -lpthread


COMPILE.c   = $(CC)  $(CFLAGS)   -c -o $@
COMPILE.cxx = $(CXX) $(CXXFLAGS) -c -o $@
LINK.o      = $(LD)  $(LDFLAGS)     -o $@
LINK.o.xx   = $(LDX) $(LDFLAGS)     -o $@


# rules

%.o: %.c
	$(COMPILE.c) $<

%.o: %.cpp
	$(COMPILE.cxx) $<

%.a:
	ar rcs $@ $^
	ranlib $@
