# source/make/common.mk

# compilers and flags

use_torc?=0
use_omp?=0

CC  := gcc
LD  := gcc
CXX := g++-8
LDX := g++-8

ifeq ($(use_torc),1)
	CFLAGS   += -D_USE_TORC_=1 `torc_cflags` 
	CXXFLAGS += -D_USE_TORC_=1 `torc_cflags` -fpermissive 
	LDLIBS   += `torc_libs`
	LDXLIBS  += `torc_libs` -ltorc
endif

ifeq ($(use_omp),1)
    CFLAGS += -D_USE_OPENMP_=1 -fopenmp
    LDLIBS += -fopenmp
	# (TODO: what about c++?)
endif


CFLAGS += -O3 -std=c99
CFLAGS += -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += -Wall -Wno-unused-function
CFLAGS += `gsl-config --cflags`

CXXFLAGS += -std=c++14
#CXXFLAGS += -O3
CXXFLAGS += -g -DDEBUG #remove debug if done and add -O3
CXXFLAGS += -Wall -Wno-unused-function
CXXFLAGS += `gsl-config --cflags`

LDLIBS  += `gsl-config --libs` -lm -lpthread
LDXLIBS += `gsl-config --libs` -lm -lpthread

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
