# compilers and flags

use_torc?=0

CC := mpicc
LD := mpicc
CXX := mpic++
LDX := mpic++

ifeq ($(use_torc),1)
	CFLAGS += -D_USE_TORC_=1 `torc_cflags`
	LDLIBS += `torc_libs`
endif

CFLAGS += -O3 -std=c99
CFLAGS += -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += -Wall -Wno-unused-function
CFLAGS += `gsl-config --cflags`

LDLIBS += `gsl-config --libs`  -lm -lpthread

CXXFLAGS += -O3
CXXFLAGS += -Wall

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
