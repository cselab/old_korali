# compilers and flags

use_torc?=1

CC := mpicc
LD := mpicc

ifeq ($(use_torc),1)
	CC = mpicc
	CFLAGS += -D_USE_TORC_=1 `torc_cflags`
	LDLIBS += `torc_libs`
endif

CFLAGS += -O3 -std=c99
CFLAGS += -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += -Wall -Wno-unused-function
CFLAGS += `gsl-config --cflags`

LDLIBS += `gsl-config --libs`  -lm -lpthread


COMPILE.c = $(CC) $(CFLAGS) -c -o $@
LINK.o    = $(LD) $(LDFLAGS) -o $@

# rules

%.o: %.c
	$(COMPILE.c) $<

%.a:; ar rcs $@ $^

