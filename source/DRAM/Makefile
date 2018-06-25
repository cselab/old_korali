include ../make/common.mk

FITFUN_OBJS ?= ./fitfun.o
FITFUN_DIR  ?= ../likelihoods
TARGET      ?= dram

all: $(TARGET)

DIR_PRIORS=../priors
LDLIBS += -L$(DIR_PRIORS) -lpriors
CFLAGS += -I$(FITFUN_DIR)

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	$(LINK.o) $^ $(LDLIBS)

# dependencies

$(TARGET): dram.o run_dram.o $(FITFUN_OBJS)
dram.o: dram.c dram.h
run_dram.o: run_dram.c dram.h

# optional: only when FITFUN_OBJS not specified
fitfun.o: fitfun.c

clean:; rm -rf  $(TARGET) *.o
clear:;	rm -rf .tags* tags *.txt *.dat actparcmaes.par *.err

.PHONY: clear clean all
