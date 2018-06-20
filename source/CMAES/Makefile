include ../make/common.mk

FITFUN_OBJS ?= ./fitfun.o
FITFUN_DIR  ?= ../likelihoods
TARGET      ?= cma

all: $(TARGET)

DIR_PRIORS=../priors
LDLIBS += -L$(DIR_PRIORS) -lpriors
CFLAGS += -I$(FITFUN_DIR)

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	$(LINK.o) $^ $(LDLIBS)

# dependencies

$(TARGET): cmaes.o engine_cmaes.o $(FITFUN_OBJS)
cmaes.o: cmaes.c cmaes.h cmaes_interface.h
engine_cmaes.o: engine_cmaes.c cmaes.h cmaes_interface.h

# optional: only when FITFUN_OBJS not specified
fitfun.o: fitfun.c


clean:; rm -rf  $(TARGETS) *.o
clear:;	rm -rf .tags* tags *.txt *.dat actparcmaes.par *.err

test: cma
	atest `find test -name main`

.PHONY: clear clean all test
