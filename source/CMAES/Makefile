include ../make/common.mk

FITFUN_OBJS ?= ./fitfun.o
FITFUN_DIR  ?= ../likelihoods
TARGET      ?= cma

all: $(TARGET)

L=../lib

DIR_PRIORS = $L/priors
DIR_CMAES  = $L/cmaes

LDLIBS += -L$(DIR_PRIORS) -lpriors
CFLAGS += -I$(DIR_PRIORS)

LDLIBS += -L$(DIR_CMAES) -lcmaes
CFLAGS += -I$(DIR_CMAES)

CFLAGS += -I$(FITFUN_DIR)

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_CMAES); make libcmaes.a)
	$(LINK.o) $^ $(LDLIBS)

# dependencies

$(TARGET): engine_cmaes.o $(FITFUN_OBJS)
engine_cmaes.o: engine_cmaes.c

# optional: only when FITFUN_OBJS not specified
fitfun.o: fitfun.c

clean:; rm -rf  $(TARGET) *.o
clear:;	rm -rf .tags* tags *.txt *.dat actparcmaes.par *.err

test: cma
	atest `find test -name main`

.PHONY: clear clean all test
