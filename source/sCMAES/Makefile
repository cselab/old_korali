include ../make/common.mk

FITFUN_OBJS ?= ./fitfun.o
FITFUN_DIR  ?= ../likelihoods
TARGET      ?= scma

all: $(TARGET)

DIR_PRIORS=../priors
DIR_CMAES=../lib/cmaes

LDLIBS += -L$(DIR_PRIORS) -lpriors
CFLAGS += -I$(DIR_PRIORS)

LDLIBS += -L$(DIR_CMAES) -lcmaes
CFLAGS += -I$(DIR_CMAES)

CFLAGS += -I$(FITFUN_DIR)

CXXFLAGS += `pkg-config libgp --cflags`
LDLIBS += `pkg-config libgp --libs`

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_CMAES); make libcmaes.a)
	$(LINK.o.xx) $^ $(LDLIBS)

# dependencies

$(TARGET): engine_scmaes.o surrogate_samples.o surrogate.o $(FITFUN_OBJS)
engine_cmaes.o: engine_cmaes.c
engine_scmaes.o: engine_scmaes.c
surrogate.o: surrogate.cpp surrogate.h
surrogate_samples.o: surrogate_samples.c surrogate_samples.h

# optional: only when FITFUN_OBJS not specified
fitfun.o: fitfun.c

clean:; rm -rf  $(TARGET) *.o
clear:;	rm -rf .tags* tags *.txt *.dat actparcmaes.par *.err

test: scma
	atest `find test -name main`

.PHONY: clear clean all test
