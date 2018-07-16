include ../make/common.mk

LIB_FITFUN ?= fitfun_tests
DIR_FITFUN ?= ../likelihoods
TARGET      ?= scma

all: $(TARGET)

L=../lib

DIR_PRIORS = $L/priors
DIR_CMAES  = $L/cmaes

LDLIBS += -L$(DIR_PRIORS) -lpriors
CFLAGS += -I$(DIR_PRIORS)

LDLIBS += -L$(DIR_CMAES) -lcmaes
CFLAGS += -I$(DIR_CMAES)

LDLIBS += -L$(DIR_FITFUN) -l$(LIB_FITFUN)
CFLAGS += -I$(DIR_FITFUN)

CXXFLAGS += `pkg-config libgp --cflags`
LDLIBS   += `pkg-config libgp --libs`

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_CMAES); make libcmaes.a)
	(cd $(DIR_FITFUN); make lib$(LIB_FITFUN).a)
	$(LINK.o.xx) $^ $(LDLIBS)

# dependencies

$(TARGET): engine_scmaes.o surrogate_samples.o surrogate.o $(FITFUN_OBJS)
engine_scmaes.o: engine_scmaes.c
surrogate.o: surrogate.cpp surrogate.h
surrogate_samples.o: surrogate_samples.c surrogate_samples.h

clean:; rm -rf  $(TARGET) *.o

test: scma
	atest `find test -name main`

.PHONY: clean all test
