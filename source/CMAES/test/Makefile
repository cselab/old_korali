# source/CMAES/test/Makefile

include ../../make/common.mk

TARGET=testxx

all: $(TARGET)

L=../../lib

LIB_FITFUN ?= fitfun_tests

DIR_PRIORS  = $L/priors
DIR_FITFUN  = ../../likelihoods
DIR_CMAES   = $L/cmaes
DIR_CMAESXX = $L/cmaesxx
DIR_CATCH2  = ../../../lib/catch2

LDLIBS := -L$(DIR_PRIORS) -lpriors $(LDLIBS)
LDLIBS := -L$(DIR_CMAES) -lcmaes $(LDLIBS)
LDLIBS := -L$(DIR_FITFUN) -l$(LIB_FITFUN) $(LDLIBS)

LDXLIBS := $(LDLIBS) -L$(DIR_CMAESXX) -lcmaesxx

CXXFLAGS += -I$(DIR_CMAES)
CXXFLAGS += -I$(DIR_CMAESXX)
CXXFLAGS += -I$(DIR_PRIORS)
CXXFLAGS += -I$(DIR_FITFUN)
CXXFLAGS += -I$(DIR_CATCH2)

testxx.o: testxx.cpp

testxx: testxx.o
	$(LINK.o.xx) $< $(LDXLIBS)

clean:; rm -rf $(TARGET) *.o

.PHONY: clean all
