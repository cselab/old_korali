# source/CMAES/test/Makefile

include ../../make/common.mk

TARGET ?= testxx

all: $(TARGET)

L=../../lib

LIB_FITFUN ?= fitfun_tests

DIR_FITFUN  = ../../likelihoods
DIR_CATCH2  = ../../../lib/catch2
DIR_PRIORS  = $L/priors
DIR_CMAES   = $L/cmaes
DIR_CMAESXX = $L/cmaesxx

LDXLIBS += -L$(DIR_CMAESXX) -lcmaesxx  
LDXLIBS += -L$(DIR_CMAES) -lcmaes
LDXLIBS += -L$(DIR_PRIORS) -lpriors
LDXLIBS += -L$(DIR_FITFUN) -l$(LIB_FITFUN)

CXXFLAGS += -I$(DIR_CMAES)
CXXFLAGS += -I$(DIR_CMAESXX)
CXXFLAGS += -I$(DIR_PRIORS)
CXXFLAGS += -I$(DIR_FITFUN)
CXXFLAGS += -I$(DIR_CATCH2)

testxx.o: testxx.cpp

$(TARGET): testxx.o
	(cd ..; make cma)
	$(LINK.o.xx) testxx.o $(LDXLIBS)

test: $(TARGET)
	./$(TARGET)

clean:; rm -rf $(TARGET) *.o

.PHONY: clean all
