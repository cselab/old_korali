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

LDXLIBS := -L$(DIR_PRIORS) -lpriors $(LDXLIBS)
LDXLIBS := -L$(DIR_CMAES) -lcmaes $(LDXLIBS)
LDXLIBS := -L$(DIR_CMAESXX) -lcmaesxx $(LDXLIBS)  
LDXLIBS := -L$(DIR_FITFUN) -l$(LIB_FITFUN) $(LDXLIBS)

CXXFLAGS += -I$(DIR_CMAES)
CXXFLAGS += -I$(DIR_CMAESXX)
CXXFLAGS += -I$(DIR_PRIORS)
CXXFLAGS += -I$(DIR_FITFUN)
CXXFLAGS += -I$(DIR_CATCH2)

$(TARGET): testxx.o
	$(LINK.o.xx) $^ $(LDXLIBS)

testxx.o: testxx.cpp

test: $(TARGET)
	./$(TARGET)

clean:; rm -rf $(TARGET) *.o

.PHONY: clean all
