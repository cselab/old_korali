# source/lib/camesxx/Makefile

include ../../make/common.mk

TARGET=libcmaesxx.a

all: $(TARGET)

L=..

DIR_CMAES  = $L/cmaes
DIR_PRIORS = $L/priors

CXXFLAGS += -I$(DIR_CMAES)
CXXFLAGS += -I$(DIR_PRIORS)

engine_cmaes.o: engine_cmaes.cpp engine_cmaes.hpp engine_cmaes_utils.hpp

libcmaesxx.a: engine_cmaes.o 

clean:; rm -rf $(TARGET) *.o *.a

.PHONY: clean all
