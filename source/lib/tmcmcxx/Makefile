# source/lib/tmcmcxx/Makefile

include ../../make/common.mk

TARGET=libtmcmcxx.a

L=..

DIR_LIKELIHOOD = ../../likelihoodsxx
DIR_PRIORS = $L/priors

CXXFLAGS += -I$(DIR_PRIORS)
CXXFLAGS += -I$(DIR_LIKELIHOOD)

all: $(TARGET)

libtmcmcxx.a: engine_tmcmc.o tmcmc_stats.o tmcmc_utils.o tmcmc_types.o

engine_tmcmc.o: engine_tmcmc.cpp

tmcmc_stats.o: tmcmc_stats.cpp

tmcmc_utils.o: tmcmc_utils.cpp

tmcmc_types.o : tmcmc_types.cpp


clean:; rm -rf $(TARGET) *.o *.a

.PHONY: clean all
