# source/lib/tmcmcxx/Makefile

include ../../make/common.mk

TARGET=libtmcmcxx.a

L=..


DIR_LIKELIHOOD = ../../likelihoodsxx
DIR_PRIORS = $L/priorsxx

LIB_FITFUNXX = fitfunxx

CXXFLAGS += -I$(DIR_PRIORS)
CXXFLAGS += -I$(DIR_LIKELIHOOD)
LDXLIBS  := -L$(DIR_LIKELIHOOD) -l$(LIB_FITFUNXX) $(LDXLIBS)

# tmcmc objective minimizers (can be combined)

CXXFLAGS += -D_USE_FMINCON_
#CXXFLAGS += -D_USE_FMINSEARCH_ 
#CXXFLAGS += -D_USE_FZEROFIND_

CXXFLAGS += -DCHECK_POSDEF

all: $(TARGET)

engine_tmcmc.o: engine_tmcmc.cpp engine_tmcmc.hpp

tmcmc_stats.o: tmcmc_obj_fmin.cpp tmcmc_obj_fmin.hpp

tmcmc_utils.o: tmcmc_utils.cpp tmcmc_utils.hpp

tmcmc_types.o: tmcmc_types.cpp tmcmc_types.hpp

libtmcmcxx.a: engine_tmcmc.o tmcmc_obj_fmin.o tmcmc_utils.o tmcmc_types.o

clean:; rm -rf $(TARGET) *.o *.a

.PHONY: clean all
