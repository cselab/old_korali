# source/TMCMC/blood/Makefile

include ../../make/common.mk

all: tmcmc_gp_blood \
#	 mtmcmc_gp_blood

LEXT = ../../../lib
DIR_STAN     ?= ${LEXT}
DIR_EIGEN    ?= ${LEXT}/stan/lib/eigen_3.3.3
DIR_BOOST    ?= ${LEXT}/stan/lib/boost_1.66.0
DIR_SUNDIALS ?= ${LEXT}/stan/lib/sundials_3.1.0/include

DIR_LIBGP    ?= ${LEXT}/libgp/libgp/
LIB_LIBGP    ?= gp
DIR_LIBGPLIB ?= $(DIR_LIBGP)build

L=../../lib
DIR_PRIORS   = $L/priors

LDXLIBS  := -L$(DIR_PRIORS) -lpriors $(LDXLIBS)
CXXFLAGS += -I$(DIR_PRIORS) -I$(DIR_PRIORS)

DIR_PRIORSXX = $L/priorsxx

LDXLIBS  := -L$(DIR_PRIORSXX) -lpriorsxx $(LDXLIBS)
CXXFLAGS += -I$(DIR_PRIORSXX) -I$(DIR_PRIORS)

LIB_FITFUNXX = fitfunxx
DIR_FITFUNXX = ../../likelihoodsxx

LDXLIBS  := -L$(DIR_FITFUNXX) -l$(LIB_FITFUNXX) $(LDXLIBS)
CXXFLAGS += -I$(DIR_FITFUNXX)

DIR_TMCMCXX  = $L/tmcmcxx

LDXLIBS  := -L$(DIR_TMCMCXX) -ltmcmcxx $(LDXLIBS)
CXXFLAGS += -I$(DIR_TMCMCXX)


STANINCLUDES := -I$(DIR_STAN)
STANINCLUDES += -I$(DIR_EIGEN)
STANINCLUDES += -I$(DIR_BOOST)
STANINCLUDES += -I$(DIR_SUNDIALS)

CXXFLAGS += $(STANINCLUDES)

CXXFLAGS += -I$(DIR_LIBGP)
LDXLIBS  += -L$(DIR_LIBGPLIB) -l$(LIB_LIBGP)

tmcmc_gp_blood: tmcmc_gp_blood.o	
	$(LINK.o.xx) $^ $(LDXLIBS)

tmcmc_gp_blood.o: tmcmc_gp_blood.cpp

mtmcmc_gp_blood: mtmcmc_gp_blood.o	
	$(LINK.o.xx) $^ $(LDXLIBS)

mtmcmc_gp_blood.o: mtmcmc_gp_blood.cpp

clean:; rm -rf tmcmc_gp_blood mtmcmc_gp_blood *.o

.PHONY: clean all
