# source/CMAES/blood/Makefile

# run: make DATA=XXX PARAM=YYY
#CXXFLAGS += -DDATASET$(DATA)
#CXXFLAGS += -DPARAM$(PARAM)

CXXFLAGS += -DDATASET$(DATA)

include ../../make/common.mk

all: cmaes_x_blood \
	 cmaes_gp_blood


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

DIR_CMAES    = $L/cmaes
DIR_CMAESXX  = $L/cmaesxx

CXXFLAGS += -I$(DIR_CMAES)
CXXFLAGS += -I$(DIR_CMAESXX)
LDXLIBS  := -L$(DIR_CMAES) -lcmaes $(LDXLIBS)
LDXLIBS  := -L$(DIR_CMAESXX) -lcmaesxx $(LDXLIBS)

STANINCLUDES := -I$(DIR_STAN)
STANINCLUDES += -I$(DIR_EIGEN)
STANINCLUDES += -I$(DIR_BOOST)
STANINCLUDES += -I$(DIR_SUNDIALS)

CXXFLAGS += $(STANINCLUDES)

CXXFLAGS += -I$(DIR_LIBGP)
LDXLIBS  += -L$(DIR_LIBGPLIB) -l$(LIB_LIBGP)


cmaes_x_blood: cmaes_x_blood.o	
	$(LINK.o.xx) $^ $(LDXLIBS)

cmaes_gp_blood: cmaes_gp_blood.o	
	$(LINK.o.xx) $^ $(LDXLIBS)

cmaes_x_blood.o: cmaes_x_blood.cpp blood.hpp
cmaes_gp_blood.o: cmaes_gp_blood.cpp blood.hpp

clean:; rm -rf cmaes_gp_blood cmaes_x_blood *.o

.PHONY: clean all
