# source/TMCMC/gaus/Makefile

include ../../make/common.mk

LIB_FITFUN   ?= fitfun_tests
LIB_FITFUNXX ?= fitfunxx
DIR_FITFUN   ?= ../../likelihoods
DIR_FITFUNXX ?= ../../likelihoodsxx
TARGET       ?= engine_tmcmc 

LEXT = ../../../lib/
DIR_STAN     ?= ${LEXT}
DIR_EIGEN    ?= ${LEXT}/stan/lib/eigen_3.3.3
DIR_BOOST    ?= ${LEXT}/stan/lib/boost_1.66.0
DIR_SUNDIALS ?= ${LEXT}/stan/lib/sundials_3.1.0/include

all: main_gaus_example \
	 main_manifold_gaus_example

L=../../lib
DIR_PRIORS   = $L/priors
DIR_PRIORSXX = $L/priorsxx
DIR_TMCMCXX  = $L/tmcmcxx

LDLIBS    := -L$(DIR_PRIORS) -lpriors $(LDLIBS)
CFLAGS    += -I$(DIR_PRIORS)

LDXLIBS  := -L$(DIR_PRIORSXX) -lpriorsxx $(LDXLIBS)
CXXFLAGS += -I$(DIR_PRIORSXX)

LDLIBS := -L$(DIR_FITFUN) -l$(LIB_FITFUN) $(LDLIBS)
CFLAGS += -I$(DIR_FITFUN)

LDXLIBS  := -L$(DIR_FITFUNXX) -l$(LIB_FITFUNXX) $(LDXLIBS)
CXXFLAGS += -I$(DIR_FITFUNXX)

LDXLIBS  := -L$(DIR_TMCMCXX) -ltmcmcxx $(LDXLIBS)
CXXFLAGS += -I$(DIR_TMCMCXX)

STANINCLUDES := -I$(DIR_STAN)
STANINCLUDES += -I$(DIR_EIGEN)
STANINCLUDES += -I$(DIR_BOOST)
STANINCLUDES += -I$(DIR_SUNDIALS)

CXXFLAGS += $(STANINCLUDES)

main_gaus_example:
	(cd $(DIR_PRIORSXX); make libpriorsxx.a)
	(cd $(DIR_FITFUNXX); make lib$(LIB_FITFUNXX).a)
	(cd $(DIR_TMCMCXX); make libtmcmcxx.a)
	$(LINK.o.xx) $^ $(LDXLIBS)

main_manifold_gaus_example:
	(cd $(DIR_PRIORSXX); make libpriorsxx.a)
	(cd $(DIR_FITFUNXX); make lib$(LIB_FITFUNXX).a)
	(cd $(DIR_TMCMCXX); make libtmcmcxx.a)
	$(LINK.o.xx) $^ $(LDXLIBS)


# dependencies

main_gaus_example: main_gaus_example.o
main_manifold_gaus_example: main_manifold_gaus_example.o

clean:; rm -rf main_gaus_example main_manifold_gaus_example  *.o
	(cd $(DIR_PRIORS); make clean)
	(cd $(DIR_PRIORSXX); make clean)
	(cd $(DIR_FITFUNXX); make clean)
	(cd $(DIR_TMCMCXX); make clean)

clean_run:; rm *.txt

.PHONY: clean all clean_run
