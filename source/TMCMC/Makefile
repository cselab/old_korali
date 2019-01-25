# source/TMCMC/Makefile

include ../make/common.mk

LIB_FITFUN   ?= fitfun_tests
LIB_FITFUNXX ?= fitfunxx
DIR_FITFUN   ?= ../likelihoods
DIR_FITFUNXX ?= ../likelihoodsxx
TARGET       ?= engine_tmcmc

DIR_STAN     ?= ../../lib/stan
DIR_EIGEN    ?= ../../lib/stan/lib/eigen_3.3.3
DIR_BOOST    ?= ../../lib/stan/lib/boost_1.66.0
DIR_SUNDIALS ?= ../../lib/stan/lib/sundials_3.1.0/include

all: $(TARGET) tmcmcxx

L=../lib
DIR_PRIORS  = $L/priors
DIR_TMCMCXX = $L/tmcmcxx

LDLIBS    := -L$(DIR_PRIORS) -lpriors $(LDLIBS)
CFLAGS    += -I$(DIR_PRIORS)

LDXLIBS  := -L$(DIR_PRIORS) -lpriors $(LDXLIBS)
CXXFLAGS += -I$(DIR_PRIORS)

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

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_FITFUN); make lib$(LIB_FITFUN).a)
	(cd $(DIR_FITFUNXX); make lib$(LIB_FITFUNXX).a)
	(cd $(DIR_TMCMCXX); make libtmcmcxx.a)
	$(LINK.o) $^ $(LDLIBS)

tmcmcxx:
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_FITFUNXX); make lib$(LIB_FITFUNXX).a)
	(cd $(DIR_TMCMCXX); make libtmcmcxx.a)
	$(LINK.o.xx) $^ $(LDXLIBS)

# dependencies

$(TARGET): tmcmc_run.o \
	tmcmc_stats.o \
	tmcmc_db.o \
	tmcmc_aux.o \
	tmcmc_engine.o \

tmcmcxx: main.o

tmcmc_run.o: tmcmc_run.c
tmcmc_stats.o: tmcmc_stats.c
tmcmc_db.o: tmcmc_db.c 
tmcmc_aux.o: tmcmc_aux.c
tmcmc_engine.o: tmcmc_engine.c
main.o: main.cpp

clean:; rm -rf  $(TARGET) tmcmcxx *.o
	(cd $(DIR_PRIORS); make clean)
	(cd $(DIR_FITFUNXX); make clean)
	(cd $(DIR_TMCMCXX); make clean)

.PHONY: clean all
