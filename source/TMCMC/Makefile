include ../make/common.mk

LIB_FITFUN   ?= fitfun_tests
DIR_FITFUN   ?= ../likelihoods
TARGET       ?= engine_tmcmc

DIR_STAN     ?= ../../lib/stan
DIR_EIGEN    ?= ../../lib/stan/lib/eigen_3.3.3
DIR_BOOST    ?= ../../lib/stan/lib/boost_1.66.0
DIR_SUNDIALS ?= ../../lib/stan/lib/sundials_3.1.0/include

all: $(TARGET)

L=../lib

DIR_PRIORS=$L/priors
LDLIBS := -L$(DIR_PRIORS) -lpriors $(LDLIBS)

CFLAGS += -I$(DIR_PRIORS)

LDLIBS := -L$(DIR_FITFUN) -l$(LIB_FITFUN) $(LDLIBS)
CFLAGS += -I$(DIR_FITFUN)

CFLAGS += -I$(DIR_STAN)
CFLAGS += -I$(DIR_EIGEN)
CFLAGS += -I$(DIR_BOOST)
CFLAGS += -I$(DIR_SUNDIALS)


$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_FITFUN); make lib$(LIB_FITFUN).a)
	$(LINK.o.xx) $^ $(LDLIBS)

# dependencies

$(TARGET): tmcmc_run.o \
	tmcmc_stats.o \
	tmcmc_db.o \
	tmcmc_aux.o \
	tmcmc_engine.o

tmcmc_run.o: tmcmc_run.c
tmcmc_stats.o: tmcmc_stats.c
tmcmc_db.o: tmcmc_db.c 
tmcmc_aux.o: tmcmc_aux.c
tmcmc_engine.o: tmcmc_engine.c

clean:; rm -rf  $(TARGET) *.o
	(cd $(DIR_PRIORS); make clean)

.PHONY: clean all
