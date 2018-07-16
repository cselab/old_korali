include ../make/common.mk

LIB_FITFUN ?= fitfun_tests
DIR_FITFUN ?= ../likelihoods
TARGET     ?= engine_tmcmc

all: $(TARGET)

L=../lib

DIR_PRIORS=$L/priors
LDLIBS += -L$(DIR_PRIORS) -lpriors

CFLAGS += -I$(DIR_PRIORS)

LDLIBS += -L$(DIR_FITFUN) -l$(LIB_FITFUN)
CFLAGS += -I$(DIR_FITFUN)


$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_FITFUN); make lib$(LIB_FITFUN).a)
	$(LINK.o) $^ $(LDLIBS)

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

.PHONY: clean all
