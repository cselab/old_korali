include ../make/common.mk

FITFUN_OBJS ?= ./fitfun.o
FITFUN_DIR  ?= ../likelihoods
TARGET      ?= engine_tmcmc

all: $(TARGET)

L=../lib

DIR_PRIORS=$L/priors
LDLIBS += -L$(DIR_PRIORS) -lpriors

CFLAGS += -I$(DIR_PRIORS)
CFLAGS += -I$(FITFUN_DIR)

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	$(LINK.o) $^ $(LDLIBS)

# dependencies

$(TARGET): tmcmc_run.o \
	tmcmc_stats.o \
	tmcmc_db.o \
	tmcmc_aux.o \
	tmcmc_engine.o \
	fitfun.o

tmcmc_run.o: tmcmc_run.c
tmcmc_stats.o: tmcmc_stats.c
tmcmc_db.o: tmcmc_db.c 
tmcmc_aux.o: tmcmc_aux.c
tmcmc_engine.o: tmcmc_engine.c

fitfun.o: fitfun.c

clean:; rm -rf  $(TARGET) *.o
clear:;	rm -rf .tags* tags *.txt *.dat actparcmaes.par *.err

.PHONY: clear clean all
