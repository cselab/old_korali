include ../make/common.mk

FITFUN ?= ./fitfun
TARGETS := engine_tmcmc

all: $(TARGETS)

DIR_PRIORS=../priors
LDLIBS += -L$(DIR_PRIORS) -lpriors


engine_tmcmc:
	(cd $(DIR_PRIORS); make priors.a)
	$(LINK.o) $^ $(LDLIBS)

# dependencies

engine_tmcmc: tmcmc_run.o \
	tmcmc_stats.o \
	tmcmc_db.o \
	tmcmc_aux.o \
	tmcmc_engine.o \
	$(FITFUN).o

tmcmc_run.o: tmcmc_run.c
tmcmc_stats.o: tmcmc_stats.c
tmcmc_db.o: tmcmc_db.c 
tmcmc_aux.o: tmcmc_aux.c
tmcmc_engine.o: tmcmc_engine.c

$(FITFUN).o: $(FITFUN).c

clean:; rm -rf  $(TARGETS) *.o
clear:;	rm -rf .tags* tags *.txt *.dat actparcmaes.par *.err

.PHONY: clear clean all
