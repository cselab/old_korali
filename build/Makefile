BUILD = $(shell pwd)
SRC   = $(BUILD)/../source
LIK   = $(SRC)/likelihoods
CMAES = $(SRC)/CMAES
TMCMC = $(SRC)/TMCMC
DRAM  = $(SRC)/DRAM

include $(SRC)/make/common.mk

CFLAGS += -I$(LIK)

TARGETS = \
	cmaes_theta_external           \
	cmaes_theta_internal           \
	tmcmc_theta_external           \
	tmcmc_theta_internal           \
	tmcmc_psi                      \
	tmcmc_posterior_theta_external \
	tmcmc_posterior_theta_internal \
	dram_theta_internal            \
	dram_psi

all: $(TARGETS)

link.cmaes = ( cd $(CMAES); make LIB_FITFUN="$<" TARGET="$(BUILD)/$@"; )
link.tmcmc = ( cd $(TMCMC); make LIB_FITFUN="$<" TARGET="$(BUILD)/$@"; )
link.dram  = ( cd $(DRAM);  make LIB_FITFUN="$<" TARGET="$(BUILD)/$@"; )

cmaes_theta_external:           fitfun_theta_ext;           $(link.cmaes)
cmaes_theta_internal:           fitfun_theta_int;           $(link.cmaes)
tmcmc_theta_external:           fitfun_theta_ext;           $(link.tmcmc)
tmcmc_theta_internal:           fitfun_theta_int;           $(link.tmcmc)
tmcmc_psi:                      fitfun_psi;                 $(link.tmcmc)
tmcmc_posterior_theta_external: fitfun_posterior_theta_ext; $(link.tmcmc)
tmcmc_posterior_theta_internal: fitfun_posterior_theta_int; $(link.tmcmc)
dram_theta_internal:            fitfun_theta_int;           $(link.dram)
dram_psi:                       fitfun_psi;                 $(link.dram)

fitfun%:; (cd $(LIK); make lib$@.a)

clean:; rm -rf $(TARGETS)

cleanall: clean
	(cd $(CMAES); make clean)
	(cd $(TMCMC); make clean)
	(cd $(DRAM);  make clean)
	(cd $(LIK);   make clean)

.PHONY: clean cleanall
