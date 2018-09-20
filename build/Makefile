BUILD = $(shell pwd)
SRC   = $(BUILD)/../source
LIK   = $(SRC)/likelihoods
CMAES = $(SRC)/CMAES
TMCMC = $(SRC)/TMCMC
DRAM  = $(SRC)/DRAM

include $(SRC)/make/common.mk

CFLAGS += -I$(LIK)

TARGETS = \
	$(BUILD)/cmaes_theta_external           \
	$(BUILD)/cmaes_theta_internal           \
	$(BUILD)/tmcmc_theta_external           \
	$(BUILD)/tmcmc_theta_internal           \
	$(BUILD)/tmcmc_psi                      \
	$(BUILD)/tmcmc_posterior_theta_external \
	$(BUILD)/tmcmc_posterior_theta_internal \
	$(BUILD)/dram_theta_internal            \
	$(BUILD)/dram_psi

all: $(TARGETS)

link.cmaes = ( cd $(CMAES); make LIB_FITFUN="$<" TARGET="$@"; )
link.tmcmc = ( cd $(TMCMC); make LIB_FITFUN="$<" TARGET="$@"; )
link.dram  = ( cd $(DRAM);  make LIB_FITFUN="$<" TARGET="$@"; )

$(BUILD)/cmaes_theta_external:           fitfun_theta_ext;           $(link.cmaes)
$(BUILD)/cmaes_theta_internal:           fitfun_theta_int;           $(link.cmaes)
$(BUILD)/tmcmc_theta_external:           fitfun_theta_ext;           $(link.tmcmc)
$(BUILD)/tmcmc_theta_internal:           fitfun_theta_int;           $(link.tmcmc)
$(BUILD)/tmcmc_psi:                      fitfun_psi;                 $(link.tmcmc)
$(BUILD)/tmcmc_posterior_theta_external: fitfun_posterior_theta_ext; $(link.tmcmc)
$(BUILD)/tmcmc_posterior_theta_internal: fitfun_posterior_theta_int; $(link.tmcmc)
$(BUILD)/dram_theta_internal:            fitfun_theta_int;           $(link.dram)
$(BUILD)/dram_psi:                       fitfun_psi;                 $(link.dram)

fitfun%:; (cd $(LIK); make lib$@.a)

clean:; rm -rf $(TARGETS)

cleanall: clean
	(cd $(CMAES); make clean)
	(cd $(TMCMC); make clean)
	(cd $(DRAM);  make clean)
	(cd $(LIK);   make clean)

.PHONY: clean cleanall
