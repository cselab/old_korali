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



FITFUN_OBJ_EXT = $(LIK)/fitfun_theta_fast.o $(LIK)/loglike_theta_fast.o
FITFUN_OBJ_INT = $(LIK)/fitfun_theta.o      $(LIK)/loglike_theta.o      $(LIK)/spawner.o
FITFUN_OBJ_PSI = $(LIK)/fitfun_psi.o        $(LIK)/loglike_psi.o

OBJ_POST = $(LIK)/fitfun_posterior_theta.o $(LIK)/loglike_posterior_theta.o

FITFUN_OBJ_POST_INT = $(OBJ_POSTERIOR) $(LIK)/loglike_theta_fast.o
FITFUN_OBJ_POST_EXT = $(OBJ_POSTERIOR) $(LIK)/loglike_theta.o      $(LIK)/spawner.o

link.cmaes = ( cd $(CMAES); make FITFUN_OBJS="$^" TARGET="$@"; )
link.tmcmc = ( cd $(TMCMC); make FITFUN_OBJS="$^" TARGET="$@"; )
link.dram  = ( cd $(DRAM);  make FITFUN_OBJS="$^" TARGET="$@"; )

$(BUILD)/cmaes_theta_external:           $(FITFUN_OBJ_EXT);      $(link.cmaes)
$(BUILD)/cmaes_theta_internal:           $(FITFUN_OBJ_INT);      $(link.cmaes)
$(BUILD)/tmcmc_theta_external:           $(FITFUN_OBJ_EXT);      $(link.tmcmc)
$(BUILD)/tmcmc_theta_internal:           $(FITFUN_OBJ_INT);      $(link.tmcmc)
$(BUILD)/tmcmc_psi:                      $(FITFUN_OBJ_PSI);      $(link.tmcmc)
$(BUILD)/tmcmc_posterior_theta_external: $(FITFUN_OBJ_POST_EXT); $(link.tmcmc)
$(BUILD)/tmcmc_posterior_theta_internal: $(FITFUN_OBJ_POST_INT); $(link.tmcmc)
$(BUILD)/dram_theta_internal:            $(FITFUN_OBJ_INT);      $(link.dram)
$(BUILD)/dram_psi:                       $(FITFUN_OBJ_PSI);      $(link.dram)

clean:; rm -rf $(TARGETS)

cleanall: clean
	(cd $(CMAES); make clean)
	(cd $(TMCMC); make clean)
	(cd $(DRAM);  make clean)

.PHONY: clean cleanall
