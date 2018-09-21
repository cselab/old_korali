include ../make/common.mk

LIB_FITFUN ?= fitfun_tests
DIR_FITFUN ?= ../likelihoods
TARGET     ?= cma

all: $(TARGET)

L=../lib

DIR_PRIORS = $L/priors
DIR_CMAES  = $L/cmaes

LDLIBS := -L$(DIR_PRIORS) -lpriors $(LDLIBS)
CFLAGS += -I$(DIR_PRIORS)

LDLIBS := -L$(DIR_CMAES) -lcmaes $(LDLIBS)
CFLAGS += -I$(DIR_CMAES)

LDLIBS := -L$(DIR_FITFUN) -l$(LIB_FITFUN) $(LDLIBS)
CFLAGS += -I$(DIR_FITFUN)

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_CMAES);  make libcmaes.a)
	(cd $(DIR_FITFUN); make lib$(LIB_FITFUN).a)
	$(LINK.o) $^ $(LDLIBS)

# dependencies

$(TARGET): engine_cmaes.o
engine_cmaes.o: engine_cmaes.c

clean:; rm -rf  $(TARGET) *.o *.a
	(cd $(DIR_PRIORS); make clean)
	(cd $(DIR_CMAES); make clean)

test: cma
	atest `find test -name main`

.PHONY: clean all test
