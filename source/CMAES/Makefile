# source/CMAES/Makefile

include ../make/common.mk

LIB_FITFUN ?= fitfun_tests
DIR_FITFUN ?= ../likelihoods
TARGET     ?= cma

all: $(TARGET)

L=../lib

DIR_PRIORS   = $L/priors
DIR_CMAES    = $L/cmaes
DIR_CMAESXX  = $L/cmaesxx
DIR_TESTS    = ./test/

LDLIBS   := -L$(DIR_PRIORS) -lpriors $(LDLIBS)
CFLAGS   += -I$(DIR_PRIORS)

LDLIBS   := -L$(DIR_CMAES) -lcmaes $(LDLIBS)
CFLAGS   += -I$(DIR_CMAES)

LDLIBS   := -L$(DIR_FITFUN) -l$(LIB_FITFUN) $(LDLIBS)
CFLAGS   += -I$(DIR_FITFUN)

LDXLIBS := $(LDLIBS) -L$(DIR_CMAESXX) -lcmaesxx

$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_CMAES);  make libcmaes.a)
	(cd $(DIR_CMAESXX);  make libcmaesxx.a)
	(cd $(DIR_FITFUN); make lib$(LIB_FITFUN).a)
	$(LINK.o) engine_cmaes.o $(LDLIBS)

# dependencies

$(TARGET): engine_cmaes.o
engine_cmaes.o: engine_cmaes.c 

clean:; rm -rf  $(TARGET) *.o *.a
	(cd $(DIR_PRIORS); make clean)
	(cd $(DIR_CMAES); make clean)
	(cd $(DIR_CMAESXX); make clean)
	(cd $(DIR_TESTS); make clean)

test: cma
	atest `find test -name main`

testxx: cma
	(cd $(DIR_TESTS); make test)

.PHONY: clean all test testxx
