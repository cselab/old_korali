include ../make/common.mk

LIB_FITFUN ?= fitfun_tests
DIR_FITFUN ?= ../likelihoods
TARGET     ?= dram

all: $(TARGET)

L=../lib

DIR_PRIORS=$L/priors

CFLAGS += -I$(DIR_PRIORS)
LDLIBS += -L$(DIR_PRIORS) -lpriors

LDLIBS += -L$(DIR_FITFUN) -l$(LIB_FITFUN)
CFLAGS += -I$(DIR_FITFUN)


$(TARGET):
	(cd $(DIR_PRIORS); make libpriors.a)
	(cd $(DIR_FITFUN); make lib$(LIB_FITFUN).a)
	$(LINK.o) $^ $(LDLIBS)

# dependencies

$(TARGET): dram.o run_dram.o
dram.o: dram.c dram.h
run_dram.o: run_dram.c dram.h

clean:; rm -rf  $(TARGET) *.o

.PHONY: clear clean all
