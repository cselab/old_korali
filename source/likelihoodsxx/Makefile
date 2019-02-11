include ../make/common.mk

TARGET=libfitfunxx.a

all: $(TARGET)

ELIB = ../../lib

DIR_STAN     ?= $(ELIB)
DIR_EIGEN    ?= $(ELIB)/stan/lib/eigen_3.3.3
DIR_BOOST    ?= $(ELIB)/stan/lib/boost_1.66.0
DIR_SUNDIALS ?= $(ELIB)/stan/lib/sundials_3.1.0/include

STANINCLUDES := -I$(DIR_STAN)
STANINCLUDES += -I$(DIR_EIGEN)
STANINCLUDES += -I$(DIR_BOOST)
STANINCLUDES += -I$(DIR_SUNDIALS)

CXXFLAGS += $(STANINCLUDES) -Wfatal-errors

libfitfunxx.a: system_utils.o coupled_ode_system.o auto_diff_fitfun.o

# deps
system_utils.o: system_utils.cpp

coupled_ode_system.o: coupled_ode_system.cpp

auto_diff_fitfun.o: auto_diff_fitfun.cpp

clean:; rm -rf *.o *.a

.PHONY: clean
