# source/CMAES/Makefile

include ../make/common.mk

TARGET     ?= diff

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

# dependencies

$(TARGET): main.o  coupled_ode_system.o auto_diff_fitfun.o system_utils.o
	$(LINK.o.xx) $^ $(LDXLIBS)

main.o:	  main.cpp
coupled_ode_system.o: coupled_ode_system.cpp fitfun.hpp
auto_diff_fitfun.o: auto_diff_fitfun.cpp fitfun.hpp
system_utils.o: system_utils.cpp

clean:; rm -rf  $(TARGET) *.o *.a

.PHONY: clean
