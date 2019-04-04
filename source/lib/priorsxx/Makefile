# source/lib/priorsxx/Makefile

include ../../make/common.mk

TARGET := libpriorsxx.a

all: $(TARGET)

libpriorsxx.a: myrand.o mvgauss.o density.o priors.o 

myrand.o: myrand.cpp myrand.hpp

mvgauss.o: mvgauss.cpp

density.o: density.cpp density.hpp

priors.o: priors.cpp priors.hpp myrand.hpp density.hpp

test: $(TARGET)
	(cd test; make test)

clean:; rm -rf *.o *.a
	(cd test; make clean)

.PHONY: clean all

