# source/lib/priorsxx/Makefile

include ../../make/common.mk

TARGET := libpriorsxx.a

all: $(TARGET)


libpriorsxx.a: myrand.o density.o priors.o

myrand.o: myrand.cpp myrand.hpp

density.o: density.cpp density.hpp

priors.o: priors.cpp priors.hpp myrand.hpp density.hpp



clean:; rm -rf *.o *.a

.PHONY: clean all

