# source/lib/cames/Makefile

include ../../make/common.mk

TARGET=libcmaes.a

all: $(TARGET)

libcmaes.a: cmaes.o cmaes_utils.o
cmaes.o: cmaes.c cmaes.h cmaes_types.h
cmaes_utils.o: cmaes_utils.c cmaes_utils.h cmaes.h cmaes_types.h

clean:; rm -rf $(TARGET) *.o *.a

.PHONY: clean all
