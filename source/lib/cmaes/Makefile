include ../../make/common.mk

TARGET=libcmaes.a

all: $(TARGET)

libcmaes.a: cmaes.o
cmaes.o: cmaes.c cmaes.h cmaes_types.h

clean:; rm -rf $(TARGET) *.o

.PHONY: clean all
