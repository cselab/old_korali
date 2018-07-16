include ../../make/common.mk

LIB = libpriors.a
TARGETS := $(LIB) test

all: $(TARGETS)
lib: $(LIB)

OBJS=priors.o myrand.o
libpriors.a: $(OBJS)
test: main.o $(OBJS)

priors.o: priors.c priors.h myrand.h
myrand.o: myrand.c myrand.h
main.o: main.c priors.h myrand.h
