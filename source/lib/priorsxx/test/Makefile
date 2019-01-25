# source/lib/priors/test/Makefile

include ../../../make/common.mk

TARGET ?= testxx

all: $(TARGET)

L=../../

LIB_PRIORSXX = priorsxx
DIR_CATCH2   = ../../../../lib/catch2
DIR_PRIORS   = ../

LDXLIBS := -L$(DIR_PRIORS) -lpriorsxx $(LDXLIBS)

CXXFLAGS += -I$(DIR_CATCH2)
CXXFLAGS += -I$(DIR_PRIORS)

$(TARGET): testxx.o
	(cd $(DIR_PRIORS); make libpriorsxx.a) 
	$(LINK.o.xx) $^ $(LDXLIBS)

testxx.o: testxx.cpp

test: $(TARGET)

ifeq ($(use_torc),1)
	$(warning Testing while Torc activated  may result in unexpected behavior)
endif
	./$(TARGET)

clean:; rm -rf $(TARGET) *.o

.PHONY: clean all
