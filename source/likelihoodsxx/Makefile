include ../make/common.mk

TARGET=libfitfunxx.a

all: $(TARGET)

OBJ_TMP = fitfun_tmp.o

libfitfunxx.a: $(OBJ_TMP)


# deps

fitfun_tmp.o: fitfun_tmp.cpp

clean:; rm -rf *.o $(LIBS) *.a

.PHONY: clean
