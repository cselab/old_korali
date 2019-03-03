
CC = gcc
CPP = g++

CFLAGS 		=  
CPPFLAGS 	= -std=c++11 -g -O0

INCLUDE = -I../include/   -I../../eigen3/include/
LIBS = -L./../build -lgp -lm



CFLAGS += $(INCLUDE)
LFLAGS +=



TARGETS = testgp

all: $(TARGETS)


testgp:	testgp.o  gp_data.o
	$(CPP) $(LFLAGS) -o testgp gp_data.o testgp.o $(LIBS)


testgp.o: testgp.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c testgp.cpp


gp_data.o: gp_data.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c gp_data.cpp






clean:
	rm -f *.o

cleanall:
	rm -rf *.o $(TARGETS) *.txt data/grad* data/gp* data/out*
