
CC = gcc
CPP = g++

CFLAGS 		=  
CPPFLAGS 	=

INCLUDE = -I../include/   -I../../eigen3/include/
LIBS = -L./../build -lgp



CFLAGS += $(INCLUDE)
LFLAGS +=



TARGETS = testgp_reg

all: $(TARGETS)


testgp_reg:	testgp_reg.o
	$(CPP) $(LFLAGS) -o testgp_reg testgp_reg.o $(LIBS)


testgp_reg.o: testgp_reg.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c testgp_reg.cpp


clean:
	rm -f *.o

clear:
	rm -rf *.o $(TARGETS) *.txt
