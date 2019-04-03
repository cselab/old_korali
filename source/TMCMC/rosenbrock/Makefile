# source/TMCMC/rosenbrock/Makefile

include ../../make/common.mk

LIB_FITFUNXX ?= ../../fitfunxx
DIR_FITFUNXX ?= ../../likelihoodsxx

LEXT = ../../../lib/
DIR_STAN     ?= ${LEXT}
DIR_EIGEN    ?= ${LEXT}/stan/lib/eigen_3.3.3
DIR_BOOST    ?= ${LEXT}/stan/lib/boost_1.66.0
DIR_SUNDIALS ?= ${LEXT}/stan/lib/sundials_3.1.0/include

DIR_LIBGP    ?= ${LEXT}/libgp/libgp/
LIB_LIBGP    ?= gp
DIR_LIBGPLIB ?= $(DIR_LIBGP)build

all: main_rosenbrock_example \
     main_manifold_rosenbrock_example \
	 main_gp_rosenbrock_example \
	 main_gp_manifold_rosenbrock_example

L=../../lib
DIR_PRIORSXX = $L/priorsxx
DIR_TMCMCXX  = $L/tmcmcxx

LDXLIBS  := -L$(DIR_PRIORSXX) -lpriorsxx $(LDXLIBS)
CXXFLAGS += -I$(DIR_PRIORSXX)

LDXLIBS  := -L$(DIR_FITFUNXX) -l$(LIB_FITFUNXX) $(LDXLIBS)
CXXFLAGS += -I$(DIR_FITFUNXX)

LDXLIBS  := -L$(DIR_TMCMCXX) -ltmcmcxx $(LDXLIBS)
CXXFLAGS += -I$(DIR_TMCMCXX)

STANINCLUDES := -I$(DIR_STAN)
STANINCLUDES += -I$(DIR_EIGEN)
STANINCLUDES += -I$(DIR_BOOST)
STANINCLUDES += -I$(DIR_SUNDIALS)

CXXFLAGS += $(STANINCLUDES)

CXXFLAGS += -I$(DIR_LIBGP)
LDXLIBS  += -L$(DIR_LIBGPLIB) -l$(LIB_LIBGP)

main_rosenbrock_example:
	(cd $(DIR_PRIORSXX); make libpriorsxx.a)
	(cd $(DIR_FITFUNXX); make lib$(LIB_FITFUNXX).a)
	(cd $(DIR_TMCMCXX); make libtmcmcxx.a)
	$(LINK.o.xx) $^ $(LDXLIBS)

main_manifold_rosenbrock_example:
	(cd $(DIR_PRIORSXX); make libpriorsxx.a)
	(cd $(DIR_FITFUNXX); make lib$(LIB_FITFUNXX).a)
	(cd $(DIR_TMCMCXX); make libtmcmcxx.a)
	$(LINK.o.xx) $^ $(LDXLIBS)

main_gp_rosenbrock_example:
	(cd $(DIR_PRIORSXX); make libpriorsxx.a)
	(cd $(DIR_FITFUNXX); make lib$(LIB_FITFUNXX).a)
	(cd $(DIR_TMCMCXX); make libtmcmcxx.a)
	(cd $(DIR_LIBGPLIB); make)
	$(LINK.o.xx) $^ $(LDXLIBS)

main_gp_manifold_rosenbrock_example:
	(cd $(DIR_PRIORSXX); make libpriorsxx.a)
	(cd $(DIR_FITFUNXX); make lib$(LIB_FITFUNXX).a)
	(cd $(DIR_TMCMCXX); make libtmcmcxx.a)
	(cd $(DIR_LIBGPLIB); make)
	$(LINK.o.xx) $^ $(LDXLIBS)

# dependencies

main_rosenbrock_example: main_rosenbrock_example.o
main_manifold_rosenbrock_example: main_manifold_rosenbrock_example.o
main_gp_rosenbrock_example: main_gp_rosenbrock_example.o
main_gp_manifold_rosenbrock_example: main_gp_manifold_rosenbrock_example.o

main_rosenbrock_example.o: main_rosenbrock_example.cpp
main_manifold_rosenbrock_example.o: main_manifold_rosenbrock_example.cpp
main_gp_manifold_rosenbrock_example.o: main_gp_manifold_rosenbrock_example.cpp

clean:; rm -rf main_rosenbrock_example main_manifold_rosenbrock_example main_gp_rosenbrock_example main_gp_manifold_rosenbrock_example *.o
	(cd $(DIR_PRIORSXX); make clean)
	(cd $(DIR_FITFUNXX); make clean)
	(cd $(DIR_TMCMCXX); make clean)

clean_run:; rm *.txt

.PHONY: clean all clean_run
