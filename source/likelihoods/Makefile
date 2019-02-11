include ../make/common.mk

LIBS = \
	libfitfun_theta_int.a \
	libfitfun_theta_ext.a \
	libfitfun_psi.a       \
	libfitfun_posterior_theta_int.a \
	libfitfun_posterior_theta_ext.a \
	libfitfun_posterior_theta_extend_int.a \
	libfitfun_posterior_theta_extend_ext.a \
	libfitfun_tests.a

all: $(LIBS)

L=../lib
CFLAGS += -I$L/priors

CFLAGS += -I.

OBJ_THETA_INT = fitfun_theta_fast.o loglike_theta_fast.o
OBJ_THETA_EXT = fitfun_theta.o      loglike_theta.o      spawner.o
OBJ_PSI       = fitfun_psi.o        loglike_psi.o

OBJ_POST = fitfun_posterior_theta.o loglike_posterior_theta.o
OBJ_POST_THETA_INT = $(OBJ_POST) loglike_theta_fast.o
OBJ_POST_THETA_EXT = $(OBJ_POST) loglike_theta.o         spawner.o

OBJ_POST_EXTND = fitfun_posterior_theta.o loglike_posterior_theta_extended.o
OBJ_POST_THETA_EXTND_INT = $(OBJ_POST_EXTND) loglike_theta_fast.o
OBJ_POST_THETA_EXTND_EXT = $(OBJ_POST_EXTND) loglike_theta.o         spawner.o

OBJ_TESTS = fitfun_tests.o


libfitfun_theta_int.a: $(OBJ_THETA_INT)
libfitfun_theta_ext.a: $(OBJ_THETA_EXT)

libfitfun_psi.a: $(OBJ_PSI)

libfitfun_posterior_theta_int.a: $(OBJ_POST_THETA_INT)
libfitfun_posterior_theta_ext.a: $(OBJ_POST_THETA_EXT)

libfitfun_posterior_theta_extend_int.a: $(OBJ_POST_THETA_EXTND_INT)
libfitfun_posterior_theta_extend_ext.a: $(OBJ_POST_THETA_EXTND_EXT)

libfitfun_tests.a: $(OBJ_TESTS)


# deps
fitfun_posterior_theta.o: fitfun_posterior_theta.c fitfun.h loglike_posterior_theta.h
fitfun_psi.o: fitfun_psi.c fitfun.h loglike_psi.h
fitfun_tests.o: fitfun_tests.c fitfun_tests.h fitfun.h
fitfun_theta.o: fitfun_theta.c fitfun.h loglike_theta.h
fitfun_theta_fast.o: fitfun_theta_fast.c fitfun.h loglike_theta_fast.h
loglike_posterior_theta_extended.o: loglike_posterior_theta_extended.c loglike_posterior_theta.h fitfun.h
loglike_posterior_theta.o: loglike_posterior_theta.c loglike_posterior_theta.h fitfun.h
loglike_psi.o: loglike_psi.c loglike_psi.h
loglike_theta.o: loglike_theta.c loglike_theta.h fitfun.h
loglike_theta_fast.o: loglike_theta_fast.c loglike_theta_fast.h fitfun.h
spawner.o: spawner.c spawner.h

clean:; rm -rf *.o $(LIBS) *.a

.PHONY: clean
