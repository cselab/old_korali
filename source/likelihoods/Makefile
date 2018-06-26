include ../make/common.mk

LIBS = \
	libftitfun_theta_int.a \
	libftitfun_theta_ext.a \
	libftitfun_psi.a       \
	libftitfun_posterior_theta_int.a \
	libftitfun_posterior_theta_ext.a

all: $(LIBS)

L=../lib
CFLAGS += -I$L/priors
# for <fitfun.h>
CFLAGS += -I.

OBJ_THETA_INT = fitfun_theta_fast.o loglike_theta_fast.o
OBJ_THETA_EXT = fitfun_theta.o      loglike_theta.o      spawner.o
OBJ_PSI       = fitfun_psi.o        loglike_psi.o

OBJ_POST = fitfun_posterior_theta.o loglike_posterior_theta.o

OBJ_POST_THETA_INT = $(OBJ_POST) loglike_theta_fast.o
OBJ_POST_THETA_EXT = $(OBJ_POST) loglike_theta.o         spawner.o


libftitfun_theta_int.a: $(OBJ_THETA_INT)
libftitfun_theta_ext.a: $(OBJ_THETA_EXT)

libftitfun_psi.a: $(OBJ_PSI)

libftitfun_posterior_theta_int.a: $(OBJ_POST_THETA_INT)
libftitfun_posterior_theta_ext.a: $(OBJ_POST_THETA_EXT)



clean:; rm -rf *.o $(LIBS)

.PHONY: clean
