#include <stdlib.h>
#include <stdio.h>
#include <fitfun.h>
#include "loglike_psi.h"


void fitfun_initialize(int argc, char **argv) {

	loglike_psi_initialize();

}

void fitfun_finalize() {

	loglike_psi_finalize();

}

double fitfun(double *x, int n, void *output, int *info) {

	return loglike_psi( x, n, output,info );

}
