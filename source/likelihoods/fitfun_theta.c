#include <fitfun.h>
#include "loglike_theta.h"


void fitfun_initialize(int argc, char **argv) {

}



void fitfun_finalize() {

}




double fitfun(double *x, int n, void *output, int *info) {

	return loglike_theta( x, n, output, info );

}
