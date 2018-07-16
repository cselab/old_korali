#include <stdlib.h>
#include <stdio.h>
#include <fitfun.h>
#include "loglike_psi.h"


void fitfun_initialize(int argc, char **argv) {
	int n;        
        char *s;
        if (argc != 1) {
            fprintf(stderr, "expected %d arguments, got %d\n", 1, argc);
            exit(1);
        }
        s = *argv;
	n = atoi(s);
	loglike_psi_initialize(n);
}

void fitfun_finalize() {

	loglike_psi_finalize();

}

double fitfun(double *x, int n, void *output, int *info) {

	return loglike_psi( x, n, output,info );

}
