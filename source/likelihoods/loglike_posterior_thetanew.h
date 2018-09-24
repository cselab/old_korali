#ifndef _LOGLIKE_POSTERIOR_THETANEW_H_
#define _LOGLIKE_POSTERIOR_THETANEW_H_


void 	loglike_posterior_thetanew_initialize();
void 	loglike_posterior_thetanew_finalize();
double 	loglike_posterior_thetanew(double *x, int n, void *output, int *info );


#endif
