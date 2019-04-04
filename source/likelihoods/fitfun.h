#ifndef _FITFUN_H_
#define _FITFUN_H_

typedef double (*fitfun_t)(double*, int);

void fitfun_finalize();

double fitfun(double *x, int N, void *output, int *info);

#endif
