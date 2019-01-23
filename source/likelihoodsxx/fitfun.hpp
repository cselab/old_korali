#ifndef FITFUN_HPP
#define FITFUN_HPP

typedef double (*fitfun_t)(const double*, int);

void fitfun_finalize();

double fitfun(const double *x, int N, void *output, int *info);

#endif //FITFUN_HPP
