#ifndef EXPERIMENTAL_FITFUN_H
#define EXPERIMENTAL_FITFUN_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <stan/math.hpp>

typedef struct {
    double loglike;
    int error_flg;
    int posdef;
    gsl_vector* grad;
    gsl_matrix* cov;
    gsl_matrix* evec;
    gsl_vector* eval;
} return_type;

class Fitfun
{

public:
    virtual return_type* fitfun (double* x, int n, void* output, int* info) = 0;

    virtual void fitfun_initialize(int argc, const  char **argv) = 0;

    virtual void fitfun_finalize() = 0;

};

#endif//EXPERIMENTAL_FITFUN_H
