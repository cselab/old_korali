#ifndef EXPERIMENTAL_FITFUN_H
#define EXPERIMENTAL_FITFUN_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <stan/math.hpp>

class Fitfun
{

public:
    virtual double fitfun (double* x, int n, void* output, int* info) = 0;

    virtual void fitfun_initialize(int argc, const  char **argv) = 0;

    virtual void fitfun_finalize() = 0;

};

#endif//EXPERIMENTAL_FITFUN_H
