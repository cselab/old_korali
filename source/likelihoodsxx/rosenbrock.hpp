#ifndef ROSENBROCK_HPP
#define ROSENBROCK_HPP

#include <math.h>

#include "fitfun.hpp"

namespace fitfun
{

double f_minusRosenbrock(const double *x, int N)
{
    int i;
    double s = 0.;
    for (i = 0; i < N-1; ++i)
        s += 100.*pow(x[i+1]-x[i]*x[i], 2) + pow(x[i]-1., 2);
    return -s;
}

}//namespace fitfun


#endif//ROSENBROCK_HPP
