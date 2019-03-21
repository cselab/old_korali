#ifndef ROSENBROCK_HPP
#define ROSENBROCK_HPP

#include <math.h>

#include "fitfun.hpp"

namespace fitfun
{

double f_minusRosenbrock(const double *x, int N)
{
    /* 
     * a = 1, b = 100
     *
     * in 1-D:
     *      f(x,y) = -(a-x)^2 - b(y-x^2)^2
     *      global max at: f(a,a^2) = 0
    */

    double s = 0.;
    for (int i = 0; i < N-1; ++i)
        s += 100.*pow(x[i+1]-x[i]*x[i], 2) + pow(x[i]-1., 2);
    return -s;
}

}//namespace fitfun


#endif//ROSENBROCK_HPP
