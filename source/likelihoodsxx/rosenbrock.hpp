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

gsl_vector* gradX_minusRosenbrock(const double *x, int N)
{
    
    gsl_vector* grad = gsl_vector_calloc(N);

    double dxi, dxii;
    for (int i = 0; i < N-1; ++i) {
        dxi  = gsl_vector_get(grad, i) + 400*(x[i+1]-x[i]*x[i])*x[i] - 2*(x[i]-1);
        dxii = -200*(x[i+1]-x[i]*x[i]);
        gsl_vector_set(grad, i, dxi);
        gsl_vector_set(grad, i+1, dxii);
    }

    return grad;
}

gsl_matrix* hessX_minusRosenbrock(const double *x, int N)
{

    gsl_matrix* hess = gsl_matrix_calloc(N, N);
    gsl_matrix_set(hess, 0, 0, 400*x[1]-1200*x[0]*x[0]-2);

    gsl_matrix_set(hess, 0, 1, 400*x[0]);
    gsl_matrix_set(hess, 1, 0, 400*x[0]);

    double d2dxi2, d2dxidxj;
    for(int i = 1; i < N-1; ++i) {
        int j = i + 1;
        
        d2dxi2 = 400*x[j]-1200*x[i]*x[i]-202;
        gsl_matrix_set(hess, i, i, d2dxi2);

        d2dxidxj = 400*x[j];
        gsl_matrix_set(hess, i, j, d2dxidxj);
        gsl_matrix_set(hess, j, i, d2dxidxj);
    }

    gsl_matrix_set(hess, N-1, N-1, 400*x[N-2]);
    return hess;
}

    
}//namespace fitfun


#endif//ROSENBROCK_HPP
