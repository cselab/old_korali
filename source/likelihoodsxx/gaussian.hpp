#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP

#include <math.h>

namespace fitfun
{

double f_gaus(const double *x, int N)
{
    double s = 0.;
    for (int i = 0; i < N; ++i) s += x[i]*x[i];
    return log(pow(2*M_PI,-0.5*N)) - 0.5*s;
}


gsl_vector* gradX_gaus(const double *x, int N)
{
   
    gsl_vector* grad = gsl_vector_calloc(N);

    for (int i = 0; i < N; ++i) gsl_vector_set(grad,i,-x[i]);

    return grad;
}

gsl_matrix* hessX_gaus(const double *x, int N)
{
    double f = f_gaus(x,N);
    
    gsl_matrix* hess = gsl_matrix_calloc(N, N);

    for(int i = 0; i < N; ++i) gsl_matrix_set(hess, i, i, -1.0);

    return hess;
}

    
}//namespace fitfun

#endif//GAUSSIAN_HPP
