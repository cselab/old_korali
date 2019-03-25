#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP

#include <math.h>

//#include "fitfun.hpp"

namespace fitfun
{

double f_gaus(const double *x, int N)
{
    double s = 0.;
    for (int i = 0; i < N; ++i) s += x[i]*x[i];
    return pow(sqrt(2*M_PI),-0.5*N)*exp(-0.5*s);
}

gsl_vector* gradX_gaus(const double *x, int N)
{
   
    double f = f_gaus(x,N);
    gsl_vector* grad = gsl_vector_calloc(N);

    for (int i = 0; i < N; ++i) gsl_vector_set(grad,i,-x[i]*f);

    return grad;
}

gsl_matrix* hessX_gaus(const double *x, int N)
{

    double f = f_gaus(x,N);
    gsl_matrix* hess = gsl_matrix_calloc(N, N);

    for(int i = 0; i < N; ++i) gsl_matrix_set(hess, i, i, (x[i]*x[i]-1)*f);


    return hess;
}

    
}//namespace fitfun

#endif//GAUSSIAN_HPP
