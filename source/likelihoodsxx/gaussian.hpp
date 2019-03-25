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

    /*
    printf("N:%d\n",N);
    printf("s:%f\n",s);
    printf("x0:%f\n",x[0]);
    printf("x1:%f\n",x[1]);
    */

    return pow(sqrt(2*M_PI),-0.5*N)*exp(-0.5*s);
}

gsl_vector* gradX_gaus(const double *x, int N)
{
   
    double f = f_gaus(x,N);
    gsl_vector* grad = gsl_vector_calloc(N);

    for (int i = 0; i < N; ++i) gsl_vector_set(grad,i,-x[i]*f);

    //printf("grad0:%f\n", gsl_vector_get(grad,0));
    //printf("grad1:%f\n", gsl_vector_get(grad,1));
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
