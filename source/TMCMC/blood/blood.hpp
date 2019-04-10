#ifndef BLOOD_HPP
#define BLOOD_HPP

#include <stdio.h>
#include <algorithm>
#include <functional>

#include <libgp/libgp/include/core/gp.h>
#include <libgp/libgp/include/core/gp_utils.h>
#include <libgp/libgp/include/core/rprop.h>

#include "fitfun.hpp"
#include "engine_tmcmc.hpp"

using namespace libgp;
using namespace fitfun;
using namespace tmcmc;

//size_t NEXP  = 5; // num experimental data
//double F[] = { 5.714981, 11.435213, 12.377169, 22.796679, 24.760591 }; //3
//double OUT[] = { 0.029267, 0.026932, 0.023028, 0.021188, 0.021306 }; //3

//size_t NEXP  = 3; // num experimental data
//double F[] = { 6.059418, 10.820390, 20.602023 }; //4
//double OUT[] = { 0.035286, 0.032080, 0.030420 }; //4

size_t NEXP  = 5; // num experimental data
double F[]   = { 6.555291, 9.832936, 13.110581, 16.388226, 19.665872 }; // 5: experimental data (Force)
double OUT[] = { 0.039364, 0.030558, 0.032415, 0.028117, 0.028064 };    // 5: experimental data (TTf)

//size_t NEXP  = 5; // num experimental data
//double F[] = { 5.867770, 11.611925, 12.005029, 23.673504, 24.671042 }; //6
//double OUT[] = { 0.035660, 0.031749, 0.031790, 0.028413, 0.028243 }; //6

//size_t NEXP  = 2; // num experimental data
//double F[] = { 6.242403, 12.590044 }; //7
//double OUT[] = { 0.031000, 0.028000 }; //7

//GaussianProcess gpblood("gp_blood_cmaes_opt.txt");
GaussianProcess gpblood("gp_blood.txt");

// theta: Q1, Q2, Q3, Q4, mu0, sigma
// N    : 6
double gpllk(const double* theta, int N) {

    double var = theta[N-1]*theta[N-1];

    double x_gp[N-1]; // Q1, Q2, Q3, Q4, F'
    for(int i = 0; i < N-2; ++i) x_gp[i] = theta[i];
    //std::copy(theta,theta+N-2,x_gp);

    double sse = 0.0;
    for(size_t i = 0; i < NEXP; ++i) {
        // Athena's transformation
        x_gp[N-2] = (F[i]/theta[N-2] - 4.0) / 16.0; 
        sse+= pow((gpblood.f(x_gp) - OUT[i]),2);
    }

    return -0.5*NEXP*log(2*M_PI*var)-0.5*sse/var;
}

// theta: Q1, Q2, Q3, Q4, mu0, sigma
// N    : 6
gsl_vector * gpllk_grad(const double * theta, int N) {

    gsl_vector * grad = gsl_vector_calloc(N);

    double var = theta[N-1]*theta[N-1];
    double mu216 = 16.0*theta[N-2]*theta[N-2];
    
    double x_gp[N-1];
    std::copy(theta,theta+N-3,x_gp);

    double sse = 0.0;
    double diff;
    gsl_vector * tmp = gsl_vector_calloc(N);
    for(size_t i = 0; i < NEXP; ++i) {
        // Athena's transformation
        x_gp[N-2] = (F[i]/theta[N-2] - 4.0) / 16.0; 
        
        Eigen::VectorXd gpgradx = gpblood.dfdx(x_gp);
        gpgradx[N-2] *= -1.0/mu216; //derivative wrt A's transformation

        diff = gpblood.f(x_gp) - OUT[i];
        sse += diff*diff;
        for(int k = 0; k < N-1; ++k) gsl_vector_set(tmp, k, diff*gpgradx[k]);
        gsl_vector_add(grad, tmp);
    }
    gsl_vector_free(tmp);

    gsl_vector_set(grad, N-1, -1.0*NEXP*theta[N-1]+sse/theta[N-1]); //divide by var below
    gsl_vector_scale(grad,1.0/var);                                 //divide all elements by var

#ifdef DEBUG
    printf("_sse: %f\n", sse);
    printf("_sig: %f\n", theta[N-1]);
    printf("_grad:");  
    for(int i = 0; i <  N; ++i) printf("  %f  ", *(grad->data+i));
    printf("\n");
#endif

    return grad;
}

// theta: Q1, Q2, Q3, Q4, mu0, sigma
// N    : 6
gsl_matrix * gpllk_FIM(const double * theta, int N) {
    
    gsl_matrix * FIM = gsl_matrix_calloc(N,N);

    double var   = theta[N-1]*theta[N-1];
    double mu216 = 16.0*theta[N-2]*theta[N-2];
    
    double x_gp[N-1];
    std::copy(theta,theta+N-2,x_gp);

    gsl_matrix * tmp = gsl_matrix_calloc(N,N);
    for(size_t i = 0; i < NEXP; ++i) {
        // Athena's transformation
        x_gp[N-1] = (F[i]/theta[N-2] - 4.0) / 16.0;

        Eigen::VectorXd gpgradx = gpblood.dfdx(x_gp);
        gpgradx[N-2] *= -1.0/mu216; //derivative wrt A's transformation
 
        for(int j = 0; j < N-1; ++j) {
            for(int k = 0; k < j; ++k) {
                gsl_matrix_set(tmp,j,k,gpgradx[j]*gpgradx[k]);
                gsl_matrix_set(tmp,k,j,gpgradx[j]*gpgradx[k]);
            }
            gsl_matrix_set(tmp,j,j,gpgradx[j]*gpgradx[j]);
        }
        gsl_matrix_add(FIM,tmp);
    }
    gsl_matrix_free(tmp);

    gsl_matrix_set(FIM,N-1,N-1,2*NEXP);
    gsl_matrix_scale(FIM,1.0/var);

#ifdef DEBUG
    printf("_FIM:");  
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            printf("  %f  ", gsl_matrix_get(FIM,i,j));
        }
        printf("\n");
    }
    printf("\n");
#endif

    return FIM;

}

#endif//BLOOD_HPP
