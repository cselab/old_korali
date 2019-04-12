#ifndef BLOOD_HPP
#define BLOOD_HPP

#include <stdio.h>
#include <algorithm>
#include <functional>

#include <libgp/libgp/include/core/gp.h>
#include <libgp/libgp/include/core/gp_utils.h>
#include <libgp/libgp/include/core/rprop.h>

#include "fitfun.hpp"

using namespace libgp;
using namespace fitfun;


#ifdef DATASET3
size_t NEXP  = 5; // num experimental data
double F[] = { 5.714981, 11.435213, 12.377169, 22.796679, 24.760591 }; //3
double OUT[] = { 0.029267, 0.026932, 0.023028, 0.021188, 0.021306 }; //3
#endif


#ifdef DATASET4
size_t NEXP  = 3; // num experimental data
double F[] = { 6.059418, 10.820390, 20.602023 }; //4
double OUT[] = { 0.035286, 0.032080, 0.030420 }; //4
#endif


#ifdef DATASET5
size_t NEXP  = 5; // num experimental data
double F[]   = { 6.555291, 9.832936, 13.110581, 16.388226, 19.665872 }; // 5: experimental data (Force)
double OUT[] = { 0.039364, 0.030558, 0.032415, 0.028117, 0.028064 };    // 5: experimental data (TTf)
#endif


#ifdef DATASET6
size_t NEXP  = 5; // num experimental data
double F[] = { 5.867770, 11.611925, 12.005029, 23.673504, 24.671042 }; //6
double OUT[] = { 0.035660, 0.031749, 0.031790, 0.028413, 0.028243 }; //6
#endif


#ifdef DATASET7
size_t NEXP  = 2; // num experimental data
double F[] = { 6.242403, 12.590044 }; //7
double OUT[] = { 0.031000, 0.028000 }; //7
#endif


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


#endif
