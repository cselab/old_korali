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

GaussianProcess gpblood("gp_blood_cmaes_opt.txt");

double gpllk(const double* theta, int N) {

    size_t Nsamples = gpblood.get_sampleset_size();
    auto sampleset  = gpblood.get_sampleset();
    
    double var = theta[N-1]*theta[N-1];

    double x[N];
    std::copy(theta,theta+N-1,x);

    double sse = 0.0;
    for(size_t i = 0; i < Nsamples; ++i) {
        x[N-1] = sampleset->x(i).data()[N-1];
        sse+= pow((gpblood.f(x) - sampleset->y(i)),2);
    }

    return -0.5*Nsamples*log(2*M_PI*var)-0.5*sse/var;
}

gsl_vector * gpllk_grad(const double * theta, int N) {

    size_t Nsamples = gpblood.get_sampleset_size();
    auto sampleset  = gpblood.get_sampleset();
    
    gsl_vector * grad = gsl_vector_calloc(N);

    double var = theta[N-1]*theta[N-1];
    
    double x[N];
    std::copy(theta,theta+N-1,x);

    double sse = 0.0;
    for(size_t i = 0; i < Nsamples; ++i) {
        x[N-1] = sampleset->x(i).data()[N-1];
        Eigen::VectorXd gpgradx = gpblood.dfdx(x);
        sse+= pow((gpblood.f(x) - sampleset->y(i)),2);
        for(int k = 0; k < N-1; ++k) {
            gsl_vector_set(grad, k, (gpblood.f(x) - sampleset->y(i))*gpgradx[k]);
            sse += pow((gpblood.f(x) - sampleset->y(i)),2);
        }
    }

    gsl_vector_set(grad, N-1, -Nsamples*theta[N-1]+sse/theta[N-1]); //divide by var below
    gsl_vector_scale(grad,1.0/var);                                 //divide all elements by var

    return grad;
}

gsl_matrix * gpllk_FIM(const double * theta, int N) {

    size_t Nsamples = gpblood.get_sampleset_size();
    auto sampleset  = gpblood.get_sampleset();
    
    double var = theta[N-1]*theta[N-1];
    double x[N];
    
    std::copy(theta,theta+N-1,x);

    gsl_matrix * FIM = gsl_matrix_calloc(N,N);
    
    double tmp;
    for(size_t i = 0; i < Nsamples; ++i) {
        x[N-1] = sampleset->x(i).data()[N-1];
        Eigen::VectorXd gpgradx = gpblood.dfdx(x);
 
        for(int j = 0; j < N-1; ++j) {
            for(int k = 0; k < j; ++k) {
                tmp = gsl_matrix_get(FIM,j,k);
                gsl_matrix_set(FIM,j,k,tmp+gpgradx[j]*gpgradx[k]);
                gsl_matrix_set(FIM,k,j,tmp+gpgradx[j]*gpgradx[k]);
            }
            tmp = gsl_matrix_get(FIM,j,j);
            gsl_matrix_set(FIM,j,j,tmp+gpgradx[j]*gpgradx[j]);
        }
    }
    
    gsl_matrix_set(FIM,N-1,N-1,2*Nsamples);
    gsl_matrix_scale(FIM,1.0/var);
    return FIM;

}

int main(int argc, char** argv)
{

    Fitfun bloodgp = Fitfun(gpllk, gpllk_grad, gpllk_FIM);
 
    // maximize with CMA-ES
    //TmcmcEngine engine(&bloodgp, Standard, "blood__tmcmc.par", "blood_priors.par");
    TmcmcEngine engine(&bloodgp, Manifold, "blood_tmcmc.par", "blood_priors.par");
    printf("run Tmcmcmc.. \n");
    engine.run();

    
    printf("exit succefull \n");
    return 0;
}

