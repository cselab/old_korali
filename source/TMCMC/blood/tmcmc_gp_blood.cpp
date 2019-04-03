#include <stdio.h>
#include <functional>

#include <libgp/libgp/include/core/gp.h>
#include <libgp/libgp/include/core/gp_utils.h>
#include <libgp/libgp/include/core/rprop.h>

#include "fitfun.hpp"
#include "engine_tmcmc.hpp"

using namespace libgp;
using namespace fitfun;
using namespace tmcmc;

GaussianProcess gpblood("gp_blood.txt");

void testGp(GaussianProcess& gp, double low, double up, size_t N);

double func(double* hyp, int N, void*, int*) {
    Eigen::Map<Eigen::VectorXd> hypvec(hyp, N);
    gpblood.covf().set_loghyper(hypvec);
    return -gpblood.log_likelihood();
}

int main(int argc, char** argv)
{

    printf("loglike before maximization %f \n", gpblood.log_likelihood());

    // init fitfun
    auto gpllk = [] (const double *theta, int N) { return 0.0; };

    /*
    // init gradient
    auto gpllkGrad = [&gp] (const double *theta, int N) { 
                                Eigen::VectorXd gradx = gp.dfdx(theta); 
                                gsl_vector* gsl_grad  = gsl_vector_alloc(N);
                                for(int i = 0; i < N; ++i) gsl_vector_set(gsl_grad, i, gradx[i]); 
                                return gsl_grad; 
    };

    // init Fisher Information Matrix
    auto gpllkFIM = [&gp] (const double *theta, int N) {
                double var = gp.var(theta);
                Eigen::VectorXd gradx = gp.dfdx(theta); 
                gsl_matrix * FIM      = gsl_matrix_calloc(N,N);
                for(size_t i = 0; i < N; ++i) {
                    for(size_t j = 0; j < i; ++j) gsl_matrix_set(FIM, i, j, gradx[i] * gradx[j] / var);
                    gsl_matrix_set(FIM, i, i, gradx[i] * gradx[i] / var);
                }
                return FIM; 
    };
    */

    Fitfun bloodgp = Fitfun(gpllk);
 
    // maximize with CMA-ES
    TmcmcEngine engine(&bloodgp, Standard, "gp_tmcmc.par", "gp_priors.par");
    printf("run Tmcmcmc.. \n");
    engine.run();

    
    printf("exit succefull \n");
    return 0;
}

