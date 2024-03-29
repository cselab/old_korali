#include <stdio.h>

#include <libgp/libgp/include/core/gp.h>
#include <libgp/libgp/include/core/gp_utils.h>
#include <libgp/libgp/include/core/rprop.h>

#include "fitfun.hpp"
#include "rosenbrock.hpp"
#include "engine_tmcmc.hpp"

using namespace libgp;
using namespace tmcmc;
using namespace fitfun;

double low = -6;
double up  = 6;

size_t nTraining = 100;
size_t nTest     = 100;
size_t nOptSteps = 5000;

void testGp(GaussianProcess& gp, double low, double up, size_t N);

int main(int argc, char *argv[])
{    
   
    // initialize Gaussian process for 2-D input using the squared exponential 
    printf("initializing Gaussian process .. \n");
    GaussianProcess gp(2, "CovSum ( CovSEiso, CovNoise)");
    // initialize hyper parameter vector
    Eigen::VectorXd params(gp.covf().get_param_dim()); // 3 (2 for exponential + 1 for noise)
    params << 0.0, 0.0, -1.0;
    // set parameters of covariance function
    gp.covf().set_loghyper(params);

    // add training patterns
    printf("add training data to Gp ( N = %zu )..  \n", nTraining);
    double y;
    double x[2];
    for(size_t i = 0; i < nTraining; ++i) {
      x[0] = drand48()*(up-low)+low;
      x[1] = drand48()*(up-low)+low;
      y = f_minusRosenbrock(x, 2);
      gp.add_pattern(x, y);
    }
    printf("loglike before maximization %f \n", gp.log_likelihood());
    testGp(gp, low, up, nTest);

    // maximize with RProp
    RProp rprop;
    rprop.init(1e-4);
    printf("maximize hyper params of Gp ( max iterations: %zu) .. \n", nOptSteps);
    rprop.maximize(&gp, nOptSteps, false);
    printf("loglike after maximization %f \n", gp.log_likelihood());
    testGp(gp, low, up, nTest);

    // init fitfun
    auto gpllk = [&gp] (const double *theta, int N) { return gp.f(theta); };

    auto gpllkGrad = [&gp] (const double *theta, int N) { 
                                Eigen::VectorXd gradx = gp.dfdx(theta); 
                                gsl_vector* gsl_grad  = gsl_vector_alloc(N);
                                for(int i = 0; i < N; ++i) gsl_vector_set(gsl_grad, i, gradx[i]); 
                                return gsl_grad; 
    };

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
    
    Fitfun rosenbrockgp = Fitfun(gpllk, gpllkGrad, gpllkFIM);

    //tmcmcm
    printf("run Tmcmcmc.. \n");
    TmcmcEngine engine(&rosenbrockgp, Manifold, "tmcmc_rosenbrock.par", "priors_rosenbrock.par");
    engine.run();
    
    //exit
    printf("exit succefull \n");
    return 0;
}

void testGp(GaussianProcess& gp, double low, double up, size_t N){
    printf("calculating MSE ( N = %zu ) .. \n", N);
    double max = -std::numeric_limits<double>::max();
    double e, se = 0;
    double x[2] , xmax[2];
    for(size_t i = 0; i < N; ++i) {
        x[0] = drand48()*(up-low)+low;
        x[1] = drand48()*(up-low)+low;
        
        e = pow(f_minusRosenbrock(x,2) - gp.f(x), 2);
        if(e > max) { max = e; xmax[0] = x[0]; xmax[1] = x[1]; }
        
        se += pow(f_minusRosenbrock(x, 2)-gp.f(x),2);
    }
    printf("mse: %f (max squared error %f at (%f, %f) ) \n", se/N, max, xmax[0], xmax[1] );
}
