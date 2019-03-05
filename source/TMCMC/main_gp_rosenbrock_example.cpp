#include <stdio.h>

#include <libgp/libgp/include/gp.h>
#include <libgp/libgp/include/gp_utils.h>
#include <libgp/libgp/include/rprop.h>

#include "fitfun.hpp"
#include "rosenbrock.hpp"
#include "engine_tmcmc.hpp"

using namespace libgp;
using namespace tmcmc;
using namespace fitfun;

int main(int argc, char *argv[])
{    
    int n = 50;
    double low = -3;
    double up  = 3;

    // initialize Gaussian process for 1-D input using the squared exponential 
    printf("initializing Gaussian process .. \n");
    GaussianProcess gp(2, "CovSum ( CovSEiso, CovNoise)");
    // initialize hyper parameter vector
    Eigen::VectorXd params(gp.covf().get_param_dim()); // 3 (2 for exponential + 1 for noise)
    params << 0.1, 0.0, 1.0;
    // set parameters of covariance function
    gp.covf().set_loghyper(params);

    // add training patterns
    printf("add training data to Gp ..  \n");
    double y;
    double x[2];
    for(int i = 0; i < n; ++i) {
      x[0] = drand48()*(up-low)-low;
      x[1] = drand48()*(up-low)-low;
      y = f_minusRosenbrock(x, 2);
      gp.add_pattern(x, y);
    }

    // maximize with RProp
    /*
    RProp rprop;
    rprop.init(1e-4);
    bool verbose = true;
    int nit = 50;
    printf("maximize hyper params of Gp (verbose: %d, iterations: %d) .. \n", verbose, nit);
    rprop.maximize(&gp, 50, nit);
    */

    // init fitfun
    auto gplambda = [&gp] (const double *theta, int N) { gp.covf().set_loghyper(theta); return gp.log_likelihood(); };
    Fitfun rosenbrock = Fitfun(gplambda);

    double theta[3];
    for(int i = 0; i < n; ++i) {
        theta[0] = drand48()*6-3;
        theta[1] = drand48()*6-3;
        theta[2] = drand48()*6-3;
        printf("test fitfun: %f\n", rosenbrock.evaluate(theta,3, NULL, NULL));
    }

    //tmcmcm
    printf("run Tmcmcmc.. \n");
    TmcmcEngine engine(&rosenbrock, Standard, "tmcmc_rosenbrock.par", "priors_rosenbrock.par");
    engine.run();
    
    //exit
    printf("exit succefull \n");
    return 0;
}