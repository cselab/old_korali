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
    double low = -20;
    double up  = 20;

    int n=1000, m=5000;
    // initialize Gaussian process for 3-D input using the squared exponential 
    GaussianProcess gp(1, "CovSum ( CovSEiso, CovNoise)");
    //GaussianProcess gp(3, "CovSEiso");
    // initialize hyper parameter vector
    Eigen::VectorXd params(gp.covf().get_param_dim());
    params << 0.0, 0.0, 1.0;
    // set parameters of covariance function
    gp.covf().set_loghyper(params);

    GaussianProcess gpOpt(gp);

    // add training patterns
    double y;
    double x[1];
    for(int i = 0; i < n; ++i) {
      x[0] = i*(up-low)/n - low;
      y = f_minusRosenbrock(x, 1);
      gp.add_pattern(x, y);
    }


    Fitfun rosenbrock = Fitfun(*f_minusRosenbrock);

    TmcmcEngine engine(&rosenbrock);
    engine.run();
    return 0;
}
