#include <stdio.h>

#include <libgp/libgp/include/core/gp.h>
#include <libgp/libgp/include/core/gp_utils.h>
#include <libgp/libgp/include/core/rprop.h>

using namespace libgp;

void testGp(GaussianProcess& gp, double low, double up, size_t N);

int main(int argc, char** argv)
{
    size_t nOptSteps = 10000;

    printf("initializing Gaussian process .. \n");

    GaussianProcess gpblood("gp_blood.txt");
    printf("loglike before maximization %f \n", gpblood.log_likelihood());

    // maximize with RProp
    RProp rprop;
    rprop.init(1e-6);
    printf("maximize hyper params of Gp ( max iterations: %zu ) .. \n", nOptSteps);
    rprop.maximize(&gpblood, nOptSteps, true);
    printf("loglike after maximization %f \n", gpblood.log_likelihood());

    gpblood.write("gp_blood_opt.txt");
}

