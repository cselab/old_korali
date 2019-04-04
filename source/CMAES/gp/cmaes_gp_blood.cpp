#include <stdio.h>
#include <functional>

#include <libgp/libgp/include/core/gp.h>
#include <libgp/libgp/include/core/gp_utils.h>
#include <libgp/libgp/include/core/rprop.h>

#include "fitfun.hpp"
#include "rosenbrock.hpp"
#include "engine_cmaes.hpp"

using namespace libgp;
using namespace fitfun;

GaussianProcess gpblood("gp_blood.txt");

void print_best(CmaesEngine& engine);
void testGp(GaussianProcess& gp, double low, double up, size_t N);

double func(double* hyp, int N, void*, int*) {
    Eigen::Map<Eigen::VectorXd> hypvec(hyp, N);
    gpblood.covf().set_loghyper(hypvec);
    return -gpblood.log_likelihood();
}

int main(int argc, char** argv)
{

    printf("loglike before maximization %f \n", gpblood.log_likelihood());

    
    // maximize with CMA-ES
    auto engine = CmaesEngine(&func, "./run/");
    engine.run();

    print_best(engine);
    
    printf("loglike after maximization %f \n", gpblood.log_likelihood());

    gpblood.write("gp_blood_cmaes_opt.txt");
}


void print_best(CmaesEngine& engine)
{

    auto evo = engine.getEvo();
    //double bestFunVal = cmaes_Get(evo,"fbestever");

    //printf("fbestever: %f\n\n",bestFunVal);

    for(int i = 0; i < evo->sp.N; ++i) {
        printf("rgxbestever[%d]: %f\n", i, evo->rgxbestever[i]);
    }
}
