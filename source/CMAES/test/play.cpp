#include <vector>

#include "catch.hpp"
#include "engine_cmaes.hpp"

extern "C" {
#include "fitfun_tests.h"
}

#define EPSILON_TEST 1e-5

void print_best(CmaesEngine& engine);

double fitfun(double *x, int N, void *output, int *info)
{
    return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

int main(int argc, char** argv)
{

    auto engine = CmaesEngine(&fitfun, "./playground/");
    engine.run();

    print_best(engine);
    return 0;
}


void print_best(CmaesEngine& engine)
{

    auto evo = engine.getEvo();
    double bestFunVal = cmaes_Get(evo,"fbestever");

    printf("fbestever: %f\n\n",bestFunVal);

    for(int i = 0; i < evo->sp.N; ++i) {
        printf("rgxbestever[%d]: %f\n", i, evo->rgxbestever[i]);
    }
}
