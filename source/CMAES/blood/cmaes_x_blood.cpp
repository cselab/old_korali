#include <fstream>
#include <stdio.h>
#include <functional>

#include <libgp/libgp/include/core/gp.h>
#include <libgp/libgp/include/core/gp_utils.h>
#include <libgp/libgp/include/core/rprop.h>

#include "blood.hpp"
#include "fitfun.hpp"
#include "engine_cmaes.hpp"

using namespace libgp;
using namespace fitfun;

void write_best(CmaesEngine& engine);

int main(int argc, char** argv)
{
    auto gpf = [] (double* theta, int N, void*, int*) { return -gpllk(theta, N); };
    auto engine = CmaesEngine(gpf, "./val/");
    engine.run();
    write_best(engine); 
}


void write_best(CmaesEngine& engine)
{
    FILE * pFile = fopen("opt.txt", "a");
    if (pFile == NULL) { printf("ERROR: could not open file!!!"); return; }
    
    auto evo = engine.getEvo();
    double bestFunVal = cmaes_Get(evo,"fbestever");

    for(int i = 0; i < evo->sp.N; ++i) {
        fprintf(pFile, "%.15f\t", evo->rgxbestever[i]);
    }
    fprintf(pFile, "%.15f\n", bestFunVal);

    fclose (pFile); 
    return;
}
