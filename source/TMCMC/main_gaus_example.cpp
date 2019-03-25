#include <stdio.h>

#include "fitfun.hpp"
#include "gaussian.hpp"
#include "engine_tmcmc.hpp"

using namespace tmcmc;
using namespace fitfun;

int main(int argc, char *argv[])
{
    Fitfun gaus = Fitfun(*f_gaus, *gradX_gaus, *hessX_gaus);

    TmcmcEngine engine(&gaus, Standard, "gaus_tmcmc.par", "gaus_priors.par");
    //TmcmcEngine engine(&gaus, Manifold, "gaus_tmcmc.par", "gaus_priors.par");
    
    engine.run();
    return 0;
}
