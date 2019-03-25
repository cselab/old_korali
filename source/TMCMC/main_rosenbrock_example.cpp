#include <stdio.h>

#include "fitfun.hpp"
#include "rosenbrock.hpp"
#include "engine_tmcmc.hpp"

using namespace tmcmc;
using namespace fitfun;

int main(int argc, char *argv[])
{
    Fitfun rosenbrock = Fitfun(*f_minusRosenbrock);

    TmcmcEngine engine(&rosenbrock, Standard, "tmcmc_rosenbrock.par", "priors_rosenbrock.par");
    engine.run();
    return 0;
}
