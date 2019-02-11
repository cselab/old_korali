#include <stdio.h>

#include "auto_diff_fitfun.hpp" //this must go before engine_tmcmc.hpp (help!)
#include "fitfun.hpp"
#include "engine_tmcmc.hpp"

using namespace tmcmc;
using namespace fitfun;

int main(int argc, char *argv[])
{
    vec_d theta = {-0.8, 1, 0.5};
    
    vec_d t = { 0.0, 9.0 };
    
    std::vector<vec_d> obs(t.size());
    obs[0] = { 2.0 };
    obs[1] = { 4.0 };

    AutoFitfun af(3,2,true); //false := no MALA
    af.setParams(theta);
    af.setObservations(t,obs);

    TmcmcEngine engine(&af);
    engine.run();
    return 0;
}
