#include <stdio.h>

#include "auto_diff_fitfun.hpp" //this must go before engine_tmcmc.hpp (help!)
#include "fitfun.hpp"
#include "engine_tmcmc.hpp"

using namespace tmcmc;
using namespace fitfun;

int main(int argc, char *argv[])
{
    
    vec_d t = { 0.0, 9.0, 12.0, 17.0, 25.0 };
    
    std::vector<vec_d> obs(t.size());
    obs[0] = { 0.8650 };
    obs[1] = { 1.0846 };
    obs[2] = { 1.0532 };
    obs[3] = { 0.8838 };
    obs[4] = { -0.7374 };

    AutoFitfun af(3,2,true); //false := no MALA
    //af.setParams(theta);
    af.setObservations(t,obs);

    TmcmcEngine engine(&af);
    engine.run();
    return 0;
}
