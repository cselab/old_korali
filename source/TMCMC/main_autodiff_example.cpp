#include <stdio.h>

#include "auto_diff_fitfun.hpp" //this must go before engine_tmcmc.hpp (help!)
#include "fitfun.hpp"
#include "engine_tmcmc.hpp"

using namespace tmcmc;
using namespace fitfun;

int main(int argc, char *argv[])
{
    
    //vec_d t = { 0.0, 9.0, 12.0, 17.0, 25.0 };
    vec_d t = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    
    std::vector<vec_d> obs(t.size());
    obs[0] = { 114.46 };
    obs[1] = { 148.65 };
    obs[2] = { 193.09 };
    obs[3] = { 269.63 };
    obs[4] = { 386.05 };
    
    AutoFitfun af(0.0, 2, 2, true); //false := no MALA
    af.setObservations(t,obs);

    TmcmcEngine engine(&af, Manifold);
    engine.run();
    return 0;
}
