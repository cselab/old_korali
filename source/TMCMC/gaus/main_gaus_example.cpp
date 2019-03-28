#include <stdio.h>

#include "fitfun.hpp"
#include "gaussian.hpp"
#include "engine_tmcmc.hpp"

#include "error_helpers.hpp"

using namespace tmcmc;
using namespace fitfun;

int main(int argc, char *argv[])
{
    int NRUNS = 100;

    int N = 5;
    double LB = -20;
    double UB = 20;

    double *le     = new double[NRUNS];

    for(int i = 0; i < NRUNS; ++i) 
    {
        Fitfun gaus = Fitfun(*f_gaus, *gradX_gaus, *hessX_gaus);
        TmcmcEngine engine(&gaus, Standard, "gaus_tmcmc.par", "gaus_priors.par");
        engine.run();
        le[i] = engine.getLogEvidence();
    }
    
    print_err_le(le, N, LB, UB, NRUNS);
    //print_err1(mean, SS, LB, UB, N);
        
    delete [] le;

    return 0;
}
