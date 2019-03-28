#include <stdio.h>

#include "fitfun.hpp"
#include "gaussian.hpp"
#include "engine_tmcmc.hpp"

using namespace tmcmc;
using namespace fitfun;

int main(int argc, char *argv[])
{

    int N = 5;
    double LB = -10.0;
    double UB = 10.0;

    Fitfun gaus = Fitfun(*f_gaus, *gradX_gaus, *hessX_gaus);

    TmcmcEngine engine(&gaus, Standard, "gaus_tmcmc.par", "gaus_priors.par");
    
    engine.run();

    double *mean = engine.getNewMean();
    double **SS = engine.getNewSampleCov();

    double err = 0;
    for(int i = 0; i < N; ++i) {
        err += fabs(mean[i])/N;
        for(int j = 0; j < N; ++j) {
            if (i == j) err += fabs(1 - SS[i][j])/(N*N);
            else err += fabs(SS[i][j])/(N*N);
        }
    }

    delete mean;
    for(int i = 0; i < 5; ++i) delete SS[i];
    delete SS;

    double explogevidence = -0.5*N*log(2.0*M_PI) - N*log(UB-LB);
    printf("err: %f\n",err);
    printf("expected log evidence: %f\n", explogevidence);
        
    return 0;
}
