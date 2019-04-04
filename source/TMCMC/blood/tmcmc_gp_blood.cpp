#include "blood.hpp"

int main(int argc, char** argv)
{

    Fitfun bloodgp = Fitfun(gpllk);
 
    TmcmcEngine engine(&bloodgp, Standard, "tmcmc_blood.par", "priors_blood.par");
    printf("run Tmcmcmc.. \n");
    engine.run();

    
    printf("exit succefull \n");
    return 0;
}

