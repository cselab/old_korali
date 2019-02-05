#include <stdio.h>

#include "system_utils.hpp"
#include "auto_diff_fitfun.hpp"

int main(int argc, char* argv[])
{
    int info[1] = {0};
    double * out;
    

    vec_d theta = {0.1, 0.15, 0.47, 1.4, 0.95, 2.0, 1.1, 0.01};
    double in[4] = {0.5, 0.5, 0.5, 0.5};
    
    AutoFitfun af;
    af.setParams(theta);

    vec_d t = { 0.0, 10.0 };
    std::vector<vec_d> obs(1);
    obs[0] = { 4.0, 8.0 };
    
    af.setObservations(t,obs);

    vec_s ic = af.getIC();
    printvec_s("ic",ic);

    vec_s dzOut(4,0.0);
    vec_s ics = { ic[0], ic[1], ic[2], ic[3] };
    af.step(ics, dzOut, 3.14);


    printvec_s("dzOut",dzOut);


    


    //printf("result: %lf\n", tmp->loglike);

    return 0;
}
