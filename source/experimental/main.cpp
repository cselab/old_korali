#include <stdio.h>

#include "system_utils.hpp"
#include "auto_diff_fitfun.hpp"

int main(int argc, char* argv[])
{
    int info[1] = {0};
    double * out = nullptr;
    

    vec_d theta = {0.1, 0.15, 0.47, 1.4, 0.95, 2.0, 1.1, 0.01};
    
    AutoFitfun af;
    af.setParams(theta);

    vec_d t = { 0.0, 0.5 };
    std::vector<vec_d> obs(1);
    obs[0] = { 4.0, 8.0 };
    
    af.setObservations(t,obs);
    /*
    vec_s ic = af.getIC();
    printvec_s("ic",ic);

    vec_s dzOut(ic.size(),0.0);
    vec_s ics = { ic[0], ic[1], ic[2], ic[3] };
    af.step(ics, dzOut, 3.14);


    printvec_s("dzOut",dzOut);
    */


    double x[8] = {0.1, 0.15, 0.47, 1.4, 0.95, 2.0, 1.1, 0.01};
    return_type * ptr = af.fitfun(x, 8, out, info);
 


    printf("loglike: %lf\n", ptr->loglike);

    return 0;
}
