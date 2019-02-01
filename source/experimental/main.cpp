#include <stdio.h>

#include "auto_diff_fitfun.hpp"

int main(int argc, char* argv[])
{
    int info[1] = {0};
    double * out;
    

    vec_d theta = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double in[4] = {0.5, 0.5, 0.5, 0.5};
    
    AutoFitfun af;
    af.setParams(theta);

    vec_d t = { 0.0, 10.0 };
    std::vector<vec_d> obs(1);
    obs[0] = { 4.0, 8.0 };
    
    af.setObservations(t,obs);

    auto ret = af.fitfun(in, 8, out, info);

    //printf("result: %lf\n", tmp->loglike);

    return 0;
}
