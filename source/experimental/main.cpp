#include <stdio.h>

#include "system_utils.hpp"
#include "auto_diff_fitfun.hpp"

int main(int argc, char* argv[])
{
    int info[1] = {0};
    double * out = nullptr;
    

    vec_d theta = { -0.8, 1, 0.5 };
    //vec_d theta = {0.1, 0.15, 0.47, 1.4, 0.95, 2.0, 1.1, 0.01};
    
    AutoFitfun af(3,2,true);
    af.setParams(theta);

    vec_d t = { 0.0, 9.0 };
    //vec_d t = { 12.0, 14.09, 16.17, 18.26, 20.35, 22.43 }; 
    std::vector<vec_d> obs(t.size());
    obs[0] = { 2.0 };
    obs[1] = { 4.0 };
    /*
    obs[0] = { 50.17 };
    obs[1] = { 73.0 };
    obs[2] = { 80.0 };
    obs[3] = { 73.0 };
    obs[4] = { 73.0 };
    obs[5] = { 70.0 };
    */
    af.setObservations(t,obs);

    double x[4] = { -0.8, 1, 0.5, 0.5 };
    //double x[9] = {0.1, 0.15, 0.47, 1.4, 0.95, 2.0, 1.1, 0.01, 1.0};
    return_type * ptr = af.fitfun(x, 4, out, info);
 


    printf("loglike: %lf\n", ptr->loglike);

    return 0;
}
