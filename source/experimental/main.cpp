#include <stdio.h>

#include "auto_diff_fitfun.hpp"

int main(int argc, char* argv[])
{

    double x = 1.0;
    int n = 1;
    int info = 0;

    double out[1] = {0.0};

    vec_d theta = vec_d(8);

    AutoFitfun af(theta,0);

    //auto tmp = af(&x, n, (void*) out, &info);

    //printf("result: %lf\n", tmp->loglike);

    return 0;
}
