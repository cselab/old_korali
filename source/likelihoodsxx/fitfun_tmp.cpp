#include <math.h>

#include "fitfun.hpp"

namespace fitfun{

    void fitfun_initialize(int argc, const  char **argv) {}

    double f_Rosenbrock(const double *x, int N) {
        int i;
        double s = 0.;
        for (i = 0; i < N-1; ++i)
            s += 100.*pow(x[i+1]-x[i]*x[i], 2) + pow(x[i]-1., 2);
        return s;
    }

    double fitfun(const double *x, int N, void *output, int *info) {
        return -f_Rosenbrock(x, N);
    }

    void fitfun_finalize() {}

}
