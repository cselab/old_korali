#ifndef TMCMC_OBJ_FMIN_HPP
#define TMCMC_OBJ_FMIN_HPP

#include "tmcmc_types.hpp"

namespace tmcmc {

    double tmcmc_objlogp(double x, const double *fj, int fn, 
            double pj, double zero);

    int fmincon(const double *fj, int fn, double pj, double objTol, 
            double *xmin, double *fmin, const optim_options& opt);

    int fminsearch(double const *fj, int fn, double pj, double objTol, 
            double *xmin, double *fmin, const optim_options& opt);
 
    int fzerofind(double const *fj, int fn, double pj, double objTol, 
        double *xmin, double *fmin, const optim_options& opt);

}//namespace tmcmc


#endif //TMCMC_OBJ_FMIN_HPP

