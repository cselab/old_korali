#ifndef _AUTO_DIFF_FITFUN_HPP_
#define _AUTO_DIFF_FITFUN_HPP_

#include "coupled_ode_system.hpp"

class AutoFitfun :  public CoupledOdeSystem
{

public:
    AutoFitfun(vec_d params, int nobs) : CoupledOdeSystem (4,8, params, 0.0, nobs) {};
    return_type* operator () (double *x, int n, void* output, int * info);

private:

    vec_s getModelIC(const vec_s & theta) const;
    void evalModel(vec_s & dyOut, double t, const vec_s & y);

};


#endif// _AutoDiffFitFun_HPP_
