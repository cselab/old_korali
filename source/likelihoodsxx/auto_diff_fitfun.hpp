#ifndef _AUTO_DIFF_FITFUN_HPP_
#define _AUTO_DIFF_FITFUN_HPP_

#include "coupled_ode_system.hpp"

namespace fitfun {

class AutoFitfun :  public CoupledOdeSystem
{

public:

    AutoFitfun(double t0, int numparam, int odedim, int mala) :
        CoupledOdeSystem (t0, numparam, odedim, mala) {};

private:

    vec_s getModelIC_s(const vec_s & params) const;
    vec_s calculateObservable(const vec_s & solution) const;
    void evalModel_s(vec_s & dyOut, const vec_s & y, const vec_s & params, double t);
};

}//namespace fitfun

#endif// _AutoDiffFitFun_HPP_
