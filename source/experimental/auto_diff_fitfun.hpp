#ifndef _AUTO_DIFF_FITFUN_HPP_
#define _AUTO_DIFF_FITFUN_HPP_

#include "coupled_ode_system.hpp"

class AutoFitfun :  public CoupledOdeSystem
{

public:
    AutoFitfun(int numparam, int odedim, int mala) : 
        CoupledOdeSystem (numparam, odedim, mala) {};
    //return_type* operator () (double *x, int n, void* output, int * info);
private:

    vec_s getModelIC_s(const vec_s & params) const;
    void evalModel_s(vec_s & dyOut, const vec_s & y, const vec_s & params, double t);
    
    vec_d getModelIC(const vec_d & params) const;
    void evalModel(vec_d & dyOut, const vec_d & y, const vec_d & params, double t);
    
    vec_s calculateObservable(const vec_s & solution);
};


#endif// _AutoDiffFitFun_HPP_
