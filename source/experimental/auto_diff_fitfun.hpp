#ifndef _AUTO_DIFF_FITFUN_HPP_
#define _AUTO_DIFF_FITFUN_HPP_

#include "coupled_ode_system.hpp"

class AutoFitfun :  public CoupledOdeSystem
{

public:
    AutoFitfun() : CoupledOdeSystem (8, 4, 0.0 ) {};
    //return_type* operator () (double *x, int n, void* output, int * info);
private:

    vec_s getModelIC() const;
    void evalModel(vec_s & dyOut, const vec_s & y, double t);

};


#endif// _AutoDiffFitFun_HPP_
