#include <math.h>

#include "coupled_ode_system.hpp"
#include "auto_diff_fitfun.hpp"

vec_s AutoFitfun::getModelIC(const vec_s & theta) const
{
    // here IC is independent of theta
    int dim = getDim();
    vec_s ic = vec_s(dim);
    ic[0] = 0.0;
    ic[1] = _obs[0][0]*_params[6];
    ic[2] = _obs[0][0]*(1.0-_params[6]);
    ic[3] = 0.0;
    return ic;
}

void AutoFitfun::evalModel(vec_s & dy_out, const vec_s & y_in, double t )
{
    // here model is independent of t
    dy_out[0] = -_params[0]*y_in[0];
    dy_out[1] = _params[3]*y_in[1]*( 1.0 - y_in[1]*y_in[2]*y_in[3]/100.0 ) +
                _params[4]*y_in[3] - _params[2]*y_in[1] -
                _params[0]*_params[1]*y_in[0]*y_in[1];
    dy_out[2] = _params[2]*y_in[1] - _params[0]*_params[1]*y_in[0]*y_in[2];
    dy_out[3] = _params[0]*_params[1]*y_in[0]*y_in[2] - _params[4]*y_in[3] -
                _params[5]*y_in[3];
}

