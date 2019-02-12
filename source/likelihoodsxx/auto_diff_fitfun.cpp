#include <math.h>

#include "coupled_ode_system.hpp"
#include "auto_diff_fitfun.hpp"

namespace fitfun {

vec_s AutoFitfun::getModelIC_s(const vec_s & params) const
{
    vec_s ic = vec_s(getDim());
    
    /* Simple Model */
    ic[0] = params[1];
    ic[1] = 0;
    
    return ic;
}

void AutoFitfun::evalModel_s(vec_s & dy_out, const vec_s & y_in, const vec_s & params, double t )
{
    /* Simple Model */
    dy_out[0] = params[0] * stan::math::cos(t) - 2.0*params[1]*params[2]*y_in[1];
    dy_out[1] = 2.0*params[2]*(y_in[0]-params[0]*stan::math::sin(t))/params[1];   
}

vec_s AutoFitfun::calculateObservable(const vec_s & equation_solution) const
{
    vec_s observable(getObsDim());
    
    // Simple Model
    observable[0] = equation_solution[0] + equation_solution[1];
    
    return observable;
}

/* FUNCTION DEFINITIONS FOR 
 * TUMOR GROWTH MODEL

vec_s AutoFitfun::getModelIC_s(const vec_s & params) const
{
    vec_s ic = vec_s(getDim());
    
    ic[0] = 1.0;
    ic[1] = _obs[0][0]*_params[6];
    ic[2] = _obs[0][0]*(1.0-_params[6]);
    ic[3] = 0.0;
    
    return ic;
}

void AutoFitfun::evalModel_s(vec_s & dy_out, const vec_s & y_in, const vec_s & params, double t )
{
    dy_out[0] = -_params[0]*y_in[0];
    dy_out[1] = _params[3]*y_in[1]*( 1.0 - y_in[1]*y_in[2]*y_in[3]/100.0 ) +
                _params[4]*y_in[3] - _params[2]*y_in[1] - _params[0]*_params[1]*y_in[0]*y_in[1];
    dy_out[2] = _params[2]*y_in[1] - _params[0]*_params[1]*y_in[0]*y_in[2];
    dy_out[3] = _params[0]*_params[1]*y_in[0]*y_in[2] - _params[4]*y_in[3] - _params[5]*y_in[3];
}

vec_s AutoFitfun::calculateObservable(const vec_s & equation_solution) const
{
    vec_s observable(getObsDim());
    
    observable[0] = equation_solution[1] + equation_solution[2] + equation_solution[3];
    
    return observable;
}
*/


}//namespace fitfun
