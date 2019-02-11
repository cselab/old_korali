#ifndef _SYSTEM_UTILS_HPP_
#define _SYSTEM_UTILS_HPP_

#include <stan/math.hpp>

namespace fitfun {

typedef stan::math::var scalar_t;
typedef std::vector<double>  vec_d;
typedef std::vector<scalar_t>  vec_s;


void decouple(
    const std::vector<vec_d>& sol_coupled_sys,
    std::vector<vec_d>& equation_solution,
    std::vector<vec_d>& equation_sensitivities,
    int nr_equations,
    int parameter_dimension);

void printvec_s(const char * name, const vec_s & v);
void printvec_d(const char * name, const vec_d & v);

}//namespace fitfun

#endif // _SYSTEM_UTILS_HPP_
