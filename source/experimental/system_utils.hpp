#ifndef _SYSTEM_UTILS_HPP_
#define _SYSTEM_UTILS_HPP_

#include "coupled_ode_system.hpp"

void decouple(
    const std::vector<vec_d>& sol_coupled_sys,
    std::vector<vec_d>& equation_solution,
    std::vector<vec_d>& equation_sensitivities,
    int nr_equations,
    int parameter_dimension);

void printvec_s(const char * name, const vec_s & v);
void printvec_d(const char * name, const vec_d & v);

#endif // _SYSTEM_UTILS_HPP_
