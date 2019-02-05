#ifndef _SYSTEM_UTILS_HPP_
#define _SYSTEM_UTILS_HPP_

#include "coupled_ode_system.hpp"

void decouple(
    std::vector<vec_d>& sol_coupled_sys,
    std::vector<vec_d>& sensitivities,
    std::vector<vec_d>& equation_solution,
    int nr_equations,
    int parameter_dimension);

void printvec_s(const char * name, const vec_s & v);

#endif // _SYSTEM_UTILS_HPP_
