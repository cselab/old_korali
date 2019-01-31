#ifndef _SYSTEM_UTILS_HPP_
#define _SYSTEM_UTILS_HPP_

#include "coupled_ode_system.hpp"

void decouple(
    std::vector<vec_d>& sol_coupled_sys,
    std::vector<vec_d>& sensitivities,
    std::vector<vec_d>& equation_solution,
    const int nr_equations,
    const int parameter_dimension);

#endif // _SYSTEM_UTILS_HPP_
