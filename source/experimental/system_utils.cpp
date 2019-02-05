#include "system_utils.hpp"

void decouple(
    std::vector<vec_d>& sol_coupled_sys,
    std::vector<vec_d>& sensitivities,
    std::vector<vec_d>& equation_solution,
    int nr_equations,
    int parameter_dimension)
{

    int total_nrObservations = sol_coupled_sys.size();
    sensitivities.resize(total_nrObservations);
    equation_solution.resize(total_nrObservations);

    for(size_t i = 0; i < total_nrObservations; ++i) {
        equation_solution[i] = vec_d( sol_coupled_sys[i].begin(),
                                      sol_coupled_sys[i].begin()+nr_equations);

        sensitivities[i] = vec_d( sol_coupled_sys[i].begin()+nr_equations,
                                  sol_coupled_sys[i].end());
    }
}

void printvec_s(const char * name, vec_s v) {
    printf("\n");
    for(int i = 0; i < v.size(); ++i) 
        printf("%s[%d]: %lf\n", name, i, v[i].val());
    printf("\n");
}


