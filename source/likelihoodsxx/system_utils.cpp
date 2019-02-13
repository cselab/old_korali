#include <vector>
#include "system_utils.hpp"

namespace fitfun{

void decouple(
    const std::vector<vec_d> & sol_coupled_sys,
    std::vector<vec_d>& equation_solution,
    std::vector<vec_d>& equation_sensitivities,
    int nr_equations,
    int parameter_dimension)
{
    size_t total_nrObservations = sol_coupled_sys.size();
    equation_sensitivities.resize(total_nrObservations);
    equation_solution.resize(total_nrObservations);

    for(size_t i = 0; i < total_nrObservations; ++i) {
        equation_solution[i] = vec_d( sol_coupled_sys[i].begin(),
                                      sol_coupled_sys[i].begin()+nr_equations);

        equation_sensitivities[i] = vec_d( sol_coupled_sys[i].begin()+nr_equations,
                                           sol_coupled_sys[i].end());
    }
}

void printvec_s(const char * name, const vec_s & v)
{
    printf("\n");
    for(int i = 0; i < v.size(); ++i)
        printf("%s[%d]: %lf\n", name, i, v[i].val());
    printf("\n");
}

void printvec_d(const char * name, const vec_d & v)
{
    printf("\n");
    for(int i = 0; i < v.size(); ++i)
        printf("%s[%d]: %lf\n", name, i, v[i]);
    printf("\n");
}

}//namespace fitfun
