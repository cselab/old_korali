#ifndef CMAES_UTILS_H
#define CMAES_UTILS_H

#include <cmaes.h>

double get_time();

void set_bounds( int verbose, const char *fname, double **p_lower_bound, double **p_upper_bound, int dim);
double load_pop_from_file(int verbose, int step, double * const* pop, double *arFunvals, int dim, int lambda, int *checkp);
void make_all_points_feasible( cmaes_t *evo, double * const *pop, double * lower_bound, double * upper_bound );
void print_the_best( cmaes_t evo, int step );
void write_pop_to_file( cmaes_t evo, double *arFunvals, double * const* pop, int step );
int is_there_enough_time( long job_max_time, double gt0, double dt );

#endif // CMAES_UTILS_H
