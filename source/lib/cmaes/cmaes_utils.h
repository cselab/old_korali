#ifndef CMAES_UTILS_H
#define CMAES_UTILS_H

#include "cmaes.h"

double get_time();

int	   cmaes_utils_file_exists( const char *fname);
void   cmaes_utils_read_bounds( int verbose, const char *fname, double **p_lower_bound, double **p_upper_bound, int dim);
double cmaes_utils_load_pop_from_file(int verbose, int step, double * const* pop, double *arFunvals, int dim, int lambda, int *checkp);
void   cmaes_utils_make_all_points_feasible( cmaes_t *evo, double * const *pop, double * lower_bound, double * upper_bound );
void   cmaes_utils_print_the_best( cmaes_t evo, int step );
void   cmaes_utils_write_pop_to_file( cmaes_t evo, double *arFunvals, double * const* pop, int step );
int    cmaes_utils_is_there_enough_time( long job_max_time, double gt0, double dt );

#endif // CMAES_UTILS_H
