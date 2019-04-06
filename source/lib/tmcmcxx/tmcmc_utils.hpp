#ifndef TMCMC_UTILS_HPP
#define TMCMC_UTILS_HPP

#include "priors.hpp"
#include "tmcmc_types.hpp"

namespace tmcmc
{

#include <gsl/gsl_rng.h>

//TODO: align with cmaes (cmaes_utils) (DW)
#ifdef _USE_TORC_

#include <mpi.h>
#include <torc.h>

#else

#include <pthread.h>
#include <sys/time.h>


inline double torc_gettime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}

inline int torc_node_id()
{
    return 0;
}
inline int torc_num_nodes()
{
    return 1;
}

#ifdef _USE_OPENMP_
#include <omp.h>
inline int torc_i_worker_id()
{
    return omp_get_thread_num();
}
inline int torc_i_num_workers()
{
    return omp_get_max_threads();
}
inline int torc_worker_id()
{
    return omp_get_thread_num();
}
#else
inline int torc_i_worker_id()
{
    return 0;
}
inline int torc_i_num_workers()
{
    return 1;
}
inline int torc_worker_id()
{
    return 0;
}
#endif

#endif

// TODO: where should this go (DW)
static int g_nfeval;
static int l_nfeval = 0;
static int t_nfeval = 0;

static pthread_mutex_t feval_m = PTHREAD_MUTEX_INITIALIZER;

int get_nfc();
int get_tfc();
void inc_nfc();
void reset_nfc();

bool compar_desc(const sort_t p1, const sort_t p2);

double compute_sum(double *v, int n);

double compute_dot_product(double row_vector[], double vector[], int dim);

void compute_mat_product_vect(double *mat/*2D*/, double vect[], double res_vect[], double coef, int dim);

void inv_matrix(double *current_hessian/*2D*/, double *inv_hessian/*2D*/, int dim);

double scale_to_box(const double* point, double sc, const double* add_vec, const double *elbds, const double *eubds, int dims);

int in_rect(double *v1, double *v2, double *diam, double sc, int D);

int make_posdef(double *mat, int dim, int method);

void print_matrixi(const char *name, int *x, int n);

void print_matrix(const char *name, double *x, int n);

void print_matrix(const char *name, double *x, int n1, int n2);

void print_matrix_i(char *name, int *x, int n);

void print_matrix_2d(const char *name, double **x, int n1, int n2);

} // namespace tmcmc

#endif //TMCMC_UTILS_HPP
