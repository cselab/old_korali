#ifndef TMCMC_UTILS_HPP
#define TMCMC_UTILS_HPP

#include <gsl/gsl_rng.h>

#include "tmcmc_types.hpp"

// TODO: where to go with this (singleton?) (DW)
extern gsl_rng   				**r;
extern int   					*local_seed;

int mvnrnd(double *mean, double *sigma, double *out, int N);

double uniformrand(double a, double b);

void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[]);

typedef struct sort_s {
    int idx;
    int nsel;
    double F;
} sort_t;

int compar_desc(const void *p1, const void *p2);


//TODO: align with cmaes (cmaes_utils) (DW)
#ifdef _USE_TORC_

	#include <mpi.h>
    #include <torc.h>

#else

    #include <pthread.h>
	int torc_node_id();
	int torc_num_nodes();

	#ifdef _USE_OPENMP_
		int torc_i_worker_id();
		int torc_i_num_workers();
		int torc_worker_id();
	#else
		int torc_i_worker_id();
		int torc_i_num_workers();
		int torc_worker_id();
	#endif

	double torc_gettime();

#endif

    int in_rect(double *v1, double *v2, double *diam, double sc, int D);
    
    void print_matrix(const char *name, double *x, int n);
    //void print_matrix_i(char *name, int *x, int n);
    void print_matrix_2d(const char *name, double **x, int n1, int n2);


#endif //TMCMC_UTILS_HPP
