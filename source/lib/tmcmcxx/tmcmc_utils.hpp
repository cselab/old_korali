#ifndef TMCMC_UTILS_HPP
#define TMCMC_UTILS_HPP

#include "priors.h"

namespace tmcmc {

    #include <gsl/gsl_rng.h>

//TODO: align with cmaes (cmaes_utils) (DW)
#ifdef _USE_TORC_

    #include <mpi.h>
    #include <torc.h>

#else

    #include <pthread.h>
	#include <sys/time.h>
 

    inline double torc_gettime(){
    	struct timeval t;
    	gettimeofday(&t, NULL);
    	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
	}   
    
    inline int torc_node_id() { return 0; }
	inline int torc_num_nodes() { return 1; }
    
#ifdef _USE_OPENMP_
    #include <omp.h>
    inline int torc_i_worker_id() { return omp_get_thread_num(); }
    inline int torc_i_num_workers() { return omp_get_max_threads(); }
    inline int torc_worker_id() { return omp_get_thread_num(); }
#else
    inline int torc_i_worker_id() { return 0; }
    inline int torc_i_num_workers() { return 1; }
    inline int torc_worker_id() { return 0; }
#endif 

#endif

    // TODO: where to go with this (singleton?) (DW)
    static const gsl_rng_type *gsl_rng_t;
    static gsl_rng            **r;
    static int   	          *local_seed;

    static int g_nfeval;
    static int l_nfeval = 0;
    static int t_nfeval = 0;

    static pthread_mutex_t feval_m = PTHREAD_MUTEX_INITIALIZER;

    int  get_nfc();
    void inc_nfc();
    void reset_nfc();

    void gsl_rand_init(int seed);
    void call_gsl_rand_init(int seed);
    void spmd_gsl_rand_init(int seed);

    typedef struct sort_s {
        int idx;
        int nsel;
        double F;
    } sort_t;

    int compar_desc(const void *p1, const void *p2);

    int in_rect(double *v1, double *v2, double *diam, double sc, int D);

    void print_matrix(const char *name, double *x, int n);
    //void print_matrix_i(char *name, int *x, int n);
    void print_matrix_2d(const char *name, double **x, int n1, int n2);
 
    void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[]);
    
    int mvnrnd(double *mean, double *sigma, double *out, int N);
   
   
   /*
    double eval_random( Density d );
    double uniformrand(double a, double b);
    */
 
} // namespace tmcmc

#endif //TMCMC_UTILS_HPP
