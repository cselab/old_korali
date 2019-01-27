#ifndef TMCMC_UTILS_HPP
#define TMCMC_UTILS_HPP

#include "priors.hpp"
#include "tmcmc_types.hpp"

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
    
    // TODO: where should this go (DW)
    static int g_nfeval;
    static int l_nfeval = 0;
    static int t_nfeval = 0;

    static pthread_mutex_t feval_m = PTHREAD_MUTEX_INITIALIZER;

    int  get_nfc();
    void inc_nfc();
    void reset_nfc();
    
    int compar_desc(const void *p1, const void *p2);

    int in_rect(double *v1, double *v2, double *diam, double sc, int D);

    void print_matrix(const char *name, double *x, int n);
    //void print_matrix_i(char *name, int *x, int n);
    void print_matrix_2d(const char *name, double **x, int n1, int n2);
 
} // namespace tmcmc

#endif //TMCMC_UTILS_HPP
