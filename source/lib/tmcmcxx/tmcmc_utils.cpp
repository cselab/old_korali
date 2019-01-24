#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

#include "tmcmc_utils.hpp"

namespace tmcmc {

    int in_rect(double *v1, double *v2, double *diam, double sc, int D) {
        int d;
        for (d = 0; d < D; ++d) {
            if (fabs(v1[d]-v2[d]) > sc*diam[d]) return 0;
        }
        return 1;
    }


    int compar_desc(const void* p1, const void* p2)
    {
        int dir = +1;   /* -1: ascending order, +1: descending order */
        sort_t *s1 = (sort_t *) p1;
        sort_t *s2 = (sort_t *) p2;

        if (s1->nsel < s2->nsel) return dir;
        if (s1->nsel > s2->nsel) return -dir;
        /*    if (s1->nsel == s2->nsel) return 0;*/
        return 0;
    }


    void get_nfc_task(int *x)
    {
        *x = l_nfeval;
    }


    int get_nfc() {
        int c[1024]; /* MAX_NODES*/
#ifdef _USE_TORC_
        for (int i = 0; i < torc_num_nodes(); ++i) {
            torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())get_nfc_task, 1,
                    1, MPI_INT, CALL_BY_RES, &c[i]);
        }
        torc_waitall();
#else
        get_nfc_task(&c[0]);
#endif

        unsigned int s = 0;
#ifdef VERBOSE
        printf("get_nfc:");
#endif
        for (int i = 0; i < torc_num_nodes(); ++i) {
            s += c[i];
#ifdef VERBOSE
            printf("+%d", c[i]);
#endif
        }
        g_nfeval = s;
#ifdef VERBOSE
        printf("=%d\n", s);
#endif
        t_nfeval += g_nfeval;
        return g_nfeval;
    }


    void reset_nfc_task()
    {
        l_nfeval = 0;
    }


    void reset_nfc() {
#ifdef _USE_TORC_
        for (int i = 0; i < torc_num_nodes(); ++i) {
            torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())reset_nfc_task, 0);
        }
        torc_waitall();
#else
        reset_nfc_task();
#endif
    }

    void gsl_rand_init(int seed){

        int local_workers = torc_i_num_workers();
        gsl_rng_env_setup();
        gsl_rng_t = gsl_rng_default;

        r = (gsl_rng **)malloc(local_workers*sizeof(gsl_rng *));


        for (int i = 0; i < local_workers; ++i) {
            r[i] = gsl_rng_alloc (T);
            //printf("...... %p \n", r[i] );
        }

        if (seed == 0) seed = time(0);

        for (int i = 0; i < local_workers; ++i) {
#if VERBOSE
            printf("node %d: initializing rng %d with seed %d\n", torc_node_id(), i, seed+i+local_workers*torc_node_id());
#endif
            gsl_rng_set(r[i], seed+i+local_workers*torc_node_id());
        }

        local_seed = (int *)malloc(local_workers*sizeof(int));
        for (int i = 0; i < local_workers; ++i) {
            local_seed[i] = seed+i+local_workers*torc_node_id();
        }
    }

    void spmd_gsl_rand_init(int seed)
    {
#ifdef _USE_TORC_
    	for (int i = 0; i < torc_num_nodes(); ++i) {
        	torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())call_gsl_rand_init, 0);
    	}
    	torc_waitall();
#else
		call_gsl_rand_init(seed);
#endif
    }

    void call_gsl_rand_init(int seed) {
#ifdef VERBOSE
        // printf("CALLING gsl_rand_init() on node %d\n", torc_node_id()); fflush(0);
#endif
        gsl_rand_init(seed);
    }


    void print_matrix(const char *name, double *x, int n) 
    {
        printf("\n%s =\n\n", name);
        for (int i = 0; i < n; ++i) printf("   %20.15lf\n", x[i]);
        printf("\n");
    }
    
    void print_matrix_2d(const char *name, double **x, int n1, int n2) 
    {
        printf("\n%s =\n\n", name);
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                printf("   %20.15lf", x[i][j]);
            }
            printf("\n");
        }
        printf("\n");

    }


    int mvnrnd(double *mean, double *sigma, double *out, int N) {

        gsl_vector_view mean_view 	= gsl_vector_view_array(mean, N);
        gsl_matrix_view sigma_view 	= gsl_matrix_view_array(sigma, N,N);
        gsl_vector_view out_view 	= gsl_vector_view_array(out, N);

        int me = torc_i_worker_id();

        gsl_matrix *L = gsl_matrix_alloc(N,N);
        gsl_matrix_memcpy( L, &sigma_view.matrix);
        gsl_linalg_cholesky_decomp( L );


        int res = gsl_ran_multivariate_gaussian( r[me], &mean_view.vector, L, &out_view.vector);

        return res;
    }


    // TODO: rename torc sth (DW)
    double uniformrand(double a, double b)
    {
        double res;

        int me = torc_i_worker_id();
        res = gsl_ran_flat(r[me], a, b);

        return res;
    }


    double uniform_pdf(double x, double *p){

        return gsl_ran_flat_pdf( x, p[0] , p[1] );

    }


    double uniform_log_pdf(double x, double *p){

        if( x>=p[0] && x<=p[1] )
            return -log(p[1]-p[0]);
        else
           return -INFINITY;	
    }


    double uniform_rnd( double *p ){
        
        double res;
        
        int me = torc_i_worker_id();
        res = gsl_ran_flat( r[me], p[0], p[1] );

        return res;
    }


    // TODO: rename torc sth (DW)
    void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[])
    {
        int me = torc_i_worker_id();
        gsl_ran_multinomial (r[me], K, N, q, nn);

        return;
    }


} //namespace
