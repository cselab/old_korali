#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

#include "tmcmc_utils.hpp"

namespace tmcmc {

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


    // TODO: rename torc sth (DW)
    void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[])
    {
        int me = torc_i_worker_id();
        gsl_ran_multinomial (r[me], K, N, q, nn);

        return;
    }



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

} //namespace
