#include <math.h>

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


    void inc_nfc() {
        pthread_mutex_lock(&feval_m);
        l_nfeval++;
        pthread_mutex_unlock(&feval_m);
    }


    void reset_nfc_task() {
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
