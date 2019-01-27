#ifndef TMCMC_TYPES_HPP
#define TMCMC_TYPES_HPP

#include <stdio.h>


#ifdef _USE_TORC_
    #include <mpi.h>
    #include <torc.h>
#else
    #include <pthread.h>
#endif 

namespace tmcmc {

    struct optim_options {
        int    MaxIter;
        double Tol;
        int    Display;
        double Step;
    };


    typedef struct data_t {
        data_t(const char * fname = "tmcmc.par");
        int Nth;       
        int MaxStages;
        int PopSize; 

        double *lowerbound;    /*[PROBDIM];*/
        double *upperbound;    /*[PROBDIM];*/

        double *compositeprior_distr; /*[PROBDIM]*/

        double *prior_mu;
        double *prior_sigma;

        int    auxil_size;
        double *auxil_data;

        int MinChainLength;
        int MaxChainLength;

        double lb, ub;        /*generic lower and upper bound*/

        double TolCOV;
        double bbeta;
        long   seed;
        int    burn_in;

        optim_options options;

        int prior_type;     /* 0: uniform, 1: gaussian, 3: composite */
        int load_from_file;

        int icdump;
        int ifdump;

        int *Num;        /*[MAXGENS];*/
        int  LastNum;

        int    use_proposal_cma;
        double **init_mean;    /* [DATANUM][PROBDIM] */

        double **local_cov;    /* [DATANUM][PROBDIM*PROBDIM] */
        int    use_local_cov;
        double local_scale;

        int stealing;
        int restart;
    } data_t;


    typedef struct runinfo_t {
        static void init(runinfo_t& runinfo, int nth, int maxstages);
        static void save(const runinfo_t& runinfo, int nth, int maxstages, 
                        const char * fname = "runinfo.txt");
        static bool load(runinfo_t& runinfo, int nth, int maxstages, 
                        const char * fname = "runinfo.txt");
        int    Gen;
        double *CoefVar;        /*[MAXGENS];*/
        double *p;            /*[MAXGENS];        // cluster-wide*/
        int    *currentuniques;    /*[MAXGENS];*/
        double *logselections;        /*[MAXGENS];*/
        double *acceptance;        /*[MAXGENS];*/
        double **SS;            /*[PROBDIM][PROBDIM];    // cluster-wide*/
        double **meantheta;         /*[MAXGENS][PROBDIM]*/
    } runinfo_t;


    typedef struct cgdbp_s {
        double *point; /*[PROBDIM];*/
        double F;
        double prior;

        int counter;    /* not used (?)*/
        int nsel;    /* for selection of leaders only*/
        int queue;    /* for submission of leaders only*/
        int surrogate;
        double error;
    } cgdbp_t;


    typedef struct cgdb_s {
        int     entries;
        cgdbp_t *entry; /*[MAX_DB_ENTRIES];*/
        pthread_mutex_t m;
    } cgdb_t;


    typedef struct dbp_s {
        double *point; /*[PROBDIM];*/
        double F;
        int    nG;
        double G[64];    /* maxG*/
        int    surrogate;
    } dbp_t;


    typedef struct db_s {
        int   entries;
        dbp_t *entry; /*[MAX_DB_ENTRIES];*/        /* */
        pthread_mutex_t m;
    } db_t;


    typedef struct resdbp_s {
        double *point;    /*[EXPERIMENTAL_RESULTS+1]; // +1 for the result (F)*/
        double F;
        int counter;    /* not used (?)*/
        int nsel;    /* for selection of leaders only*/
    } resdbp_t;


    typedef struct resdb_s {
        int      entries;
        resdbp_t *entry; /*[MAX_DB_ENTRIES];*/
        pthread_mutex_t m;
    } resdb_t;

    
    typedef struct fparam_s {
        const double *fj;
        int           fn;
        double        pj;
        double        tol;
    } fparam_t;

    
    typedef struct sort_s {
        int idx;
        int nsel;
        double F;
    } sort_t;


} //namespace tmcmc

#endif//TMCMC_TYPES_HPP
