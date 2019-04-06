#ifndef TMCMC_TYPES_HPP
#define TMCMC_TYPES_HPP

#include <stdio.h>


#ifdef _USE_TORC_
#include <mpi.h>
#include <torc.h>
#else
#include <pthread.h>
#endif

namespace tmcmc
{

enum Method {
    Standard,
    Manifold
};

struct optim_options {
    int    MaxIter;             /* Max number of search iterations */
    double Tol;                 /* Tolerance for root finding */
    int    Display;             /* Print output */
    double Step;                /* Search stepsize */
    double LowerBound;          /* Lower bound for root finding (fmincon &
                                        fzerosearch)*/
    double UpperBound;          /* Upper bound for root finding (fmincon &
                                        fzerosearch)*/
    double Zdump;               /* Dump fzerosearch  stats to file */
};


struct manifold_options {
    double eps;                 /* manifold time step discretization */
    bool   adapt;               /* do COV adaptions? */
    double conf;                /* confidence region to fit bounds (or ebds) */
    double chi2;                /* chi square value calcualted from conf and dim */
    bool   use_ebds;           /* use extended bounds for cov adaption? */
    double pct_elb, pct_eub;    /* extended bounds in pct of domain widths */
    double *elbds, *eubds;      /* actual values for extended bounds ([PROBDIM]) */
};


typedef struct data_t {
    data_t(const char * fname = "tmcmc.par");
    ~data_t();
    int Nth;                        /* PROBDIM */
    int MaxStages;                  /* Max number of tmcmc generations */
    int PopSize;                    /* Number of samples */

    double *lowerbound;             /*[PROBDIM];*/
    double *upperbound;             /*[PROBDIM];*/

    double *compositeprior_distr;   /*[PROBDIM]*/

    double *prior_mu;
    double *prior_sigma;

    int    auxil_size;
    double *auxil_data;

    int MinChainLength;         /* MinChainLength >= 0: setting MinChainLength */
    int MaxChainLength;         /* MaxChainLength >= 0: splitting long chains */

    double lb, ub;              /* generic lower and upper bound*/

    double TolCOV;              /* Target coefficient of variation of weights */
    double MinStep;             /* Min update of rho */
    double beta2;               /* Covariance scaling parameter */
    long   seed;
    int    burn_in;             /* Number of burn in iterations */

    optim_options options;      /* Optimization options (see above) */
    manifold_options moptions;  /* Options for mTMCMC (see above) */

    int prior_type;             /* 0: uniform, 1: gaussian, 3: composite */
    int load_from_file;

    int icdump;
    int ifdump;

    int *Num;               /*[MAXGENS];*/

    int    use_proposal_cma;
    double **init_mean;     /* [DATANUM][PROBDIM] */

    double **local_cov;     /* [DATANUM][PROBDIM*PROBDIM] */
    bool use_local_cov;
    double local_scale;

    int stealing;
    int restart;
} data_t;


typedef struct runinfo_t {
    ~runinfo_t();
    static void init(runinfo_t& runinfo, int nth, int maxstages);
    static void save(const runinfo_t& runinfo, int nth, int maxstages,
                     const char * fname = "runinfo.txt");
    static bool load(runinfo_t& runinfo, int nth, int maxstages,
                     const char * fname = "runinfo.txt");
    int    Gen;
    double *CoefVar;        /*[MAXGENS];*/
    double *p;              /*[MAXGENS];*/
    int    *currentuniques; /*[MAXGENS];*/
    double *logselections;  /*[MAXGENS];*/
    double *acceptance;     /*[MAXGENS];*/
    double **SS;            /*[PROBDIM][PROBDIM];*/
    double **meantheta;     /*[MAXGENS][PROBDIM]*/

    int *outside;           /*[MAXGENS]*/
    int *corrections;       /*[MAXGENS]*/
    int *failedcorrections; /*[MAXGENS]*/

} runinfo_t;


typedef struct cgdbp_s {
    double *point;      /*[PROBDIM];*/
    double F;
    double prior;

    int counter;        /* not used (?)*/
    int nsel;           /* for selection of leaders only*/
    int queue;          /* for submission of leaders only*/
    int surrogate;      //TODO: used? (DW)
    double error;       //TODO: used? (DW)

    int error_flg; 
    int posdef;         //TODO: can we combine this with error_flg? (DW)
    double *gradient;   /*[PROBDIM]*/
    double *cov;        /*[PROBDIM*PROBDIM]*/
    double *evec;       /*[PROBDIM*PROBDIM]*/
    double *eval;       /*[PROBDIM]*/
} cgdbp_t;


typedef struct cgdb_s {
    size_t  entries;    
    cgdbp_t *entry;     /*[MAX_DB_ENTRIES];*/
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
    int idx  = -1;
    int nsel = 0;
} sort_t;


} //namespace tmcmc

#endif//TMCMC_TYPES_HPP
