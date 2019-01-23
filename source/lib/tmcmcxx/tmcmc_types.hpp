#ifndef TMCMC_TYPES
#define TMCMC_TYPES

#include <stdio.h>
#include <string>

#include "tmcmc_utils.hpp"

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


data_t::data_t(const char * fname) {
    Nth       = -1;
    MaxStages = -1;
    PopSize   = -1;

    MinChainLength = 0;
    MaxChainLength = 0;

    TolCOV  = -1;
    bbeta   = -1;
    seed    = 280675;
    burn_in = -1;

    options.MaxIter = -1;
    options.Tol     = -1;
    options.Display = 0;
    options.Step    = -1;

    load_from_file = 0;

    icdump = 1;
    ifdump = 0;

    stealing = 0;
    restart  = 0;

    FILE *f = fopen(fname, "r");
    if (f == NULL) {
		printf("\nThe input file 'tmcmc.par' is missing. Exit...\n");
		exit(1);
    }

    char line[256];
    // TODO: can this initialization be aligned with CMAES?
    int line_no = 0;
    while (fgets(line, 256, f)!= NULL) {
        line_no++;
        if ((line[0] == '#')||(strlen(line)==0)) {
            continue;
        }

        if (strstr(line, "Nth")) {
            sscanf(line, "%*s %d", &Nth);
        }
        else if (strstr(line, "MaxStages")) {
            sscanf(line, "%*s %d", &MaxStages);
        }
        else if (strstr(line, "PopSize")) {
            sscanf(line, "%*s %d", &PopSize);
        }
        else if (strstr(line, "TolCOV")) {
            sscanf(line, "%*s %lf", &TolCOV);
        }
        else if (strstr(line, "bbeta")) {
            sscanf(line, "%*s %lf", &bbeta);
        }
        else if (strstr(line, "seed")) {
            sscanf(line, "%*s %ld", &seed);
        }
        else if (strstr(line, "burn_in")) {
            sscanf(line, "%*s %d", &burn_in);
        }
        else if (strstr(line, "opt.MaxIter")) {
            sscanf(line, "%*s %d", &options.MaxIter);
        }
        else if (strstr(line, "opt.Tol")) {
            sscanf(line, "%*s %lf", &options.Tol);
        }
        else if (strstr(line, "opt.Display")) {
            sscanf(line, "%*s %d", &options.Display);
        }
        else if (strstr(line, "opt.Step")) {
            sscanf(line, "%*s %lf", &options.Step);
            printf("setting step = %f\n", options.Step);
        }
        else if (strstr(line, "icdump")) {
            sscanf(line, "%*s %d", &icdump);
        }
        else if (strstr(line, "ifdump")) {
            sscanf(line, "%*s %d", &ifdump);
        }
        else if (strstr(line, "use_local_cov")) {
            sscanf(line, "%*s %d", &use_local_cov);
        }
        else if (strstr(line, "stealing")) {
            sscanf(line, "%*s %d", &stealing);
        }
        else if (strstr(line, "restart")) {
            sscanf(line, "%*s %d", &restart);
        }
    }

	fclose(f);


	//XXX add: check if all parameters are ok
    lowerbound = new double[Nth];
    upperbound = new double[Nth]; 
    for(int i=0; i<Nth; ++i){
 		lowerbound[i] = std::numeric_limits<double>::min();
    	upperbound[i] = std::numeric_limits<double>::max();
	}

    Num = new int[MaxStages];
    for (int i = 0; i < MaxStages; ++i){
        Num[i] = PopSize;
    }
    
    LastNum = PopSize;

	double *LCmem = new double[PopSize*Nth*Nth];
    local_cov     = new double*[PopSize];
    for (int pos=0; pos < PopSize; ++pos){
        local_cov[pos] = LCmem + pos*Nth*Nth;
        for (int i=0; i<Nth; ++i)
            local_cov[pos][i*Nth+i] = 1;
    }

}


typedef struct runinfo_t {
    static void init(runinfo_t& runinfo, int nth, int maxstages);
    static void save(const runinfo_t& runinfo, int nth, int maxstages, 
                    const char * fname = "runinfo.txt");
    static int load(runinfo_t& runinfo, int nth, int maxstages, 
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

void runinfo_t::init(runinfo_t& runinfo, int nth, int maxstages) {
    
    runinfo.Gen = 0;

    runinfo.CoefVar        = new double[maxstages+1];
    runinfo.CoefVar[0]     = 10;
    runinfo.p              = new double[maxstages+1];
    runinfo.currentuniques = new int[maxstages];
    runinfo.logselections  = new double[maxstages];
    runinfo.acceptance     = new double[maxstages];
    
    // TODO: is this leaking? (DW)
    double *SSmem = new double[nth*nth];
    runinfo.SS = new double*[nth];
    for(int i = 0; i< nth; ++i) {
        runinfo.SS[i] = SSmem + i*nth;
    }

    runinfo.meantheta = new double*[maxstages+1];
    for(int i = 0; i < maxstages+1; ++i) {
        runinfo.meantheta[i] = new double[nth];
    }
}

void runinfo_t::save(const runinfo_t& runinfo, int nth, int maxstages, const char * fname) {

    FILE *f = fopen(fname, "w");

    fprintf(f, "Gen=\n");
    fprintf(f, "%d\n", runinfo.Gen);

    fprintf(f, "CoefVar[%d]=\n", maxstages);
    for (int i = 0; i < maxstages; ++i) fprintf(f, "%.16lf\n", runinfo.CoefVar[i]);

    fprintf(f, "p[%d]=\n", maxstages);
    for (int i = 0; i < maxstages; ++i) fprintf(f, "%.16lf\n", runinfo.p[i]);

    fprintf(f, "currentuniques[%d]=\n", maxstages);
    for (int i = 0; i < maxstages; ++i) fprintf(f, "%d\n", runinfo.currentuniques[i]);

    fprintf(f, "logselection[%d]=\n", maxstages);
    for (int i = 0; i < maxstages; ++i) fprintf(f, "%.16lf\n", runinfo.logselections[i]);

    fprintf(f, "acceptance[%d]=\n", maxstages);
    for (int i = 0; i < maxstages; ++i) fprintf(f, "%.16lf\n", runinfo.acceptance[i]);

    fprintf(f, "SS[%d][%d]=\n", nth, nth);
    for (int i = 0; i < nth; ++i)
        for (int j = 0; j < nth; ++j)
            fprintf(f, "%.16lf\n", runinfo.SS[i][j]);

    fprintf(f, "meantheta[%d][%d]\n", maxstages, nth);
    for (int i = 0; i < maxstages; ++i)
        for (int j = 0; j < nth; ++j)
            fprintf(f, "%.16lf\n", runinfo.meantheta[i][j]);

    fclose(f);
}


int runinfo_t::load(runinfo_t& runinfo, int nth, int maxstages, const char * fname) {
    char header[256];

    /* allocate and initialize runinfo */
    FILE *f = fopen("runinfo.txt", "r");
    if (f == NULL) return 1;

    fscanf(f, "%s", header);
    fscanf(f, "%d", &runinfo.Gen);

    fscanf(f, "%s", header);
    for (int i = 0; i < maxstages; ++i) fscanf(f, "%lf\n", &runinfo.CoefVar[i]);

    fscanf(f, "%s", header);
    for (int i = 0; i < maxstages; ++i) fscanf(f, "%lf\n", &runinfo.p[i]);

    fscanf(f, "%s", header);
    for (int i = 0; i < maxstages; ++i) fscanf(f, "%d\n", &runinfo.currentuniques[i]);

    fscanf(f, "%s", header);
    for (int i = 0; i < maxstages; ++i) fscanf(f, "%lf\n", &runinfo.logselections[i]);

    fscanf(f, "%s", header);
    for (int i = 0; i < maxstages; ++i) fscanf(f, "%lf\n", &runinfo.acceptance[i]);

    fscanf(f, "%s", header);
    for (int i = 0; i < nth; ++i)
        for (int j = 0; j < nth; ++j)
            fscanf(f, "%lf\n", &runinfo.SS[i][j]);

    fscanf(f, "%s", header);
    for (int i = 0; i < maxstages; ++i)
        for (int j = 0; j < nth; ++j)
            fscanf(f, "%lf\n", &runinfo.meantheta[i][j]);

    fclose(f);

    return 0;
}


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


#endif
