#include <cstring>
#include <stdlib.h>
#include <limits>

#include <gsl/gsl_cdf.h>

#include "tmcmc_types.hpp"

namespace tmcmc
{

data_t::data_t(const char * fname)
{
    Nth       = -1;
    MaxStages = -1;
    PopSize   = -1;

    MinChainLength = 0;
    MaxChainLength = 0;

    TolCOV  = -1;
    MinStep = 1e-6;
    bbeta   = -1;
    seed    = 280675;
    burn_in = -1;

    options.MaxIter    = -1;
    options.Tol        = -1;
    options.Display    = 0;
    options.Step       = -1;
    options.LowerBound = 0.0;
    options.UpperBound = 1.0;
    options.Zdump      = 0;


    moptions.use_ebds = false;
    moptions.conf     = 0.68;
    moptions.pct_elb  = 0.0;
    moptions.pct_eub  = 0.0;
    moptions.eps      = 0.1;

    load_from_file = 0;

    icdump = 1;
    ifdump = 0;

    stealing = 0;
    restart  = 0;

    int idx = 0;

    FILE *f = fopen(fname, "r");
    if (f == NULL) {
        printf("\nThe input file 'tmcmc.par' is missing. Exit...\n");
        exit(EXIT_FAILURE);
    }

    char line[256];
    // TODO: can this initialization be aligned with CMAES? (DW)
    int line_no = 0;
    while (fgets(line, 256, f)!= NULL) {
        line_no++;
        if ((line[0] == '#')||(strlen(line)==0)) {
            continue;
        }

        if (strstr(line, "Nth")) {
            sscanf(line, "%*s %d", &Nth);
            lowerbound = new double[Nth];
            upperbound = new double[Nth];
        } else if (strstr(line, "MaxStages")) {
            sscanf(line, "%*s %d", &MaxStages);
            printf("setting MaxStages = %d\n", MaxStages);
        } else if (strstr(line, "PopSize")) {
            sscanf(line, "%*s %d", &PopSize);
            printf("setting PopSize = %d\n", PopSize);
        } else if (strstr(line, "TolCOV")) {
            sscanf(line, "%*s %lf", &TolCOV);
            printf("setting TolCOV = %lf\n", TolCOV);
        } else if (strstr(line, "MinStep")) {
            sscanf(line, "%*s %lf", &MinStep);
            printf("setting Minstep = %lf\n", MinStep);
        } else if (strstr(line, "bbeta")) {
            sscanf(line, "%*s %lf", &bbeta);
            printf("setting bbeta = %lf\n", bbeta);
        } else if (strstr(line, "seed")) {
            sscanf(line, "%*s %ld", &seed);
            printf("setting seed = %ld\n", seed);
        } else if (strstr(line, "burn_in")) {
            sscanf(line, "%*s %d", &burn_in);
            printf("setting burn_in = %d\n", burn_in);
        } else if (strstr(line, "opt.MaxIter")) {
            sscanf(line, "%*s %d", &options.MaxIter);
            printf("setting options.MaxIter = %d\n", options.MaxIter);
        } else if (strstr(line, "opt.Tol")) {
            sscanf(line, "%*s %lf", &options.Tol);
            printf("setting options.Tol = %lf\n", options.Tol);
        } else if (strstr(line, "opt.Display")) {
            sscanf(line, "%*s %d", &options.Display);
            printf("setting options.Display = %d\n", options.Display);
        } else if (strstr(line, "opt.Step")) {
            sscanf(line, "%*s %lf", &options.Step);
            printf("setting options.Step = %f\n", options.Step);
        } else if (strstr(line, "opt.Lowerbound")) {
            sscanf(line, "%*s %lf", &options.LowerBound);
            printf("setting options.Lowerbound = %f\n", options.LowerBound);
        } else if (strstr(line, "opt.Upperbound")) {
            sscanf(line, "%*s %lf", &options.UpperBound);
            printf("setting options.Upperbound = %f\n", options.UpperBound);
        } else if (strstr(line, "icdump")) {
            sscanf(line, "%*s %d", &icdump);
            printf("setting icdump = %d\n", icdump);
        } else if (strstr(line, "ifdump")) {
            sscanf(line, "%*s %d", &ifdump);
            printf("setting ifdump = %d\n", ifdump);
        } else if (strstr(line, "use_local_cov")) {
            sscanf(line, "%*s %d", &use_local_cov);
            printf("setting use_local_cov = %d\n", use_local_cov);
        } else if (strstr(line, "stealing")) {
            sscanf(line, "%*s %d", &stealing);
            printf("setting stealing = %d\n", stealing);
        } else if (strstr(line, "restart")) {
            sscanf(line, "%*s %d", &restart);
            printf("setting restart = %d\n", restart);
        } else if (strstr(line, "Bound")) {
            sscanf(line, "%*s %lf %lf", &lowerbound[idx], &upperbound[idx]);
            printf("setting bounds = %lf %lf\n", lowerbound[idx], upperbound[idx]);
            idx++;
        } else if (strstr(line, "MaxChainLength")) {
            sscanf(line, "%*s %d", &MaxChainLength);
            printf("setting MaxChainLength = %d\n", MaxChainLength);
        } else if (strstr(line, "MinChainLength")) {
            sscanf(line, "%*s %d", &MinChainLength);
            printf("setting MinChainLength = %d\n", MinChainLength);
        } else if (strstr(line, "moptions.eps")) {
            sscanf(line, "%*s %lf", &moptions.eps);
            printf("setting manifold options epsilon = %lf\n", moptions.eps);
        } else if (strstr(line, "moptions.use_ebds")) {
            sscanf(line, "%*s %d", &moptions.use_ebds);
            printf("setting manifold options use_ebds = %d\n", moptions.use_ebds);
        } else if (strstr(line, "moptions.conf")) {
            sscanf(line, "%*s %lf", &moptions.conf);
            printf("setting manifold options conf = %lf\n", moptions.conf);
        }  else if (strstr(line, "moptions.pct_elb")) {
            sscanf(line, "%*s %lf", &moptions.pct_elb);
            printf("setting manifold options pct_elb = %lf\n", moptions.pct_elb);
        } else if (strstr(line, "moptions.pct_eub")) {
            sscanf(line, "%*s %lf", &moptions.pct_eub);
            printf("setting manifold options.pct_eub = %lf\n", moptions.pct_eub);
        }

    }

    fclose(f);

    if (idx != Nth) {
        printf("\nNumber of Bounds does not match Nth. Exit...\n");
        exit(1);
    }

    Num = new int[MaxStages];
    for (int i = 0; i < MaxStages; ++i) {
        Num[i] = PopSize;
    }

    double *LCmem = new double[PopSize*Nth*Nth];
    local_cov     = new double*[PopSize];
    for (int pos=0; pos < PopSize; ++pos) {
        local_cov[pos] = LCmem + pos*Nth*Nth;
        for (int i=0; i<Nth; ++i)
            local_cov[pos][i*Nth+i] = 1;
    }

    moptions.chi2 = gsl_cdf_chisq_Pinv(moptions.conf,Nth);
    moptions.elbds = new double[Nth];
    moptions.eubds = new double[Nth];
    for(int i = 0; i < Nth; ++i) {
        moptions.elbds[i] = lowerbound[i];
        moptions.eubds[i] = upperbound[i];
        if (moptions.use_ebds) {
            double width = upperbound[i] - lowerbound[i];
            moptions.elbds[i] -= width;
            moptions.eubds[i] += width;
        }
    }

}

data_t::~data_t()
{
    delete [] lowerbound;
    delete [] upperbound;
    delete [] moptions.elbds;
    delete [] moptions.eubds;
    delete [] Num;
    delete [] local_cov[0];
    delete [] local_cov;
}

void runinfo_t::init(runinfo_t& runinfo, int nth, int maxstages)
{

    runinfo.CoefVar        = new double[maxstages+1];
    runinfo.p              = new double[maxstages+1];
    runinfo.currentuniques = new int[maxstages];
    runinfo.logselections  = new double[maxstages];
    runinfo.acceptance     = new double[maxstages];

    runinfo.Gen = 0;
    runinfo.CoefVar[0] = std::numeric_limits<double>::infinity();

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


runinfo_t::~runinfo_t()
{
    delete [] CoefVar;
    delete [] p;
    delete [] currentuniques;
    delete [] logselections;
    delete [] acceptance;
    delete [] SS[0];
    delete [] SS;
    //TODO: delete content of meantheta (DW);
    delete [] meantheta;
}


void runinfo_t::save(const runinfo_t& runinfo, int nth, int maxstages, const char * fname)
{

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


bool runinfo_t::load(runinfo_t& runinfo, int nth, int maxstages, const char * fname)
{
    char header[256];

    /* allocate and initialize runinfo */
    FILE *f = fopen("runinfo.txt", "r");
    if (f == NULL) return false;

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

    return true;
}


}//namespace tmcmc
