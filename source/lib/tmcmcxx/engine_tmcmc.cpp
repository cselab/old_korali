#include <cmath>
#include <numeric>

#include <iostream> //cin test remove

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>

#include "fitfun.hpp"
#include "tmcmc_obj_fmin.hpp"
#include "engine_tmcmc.hpp"

using namespace priors;

namespace tmcmc {

    TmcmcEngine::TmcmcEngine() : data(data_t()), 
                                 nchains(data.Num[0]),
                                 out_tparam(new double[data.PopSize]),
                                 leaders(new cgdbp_t[data.PopSize]),
                                 prior("priors.par") {
        // TODO (DW)
        for (int i = 0; i< data.PopSize; ++i) 
            leaders[i].point = new double[data.Nth];

        curres_db.entries = 0;
    }

    TmcmcEngine::~TmcmcEngine() {
        // TODO (DW)
    }

    void TmcmcEngine::run() {
        init(); 

        if (data.MaxStages == 1) {
            printf("Maxstages == 1, nothing to do\n");
            return; 
        } else if (runinfo.p[runinfo.Gen] == 1.0) {
            printf("p == 1 from previous run, nothing more to do\n");
            return;
        }

        nchains = prepare_newgen(nchains, leaders);
        spmd_update_runinfo();
        if (data.options.Display) print_runinfo();


        while(++runinfo.Gen < data.MaxStages && runinfo.p[runinfo.Gen] < 1.0) {
            evalGen();
        }

        print_matrix((char *)"runinfo.p", runinfo.p, runinfo.Gen+1);
        print_matrix((char *)"runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
        print_matrix_i((char *)"runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
        print_matrix((char *)"runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
        print_matrix((char *)"runinfo.logselection", runinfo.logselections, runinfo.Gen+1);

        double logEvidence[1];
        logEvidence[0] = compute_sum(runinfo.logselections, runinfo.Gen+1);
        print_matrix((char *)"logEvidence", logEvidence, 1);

        FILE *fp;
        fp = fopen("log_evidence.txt", "w");
        fprintf(fp, "%lf\n", logEvidence[0]);
        fclose(fp);

        runinfo_t::save(runinfo, data.Nth, data.MaxStages);

        if (data.icdump){
            printf("lastgen = %d\n", runinfo.Gen);
            char cmd[256];
            sprintf(cmd, "cp curgen_db_%03d.txt final.txt", runinfo.Gen);
            system(cmd);
        }

        fitfun::fitfun_finalize();
        printf("total function calls = %d\n", get_tfc());
#if defined(_USE_TORC_)
        torc_finalize();
#endif
        return;
    }


    void TmcmcEngine::init() {
        
        gt0 = t0 = torc_gettime();

        init_curgen_db();
        init_curres_db();
        init_full_db();

        runinfo_t::init(runinfo, data.Nth, data.MaxStages);

        if (data.options.Display) print_runinfo();

        spmd_gsl_rand_init(data.seed);

        prior.print();
        
        curgen_db.entries = 0;
        
        bool loaded = false;
        if (data.restart) loaded = load_data();
        if (loaded == false) {
            sample_from_prior();
            if( data.icdump ) dump_curgen_db();
            if( data.ifdump ) dump_full_db();
        }
        
        runinfo_t::save(runinfo, data.Nth, data.MaxStages);
        if (data.restart) check_for_exit();
        
        printf("----------------------------------------------------------------\n");
        
        return;
    }

    void TmcmcEngine::evalGen() {
    
        int nsteps;
        gt0 = torc_gettime();

#ifdef _USE_OPENMP_
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
#endif
                int winfo[4];
                double in_tparam[data.Nth];
                double init_mean[data.Nth];
                double chain_cov[data.Nth*data.Nth];


                for (int i = 0; i < nchains; ++i) {
                    winfo[0] = runinfo.Gen;
                    winfo[1] = i;
                    winfo[2] = -1;    /* not used */
                    winfo[3] = -1;    /* not used */
                
                    for(int d = 0; d < data.Nth; ++d) 
                        in_tparam[d] = leaders[i].point[d];

                    nsteps = leaders[i].nsel;
                    
                    if (data.use_local_cov){
                        for (int d = 0; d < data.Nth*data.Nth; ++d)
                            chain_cov[d] = data.local_cov[i][d];

                        for (int d = 0; d < data.Nth; ++d) {
                            if (data.use_proposal_cma)
                                init_mean[d] = data.init_mean[i][d];
                            else
                                init_mean[d] = leaders[i].point[d];
                        }
                    } else {
                        for (int d = 0; d < data.Nth; ++d)
                            for (int e = 0; e < data.Nth; ++e)
                                chain_cov[d*data.Nth+e]=
                                    data.bbeta*runinfo.SS[d][e];

                        for (int d = 0; d < data.Nth; ++d)
                            init_mean[d] = in_tparam[d];
                    }

                    out_tparam[i] = leaders[i].F;    /* loglik_leader...*/


#ifdef _USE_TORC_
                    torc_create(leaders[i].queue, (void (*)())chaintask, 7,
                                data.Nth, MPI_DOUBLE, CALL_BY_COP,
                                1, MPI_INT, CALL_BY_COP,
                                1, MPI_INT, CALL_BY_COP,
                                1, MPI_DOUBLE, CALL_BY_REF,
                                4, MPI_INT, CALL_BY_COP,
                                data.Nth, MPI_DOUBLE, CALL_BY_COP,
                                data.Nth*data.Nth, MPI_DOUBLE, CALL_BY_COP,
                                in_tparam, &data.Nth, &nsteps, &out_tparam[i], winfo,
                                init_mean, chain_cov);
#else
#ifdef _USE_OPENMP_
                    #pragma omp task shared(data, out_tparam) firstprivate(i, nsteps, in_tparam, winfo, init_mean, chain_cov)
#endif
                    chaintask( in_tparam, &nsteps, &out_tparam[i], winfo, init_mean, chain_cov );
#endif
                }

#if defined(_USE_OPENMP_)
            }// single
        }// parallel
#endif

#if defined(_USE_TORC_)
        if (data.stealing) torc_enable_stealing();

        torc_waitall();
        if (data.stealing) torc_disable_stealing();
#endif

        gt1 = torc_gettime();
        int g_nfeval = get_nfc();
        printf("evalGen: Generation %d: total elapsed time = %lf secs, \
                    generation elapsed time = %lf secs for function calls = %d\n", 
                    runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
        
        reset_nfc();
		
        if (data.icdump) dump_curgen_db();
        if (data.ifdump) dump_full_db();

        runinfo_t::save(runinfo, data.Nth, data.MaxStages);
 
        if (data.restart) check_for_exit();

        curres_db.entries = 0;
        nchains = prepare_newgen(nchains, leaders);

        spmd_update_runinfo(); 

        if(data.options.Display) print_runinfo();

        printf("Acceptance rate   :  %lf \n",runinfo.acceptance[runinfo.Gen]) ;
        printf("Annealing exponent:  %lf \n",runinfo.p[runinfo.Gen]) ;
        printf("----------------------------------------------------------------\n");
        return; 
    }
    
    bool TmcmcEngine::load_data() {
        bool res = runinfo_t::load(runinfo, data.Nth, data.MaxStages);
    	if (res == 0) {
            load_curgen_db();
            printf("nchains = %d\n", data.Num[0]);
    	}
        return res;
    }

    void TmcmcEngine::sample_from_prior() {

#ifdef _USE_OPENMP_
        #pragma omp parallel
        {
            printf("Hello from thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
            #pragma omp for
#endif
            int winfo[4];
            for (int i = 0; i < nchains; ++i){
                winfo[0] = runinfo.Gen;
                winfo[1] = i;
                winfo[2] = -1;
                winfo[3] = -1;

                //store draw from prior
                double in_tparam[data.Nth]; 
                for (int d = 0; d < data.Nth; ++d)
                    in_tparam[d] = prior.rand(d); 

#ifdef _USE_TORC_
                torc_create(-1, (void (*)())initchaintask, 4,
                            data.Nth, MPI_DOUBLE, CALL_BY_COP,
                            1, MPI_INT, CALL_BY_COP,
                            1, MPI_DOUBLE, CALL_BY_RES,
                            4, MPI_INT, CALL_BY_COP,
                            in_tparam, &data.Nth, &out_tparam[i], winfo);
#else
#ifdef _USE_OPENMP_
                //#pragma omp task shared(data, out_tparam) firstprivate(i, in_tparam, winfo)
                #pragma omp task firstprivate(i, winfo, in_tparam) shared(data, out_tparam)
#endif
                {
                    initchaintask(in_tparam, &out_tparam[i], winfo);
                }
#endif
            }

#if defined(_USE_OPENMP_)
        } 
#endif


#if defined(_USE_TORC_)
        if (data.stealing) torc_enable_stealing();

        torc_waitall();

        if (data.stealing) torc_disable_stealing();
#endif

        gt1 = torc_gettime();
        
        int g_nfeval = get_nfc();
        
        printf("sample_from_prior: Generation %d: total elapsed time = %lf sec,\
                generation elapsed time = %lf secs for function calls = %d\n", 
                runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
        
        reset_nfc();
	       
        return;
    }


    void TmcmcEngine::init_full_db() {
        pthread_mutex_init(&full_db.m, NULL);
        full_db.entries = 0;
        full_db.entry   = new dbp_t[data.MaxStages*data.PopSize];
    }


    void TmcmcEngine::update_full_db(double point[], double F, double *G, int n, int surrogate) {
        pthread_mutex_lock(&full_db.m);
        int pos = full_db.entries;
        full_db.entries++;
        pthread_mutex_unlock(&full_db.m);

        if (full_db.entry[pos].point == NULL) 
            full_db.entry[pos].point = new double[data.Nth] ;
        
        for (int i = 0; i < data.Nth; ++i) full_db.entry[pos].point[i] = point[i];
        full_db.entry[pos].F = F;
        full_db.entry[pos].nG = n;

        for (int i = 0; i < n; ++i) full_db.entry[pos].G[i] = G[i];
        full_db.entry[pos].surrogate = surrogate;
    }


    void TmcmcEngine::torc_update_full_db(double point[], double F, double *G, int n, int surrogate) {
        if (torc_node_id() == 0) {
            update_full_db(point, F, G, n, surrogate);
            return;
        }

#ifdef _USE_TORC_
            if (n == 0)
                torc_create_direct( 0, (void (*)())torc_update_full_db_task, 3,
                                        /* message to the database manager (separate process?) or direct execution by server thread */
                                    data.Nth, MPI_DOUBLE, CALL_BY_VAL,
                                    1, MPI_DOUBLE, CALL_BY_COP,
                                    1, MPI_INT, CALL_BY_COP,
                                    point, &F, &surrogate);
            else
                torc_create_direct( 0, (void (*)())torc_update_full_db_task_p5, 5,
                                        data.Nth, MPI_DOUBLE, CALL_BY_VAL,
                                        1, MPI_DOUBLE, CALL_BY_COP,
                                        n, MPI_DOUBLE, CALL_BY_COP,    /* xxx: for CALL_BY_VAL: in the full-version of the library, with n=1 we had segv */
                                        1, MPI_INT, CALL_BY_COP,
                                        1, MPI_INT, CALL_BY_COP,
                                        point, &F, G, &n, &surrogate);

                torc_waitall3();
#endif
    }


    void TmcmcEngine::print_full_db() {
        int i;
        printf("=======\n");
        printf("FULL_DB\n");
        for (int pos = 0; pos < full_db.entries; ++pos) {
            printf("ENTRY %d: POINT(%20.16lf,%20.16lf) F=%20.16lf SG=%d\n",
                    pos, full_db.entry[pos].point[0], full_db.entry[pos].point[1],	/* extend it*/
                    full_db.entry[pos].F, full_db.entry[pos].surrogate);
            printf("\tG=[");
            for (i = 0; i < full_db.entry[pos].nG-1; ++i) printf("%20.16lf,", full_db.entry[pos].G[i]);
            printf("%20.16lf]\n", full_db.entry[pos].G[i]);
        }
        printf("=======\n");
    }


    void TmcmcEngine::dump_full_db() {
        FILE *fp;
        char fname[256];

        sprintf(fname, "full_db_%03d.txt", runinfo.Gen);
        fp = fopen(fname, "w");
        for (int pos = 0; pos < full_db.entries; pos++) {
            for (int i = 0; i < data.Nth; ++i) {
                fprintf(fp, "%20.16lf ", full_db.entry[pos].point[i]);
            }
            fprintf(fp, "%20.16lf\n", full_db.entry[pos].F);
        }
        fclose(fp);
    }


    void TmcmcEngine::init_curgen_db() {
        pthread_mutex_init(&curgen_db.m, NULL);
        curgen_db.entries = 0;
        curgen_db.entry   = new cgdbp_t[(data.MinChainLength+1)*data.PopSize];
    }


    void TmcmcEngine::update_curgen_db(double point[], double F, double prior) {
        pthread_mutex_lock(&curgen_db.m);
        int pos = curgen_db.entries;
        curgen_db.entries++;
        pthread_mutex_unlock(&curgen_db.m);

        if (curgen_db.entry[pos].point == NULL) curgen_db.entry[pos].point = (double *)malloc(data.Nth*sizeof(double));

        for (int i = 0; i < data.Nth; ++i) curgen_db.entry[pos].point[i] = point[i];
        curgen_db.entry[pos].F = F;
        curgen_db.entry[pos].prior = prior;
    }


    void TmcmcEngine::torc_update_curgen_db_task(double point[], double *pF, double *pprior) {
        double F = *pF;
        double prior = *pprior;
        update_curgen_db(point, F, prior);
    }


    void TmcmcEngine::torc_update_curgen_db(double point[], double F, double prior) {
        int me = torc_node_id();

        if (me == 0) {
            update_curgen_db(point, F,prior);
            return;
        }

#ifdef _USE_TORC_
            torc_create_direct(0, (void (*)())torc_update_curgen_db_task, 3,
                            /* message to the database manager (separate process?) or direct execution by server thread */
                data.Nth, MPI_DOUBLE, CALL_BY_COP,
                1, MPI_DOUBLE, CALL_BY_COP,
                1, MPI_DOUBLE, CALL_BY_COP,
                point, &F,&prior);

                torc_waitall3();    /* wait without releasing the worker */
#endif

    }


    void TmcmcEngine::dump_curgen_db() {

        FILE *fp;
        char fname[256];

        sprintf(fname, "curgen_db_%03d.txt", runinfo.Gen);
        fp = fopen(fname, "w");
        for (int pos = 0; pos < curgen_db.entries; pos++) {

            for (int i = 0; i < data.Nth; ++i) {
                fprintf(fp, "%20.16lf ", curgen_db.entry[pos].point[i]);
            }
            fprintf(fp, "%20.16lf ", curgen_db.entry[pos].F);
            fprintf(fp, "%20.16lf ", curgen_db.entry[pos].prior);
            fprintf(fp,"\n");
        }
        fclose(fp);
    }


    int TmcmcEngine::load_curgen_db() {

        FILE *fp;
        char fname[256];
        sprintf(fname, "curgen_db_%03d.txt", runinfo.Gen);
        fp = fopen(fname, "r");
        if (fp == NULL) {
            printf("DB file: %s not found!!!\n", fname);
            exit(1);
            return 1;
        }

        curgen_db.entries = 0;
        char line[1024];
        while (fgets(line, 1024, fp) != NULL)
            curgen_db.entries++;

        fclose(fp);	/* fseek...*/
        fp = fopen(fname, "r");

        int pos;
        for (pos = 0; pos < curgen_db.entries; ++pos) {
            for (int i = 0; i < data.Nth; ++i) {
                if (curgen_db.entry[pos].point == NULL) 
                    curgen_db.entry[pos].point = new double[data.Nth];
                fscanf(fp, "%lf", &curgen_db.entry[pos].point[i]);
            }
            fscanf(fp, "%lf", &curgen_db.entry[pos].F);
            fscanf(fp, "%lf", &curgen_db.entry[pos].prior);
        }
        fclose(fp);

        return 0;
    }


    void TmcmcEngine::init_curres_db() {
        pthread_mutex_init(&curres_db.m, NULL);
        curres_db.entries = 0;
        curgen_db.entry   = new cgdbp_t[(data.MinChainLength+1)*data.PopSize]; 
    }

    // TODO: is this needed (DW)?
    void TmcmcEngine::update_curres_db(double point[/* EXPERIMENTAL_RESULTS */], double F) {
        
#if (EXPERIMENTAL_RESULTS <=0)
        return; 
#endif
        pthread_mutex_lock(&curres_db.m);
        int pos = curres_db.entries;
        curres_db.entries++;
        pthread_mutex_unlock(&curres_db.m);

        if (curres_db.entry[pos].point == NULL) 
            curres_db.entry[pos].point = new double[EXPERIMENTAL_RESULTS+1];

        for (int i = 0; i < EXPERIMENTAL_RESULTS; ++i) 
            curres_db.entry[pos].point[i] = point[i];
        
        curres_db.entry[pos].F = F;	
    }


    void TmcmcEngine::taskfun(const double *x, int *pN, double *res, int winfo[4]) {
        double f;
        int N = *pN;

        inc_nfc(); // (TODO: include this again, use singleton or so (DW))   /* increment function call counter*/

        f = fitfun::fitfun(x, N, (void *)NULL, winfo);
#if (EXPERIMENTAL_RESULTS > 0)    /* peh: decide about this (results should be passed as argument to fitfun) */
            double results[EXPERIMENTAL_RESULTS];
            for (int i = 0; i < EXPERIMENTAL_RESULTS; ++i) {
                if (i < data.Nth)
                    results[i] = x[i];
                else
                    results[i] = 0.0;
            }
            torc_update_curres_db(results, f);
#endif

        *res = f;
        return;
    }


    void TmcmcEngine::evaluate_F(double point[], double *Fval, int worker_id, 
                                    int gen_id, int chain_id, int step_id, int ntasks) {
        double F;
        int winfo[4];
        int dim = data.Nth;

        winfo[0] = gen_id;
        winfo[1] = chain_id;
        winfo[2] = step_id;
        winfo[3] = 0;

#if VERBOSE
            printf("running on worker %d\n", worker_id);
#endif

        taskfun(point, &dim, &F, winfo);

        *Fval = F;
    }

    void TmcmcEngine::initchaintask(double in_tparam[], double *out_tparam, int winfo[4]) {
        int gen_id   = winfo[0];
        int chain_id = winfo[1];

        long   me = torc_worker_id();
        double point[data.Nth], fpoint;

        for (int i = 0; i < data.Nth; ++i)
            point[i] = in_tparam[i];

        evaluate_F(point, &fpoint, me, gen_id, chain_id, 0, 1);

        double logprior = prior.eval_logpdf(point); 

        /* update current db entry */
        torc_update_curgen_db( point, fpoint, logprior );
        if (data.ifdump) torc_update_full_db(point, fpoint, NULL, 0, 0);
        *out_tparam = fpoint;    /* currently not required, the result is already in the db*/

        return;
    }

    void TmcmcEngine::calculate_statistics(double flc[], int nselections, 
                                           int gen, unsigned int sel[]) {

        int display = data.options.Display;    
        
        int *num              = data.Num;
        double *coefVar       = runinfo.CoefVar;
        double *p             = runinfo.p;
        double *logselections = runinfo.logselections;
        
        /*double pflag = 0;*/
        double fmin = 0, xmin = 0;
        bool conv = 0;

#ifdef _USE_FMINCON_
            conv = fmincon(flc, curgen_db.entries, p[gen], data.TolCOV, &xmin, &fmin,
                    data.options);
            if (display) printf("calculate_statistics: \
                fmincon conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
#endif

#ifdef _USE_FMINSEARCH_
            if (!conv){
                conv = fminsearch(flc, curgen_db.entries, p[gen], data.TolCOV, &xmin, &fmin, 
                        data.options);
                if (display) printf("calculate_statistics: \
                    fminsearch conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
            }
#endif

#ifdef _USE_FZEROFIND_
            if (!conv) {
                conv = fzerofind(flc, curgen_db.entries, p[gen], data.TolCOV, &xmin, &fmin,
                        data.options);
                if (display) printf("calculate_statistics: \
                    fzerofind conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
            }
#endif


        /* gen: next generation number */
        unsigned int j = gen+1;

        if ( conv && (xmin > p[gen] /* + data.MinStep */  ) ) {
            p[j]       = xmin;
            coefVar[j] = fmin;
        } else {
            p[j]       = p[gen] + data.MinStep;
            coefVar[j] = coefVar[gen];
        }

        if (p[j] > 1) {
            /*pflag=p[j-1];*/
            p[j]       = 1;
            coefVar[j] = tmcmc_objlogp(p[j], flc, curgen_db.entries, 
                                        p[j-1], data.TolCOV);
            num[j] = 0; // TODO: data.LastNum;
        }

        /* Compute weights and normalize*/
        unsigned int i;

        double *flcp  = new double[curgen_db.entries]; 	
        for (i = 0; i < curgen_db.entries; ++i)
            flcp[i] = flc[i]*(p[j]-p[j-1]);


        const double fjmax = gsl_stats_max(flcp, 1, curgen_db.entries);
        double *weight     = new double[curgen_db.entries]; 
        for (i = 0; i < curgen_db.entries; ++i)
            weight[i] = exp( flcp[i] - fjmax );

        if (display>2) print_matrix((char *)"weight", weight, curgen_db.entries);

        double sum_weight = std::accumulate(weight, weight+curgen_db.entries, 0.0);

        double *q = new double[curgen_db.entries];
        for (i = 0; i < curgen_db.entries; ++i)
            q[i] = weight[i]/sum_weight;

        if (display>2) print_matrix((char *)"runinfo_q", q, curgen_db.entries);

        logselections[gen] = log(sum_weight) + fjmax - log(curgen_db.entries);

        if (display) print_matrix((char *)"logselections", logselections, gen+1);

        double mean_q = gsl_stats_mean(q, 1, curgen_db.entries);
        double std_q  = gsl_stats_sd_m(q, 1, curgen_db.entries, mean_q);

        coefVar[gen] = std_q/mean_q;

        if (display) print_matrix((char *)"CoefVar", coefVar, gen+1);

        unsigned int N = 1;

        unsigned int *nn = new unsigned int[curgen_db.entries];

        for (i = 0; i < curgen_db.entries; ++i) sel[i] = 0;

        if (nselections == 0) nselections = curgen_db.entries; /* n;*/
        N = nselections;
        multinomialrand (curgen_db.entries, N, q, nn);
        for (i = 0; i < curgen_db.entries; ++i) sel[i]+=nn[i];

        if (display>2) {
            printf("\n s = [");
            for (i = 0; i < curgen_db.entries; ++i) printf("%d ", sel[i]);
            printf("]\n");
        }

        /* compute SS */
        unsigned int PROBDIM = data.Nth;

        double mean_of_theta[PROBDIM];

        for (i = 0; i < PROBDIM; ++i) {
            mean_of_theta[i] = 0;
            for (j = 0; j < curgen_db.entries; ++j) mean_of_theta[i]+=curgen_db.entry[j].point[i]*q[j];

            runinfo.meantheta[gen][i] = mean_of_theta[i];
        }

        if (display) print_matrix((char *)"mean_of_theta", mean_of_theta, PROBDIM);

        double meanv[PROBDIM];
        for (i = 0; i < PROBDIM; ++i) {
            meanv[i] = mean_of_theta[i];
        }

        for (i = 0; i < PROBDIM; ++i) {
            for (j = i; j < PROBDIM; ++j) {
                double s = 0;
                for (unsigned int k = 0; k < curgen_db.entries; ++k) {
                    s += q[k]*(curgen_db.entry[k].point[i]-meanv[i])*(curgen_db.entry[k].point[j]-meanv[j]);
                }
                runinfo.SS[i][j] = runinfo.SS[j][i] = s;
            }
        }

#ifdef CHECK_POSDEF
        int fixed = make_posdef(runinfo.SS[0], PROBDIM, 2);
        if (fixed) {
            printf("WARNING: runinfo.SS was forced to become positive definite\n");
        }
#endif

        if (display) print_matrix_2d((char *)"runinfo.SS", runinfo.SS, PROBDIM, PROBDIM);

        delete [] flcp;
        delete [] weight;
        delete [] q;
        delete [] nn;
    }


    void TmcmcEngine::check_for_exit() {
        int val, exitgen = -1;
        char *s;

        s = (char *) getenv("EXIT_GEN");
        if (s != 0 && sscanf(s, "%d", &val) == 1 && val >= 0)
            exitgen = val;

        if (exitgen == runinfo.Gen) {
            printf("Read Exit Envrironment Variable!!!\n");

#ifdef _USE_TORC_
                torc_finalize();
#endif
            exit(1); // TODO: can we do this smoother? (DW)
        }

        FILE *fp;
        fp = fopen("exit.txt", "r");
        if (fp != NULL) {
            printf("Found Exit File!!!\n");
            //unlink("exit.txt"); TODO: reinsert? (DW)
#ifdef _USE_TORC_
            torc_finalize();
#endif
            exit(1); // TODO: can we do this smoother? (DW)
        }
    }


    void TmcmcEngine::precompute_chain_covariances(const cgdbp_t* leader,double** init_mean, double** chain_cov, int newchains)
    {
        printf("Precomputing covariances for the current generation...\n");

        int D = data.Nth;
        int N = curgen_db.entries;

        double my_time = clock();

        // allocate space
        int* nn_ind        = new int[newchains];
        int* nn_count      = new int[newchains];
        double* diam       = new double[D]; 
        double* chain_mean = new double[D];
        gsl_matrix* work   = gsl_matrix_alloc(D, D);

        // find diameters
        for (int d = 0; d < D; ++d) {
            double d_min = +1e6;
            double d_max = -1e6;
            for (int pos = 0; pos < N; ++pos) {
                double s = curgen_db.entry[pos].point[d];
                if (d_min > s) d_min = s;
                if (d_max < s) d_max = s;
            }
            diam[d] = d_max-d_min;
            printf("Diameter %d: %.6lf\n", d, diam[d]);
        }

        int ind, pos;
        int status = 0;
        double ds = 0.05;
        for (double scale = 0.1; scale <= 1.0; scale += ds) {
            // find neighbors in a rectangle - O(N^2)
            for (pos = 0; pos < newchains; ++pos) {
                nn_count[pos] = 0;
                double* curr = leader[pos].point;
                for (int i = 0; i < N; ++i) {
                    double* s = curgen_db.entry[i].point;
                    if (in_rect(curr, s, diam, scale, D)) {
                        nn_ind[pos*N+nn_count[pos]] = i;
                        nn_count[pos]++;
                    }
                }
            }

            // compute the covariances
            for (pos = 0; pos < newchains; ++pos) {
                for (int d = 0; d < D; ++d) {
                    chain_mean[d] = 0;
                    for (int k = 0; k < nn_count[pos]; ++k) {
                        ind = nn_ind[pos*N+k];
                        chain_mean[d] += curgen_db.entry[ind].point[d];
                    }
                    chain_mean[d] /= nn_count[pos];
                }

                for (int i = 0; i < D; ++i)
                    for (int j = 0; j < D; ++j) {
                        double s = 0;
                        for (int k = 0; k < nn_count[pos]; k++) {
                            ind = nn_ind[pos*N+k];
                            s  += (curgen_db.entry[ind].point[i]-chain_mean[i]) *
                                    (curgen_db.entry[ind].point[j]-chain_mean[j]);
                        }
                        chain_cov[pos][i*D+j] = chain_cov[pos][j*D+i] = s/nn_count[pos];
                    }

                // check if the matrix is positive definite
                for (int i = 0; i < D; ++i)
                    for (int j = 0; j < D; ++j){
                        double s = chain_cov[pos][i*D+j];
                        gsl_matrix_set(work, i, j, s);
                    }
                gsl_set_error_handler_off();
                status = gsl_linalg_cholesky_decomp(work);
                if (status == GSL_SUCCESS) break;
            }
        }

        if (status != GSL_SUCCESS) {
            for (int i = 0; i < D; ++i)
                for (int j = 0; j < D; ++j)
                    chain_cov[pos][i*D+j] = data.bbeta*runinfo.SS[i][j];
        }

        // deallocate space
        delete[] nn_ind;
        delete[] nn_count;
        delete[] diam;
        delete[] chain_mean;
        gsl_matrix_free(work);

        my_time = (clock() - my_time) / CLOCKS_PER_SEC;
        printf("Covariance computation time: %.2lf sec\n", my_time);
    }


    int TmcmcEngine::compute_candidate(double candidate[], double chain_mean[], double var) {
        double bSS[data.Nth*data.Nth];

        for (int i = 0; i < data.Nth; ++i)
            for (int j = 0; j < data.Nth; ++j)
                bSS[i*data.Nth+j]= data.bbeta*runinfo.SS[i][j];


        mvnrnd(chain_mean, (double *)bSS, candidate, data.Nth);

        int idx = 0;
        for (; idx < data.Nth; ++idx) {
            if (isnan(candidate[idx])) {
                printf("!!!!  isnan in candidate point!\n");
                exit(1);
                break;
            }
            if ((candidate[idx] < data.lowerbound[idx]) || 
                    (candidate[idx] > data.upperbound[idx])) break;
        }

        if (idx < data.Nth) return -1;

        return 0;    // all good
    }


    int TmcmcEngine::compute_candidate_cov(double candidate[], double chain_mean[], 
                                            double chain_cov[]) {

        mvnrnd(chain_mean, (double *)chain_cov, candidate, data.Nth);
        for (int i = 0; i < data.Nth; ++i) {
            if (isnan(candidate[i])) {
                printf("!!!!  isnan in candidate point!\n");
                exit(1);
                break;
            }
            if ((candidate[i] < data.lowerbound[i])||(candidate[i] > data.upperbound[i])) return -1;
        }
        return 0;
    }


    void TmcmcEngine::chaintask(double in_tparam[], int *pnsteps, double *out_tparam, int winfo[4],
            double *init_mean, double *chain_cov) {
        
        int nsteps   = *pnsteps;
        int gen_id   = winfo[0];
        int chain_id = winfo[1];

        long me = torc_worker_id();

        double leader[data.Nth], loglik_leader, logprior_leader;            /* old*/
        double candidate[data.Nth], loglik_candidate, logprior_candidate;   /* new*/

        // get initial leader and its value
        for (int i = 0; i < data.Nth; ++i) leader[i] = in_tparam[i];    
        loglik_leader   = *out_tparam;    
        logprior_leader = prior.eval_logpdf(leader); 

        double pj = runinfo.p[runinfo.Gen];

        int burn_in = data.burn_in;

        for (int step = 0; step < nsteps + burn_in; ++step) {
            double chain_mean[data.Nth];
            if (step == 0)
                for (int i = 0; i < data.Nth; ++i) chain_mean[i] = init_mean[i];
            else
                for (int i = 0; i < data.Nth; ++i) chain_mean[i] = leader[i];

            //printf("---> %d -  %ld/%d   ---  %p  \n", chain_id, me, torc_i_num_workers(), priors ); fflush(NULL);

#if 0
                int fail = compute_candidate_cov(candidate, chain_mean, chain_cov);
#else
                int fail = compute_candidate(candidate, chain_mean, 1); // I keep this for the moment, for performance reasons
#endif

            if (!fail){

                evaluate_F(candidate, &loglik_candidate, me, gen_id, chain_id, step, 1);    // this can spawn many tasks

                if (data.ifdump && step >= burn_in) torc_update_full_db(candidate, loglik_candidate, NULL, 0, 0);
                                                    // last argument should be 1 if it is a surrogate

                // decide
                logprior_candidate = prior.eval_logpdf(candidate); 
                double L = exp((logprior_candidate-logprior_leader)+(loglik_candidate-loglik_leader)*pj);

                if (L > 1) L = 1;
                double P = uniformrand(0,1);
                if (P < L) {
                    // new leader
                    for (int i = 0; i < data.Nth; ++i) leader[i] = candidate[i]; 
                    loglik_leader = loglik_candidate;
                    if (step >= burn_in) {
                        logprior_leader = prior.eval_logpdf(leader); 
                        torc_update_curgen_db(leader, loglik_leader, logprior_leader);
                    }
                }
                else {
                    // increase counter or add the leader again in curgen_db
                    if (step >= burn_in) {
                        logprior_leader = prior.eval_logpdf(leader); 
                        torc_update_curgen_db(leader, loglik_leader, logprior_leader);
                    }

                }
            }
            else{
                // increase counter or add the leader again in curgen_db
                if (step >= burn_in){
                    logprior_leader = prior.eval_logpdf(leader);
                    torc_update_curgen_db(leader, loglik_leader, logprior_leader);
                }
            }
        }

        return;
    }

    int TmcmcEngine::prepare_newgen(int nchains, cgdbp_t *leaders) {
        /* process curgen_db -> calculate statitics */
        /* compute probs based on F values */
        /* draw new samples (nchains or user-specified) */
        /* find unique samples: fill the (new) leaders table */
        /* count how many times they appear -> nsteps */
        /* return the new sample size (number of chains) */

        int i, p;

        int n = curgen_db.entries;

        double *fj 		  = new double[n];
        unsigned int *sel = new unsigned int[n];

        double **g_x = new double*[data.Nth];
        for (i = 0; i < data.Nth; ++i) g_x[i] = new double[n]; 

        if(0)
        {
            double **x = g_x;

            for (p = 0; p < data.Nth; p++) {
                for (i = 0; i < n; ++i) {
                    x[p][i] = curgen_db.entry[i].point[p];
                }
            }

            double meanx[data.Nth], stdx[data.Nth];
            for (p = 0; p < data.Nth; p++) {
                meanx[p] = gsl_stats_mean(x[p], 1, n);
                stdx[p]  = gsl_stats_sd_m(x[p], 1, n, meanx[p]);
            }

            if(data.options.Display){
                printf("prepare_newgen: CURGEN DB (COMPLE) %d\n", runinfo.Gen);
                print_matrix("means", meanx, data.Nth);
                print_matrix("std", stdx, data.Nth);
            }

        }



        if (1)
        {
            double **uniques = g_x;
            int un = 0, unflag, j;

            for( p = 0; p < data.Nth; ++p )
                uniques[p][un] = curgen_db.entry[0].point[p];

            un++;
            for (i = 1; i < n; ++i){
                double xi[data.Nth];
                for (p = 0; p < data.Nth; ++p)
                    xi[p] = curgen_db.entry[i].point[p];

                unflag = 1;                 /* is this point unique? */
                for (j = 0; j < un; ++j){   /* compare with  previous uniques */
                    for (p = 0; p < data.Nth; ++p){
                        if (fabs(xi[p]-uniques[p][j]) > 1e-6) {
                            break;          /* they differ, check next  */
                        }
                        unflag = 0;         /* not unique */
                    }

                    if (unflag == 0) break; /* not unique - stop comparison */
                }

                if (unflag){                /* unique, put it in the table */
                    for (p = 0; p < data.Nth; ++p) uniques[p][un] = xi[p];
                    un++;
                }
            }

            runinfo.currentuniques[runinfo.Gen] = un;
            runinfo.acceptance[runinfo.Gen]     = (1.0*runinfo.currentuniques[runinfo.Gen])/data.Num[runinfo.Gen]; /* check this*/

            double meanu[data.Nth], stdu[data.Nth];
            for (p = 0; p < data.Nth; ++p) {
                meanu[p] = gsl_stats_mean(uniques[p], 1, n);
                stdu[p]  = gsl_stats_sd_m(uniques[p], 1, n, meanu[p]);
            }

            printf("prepare_newgen: CURGEN DB (UNIQUE) %d: [un = %d]\n", runinfo.Gen, un);
            if(data.options.Display){
                print_matrix((char *)"uniques mean", meanu, data.Nth);
                print_matrix((char *)"uniques std", stdu, data.Nth);
            }

        } /* end block*/

        {
            double t0 = torc_gettime();
            for (i = 0; i < n; ++i)
                fj[i] = curgen_db.entry[i].F;    /* separate point from F ?*/
            double t1 = torc_gettime();
            calculate_statistics(fj, data.Num[runinfo.Gen], runinfo.Gen, sel);
            double t2 = torc_gettime();
            printf("prepare_newgen: init + calc stats : %lf + %lf = %lf seconds\n", t2-t1, t1-t0, t2-t0);
        }

        int newchains = 0;
        sort_t *list = new sort_t[n];
        for (i = 0; i < n; ++i) {
            list[i].idx  = i;
            list[i].nsel = sel[i];
            list[i].F    = curgen_db.entry[i].F;
            if (sel[i] != 0) newchains++;
        }

#if VERBOSE
        printf("Points before qsort\n");
        for (i = 0; i < n; ++i) 
            printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
        
#endif

        qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
        printf("Points after qsort\n");
        for (i = 0; i < n; ++i) 
            printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
        
#endif

#if 0 // TODO: check what is going on here (DW)   
        /* UPPER THRESHOLD */
        /* peh:check this */
        /* breaking long chains */
        int initial_newchains = newchains;
        int h_threshold = data.MaxChainLength;    /* peh: configuration file + more balanced lengths */
        for (i = 0; i < initial_newchains; ++i) {
            if (list[i].nsel > h_threshold) {
                while (list[i].nsel > h_threshold) {
                    list[newchains] = list[i];
                    list[newchains].nsel = h_threshold;
                    list[i].nsel = list[i].nsel - h_threshold;
                    newchains++;
                }
            }
        }

        qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
        printf("Points broken\n");
        for (i = 0; i < n; ++i) 
            printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
        
#endif

#endif

#if 0 // TODO: check what is going on here (DW)   
        /* LOWER THRESHOLD */
        /* single to double step chains */
        int l_threshold = data.MinChainLength;    /* peh: configuration file + more balanced lengths */
        for (i = 0; i < newchains; ++i) {
            if ((list[i].nsel > 0)&&(list[i].nsel < l_threshold)) {
                list[i].nsel = l_threshold;
            }
        }

        qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
        printf("Points advanced\n");
        for (i = 0; i < n; ++i) {
            printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
        }
#endif

#endif

        int ldi = 0;                    /* leader index */
        for (i = 0; i < n; ++i) {       /* newleader */
            if (list[i].nsel != 0) {
                int idx = list[i].idx;
                for (p = 0; p < data.Nth ; p++) {
                    leaders[ldi].point[p] = curgen_db.entry[idx].point[p];
                }
                leaders[ldi].F = curgen_db.entry[idx].F;
                leaders[ldi].nsel = list[i].nsel;
                ldi++;
            }
        }

        delete[] list;

        for (i = 0; i < newchains; ++i) leaders[i].queue = -1;    /* rr*/


#ifdef _USE_TORC_

#if VERBOSE
        printf("Leaders before partitioning\n");
        for (i = 0; i < newchains; ++i) {
            printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
        }
#endif

        /* cool and greedy partitioning ala Panos-- ;-) */

        int nworkers  = torc_num_workers();
        int *workload = new int[nworkers];

        for (i = 0; i < newchains; ++i) {
            int least_loader_worker = compute_min_idx_i(workload, nworkers);
            leaders[i].queue = least_loader_worker;
            workload[least_loader_worker] += leaders[i].nsel;
        }

        print_matrix_i("initial workload", workload, nworkers);
        delete[] workload;

#if VERBOSE
        printf("Leaders after partitioning\n");
        for (i = 0; i < newchains; ++i) {
            printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
        }
#endif

#endif
        if (1)
        {
            double **x = g_x;
            for (i = 0; i < newchains; ++i) {
                for (p = 0; p < data.Nth; p++) {
                    x[p][i] = leaders[i].point[p];
                }
            }

            double meanx[data.Nth], stdx[data.Nth];
            for (p = 0; p < data.Nth; p++) {
                meanx[p] = gsl_stats_mean(x[p], 1, newchains);
                stdx[p]  = gsl_stats_sd_m(x[p], 1, newchains, meanx[p]);
            }

            printf("prepare_newgen: CURGEN DB (LEADER) %d: [nlead=%d]\n", runinfo.Gen, newchains);
            if(data.options.Display){
                print_matrix("means", meanx, data.Nth);
                print_matrix("std", stdx, data.Nth);
            }
        }

        if (data.use_local_cov)
            precompute_chain_covariances(leaders, data.init_mean, data.local_cov, newchains);

        curgen_db.entries = 0;
        printf("prepare_newgen: newchains=%d\n", newchains);

        for (i = 0; i < data.Nth; ++i) delete g_x[i]; 
        delete[] g_x;
        delete[] fj;
        delete[] sel;

        return newchains;
    }

    void TmcmcEngine::print_runinfo() {
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
        printf("runinfo.Gen = \n\n   %d\n\n", runinfo.Gen);
        print_matrix("runinfo.p", runinfo.p, runinfo.Gen+1); 
        print_matrix_2d("runinfo.SS", runinfo.SS, data.Nth, data.Nth);
        printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    }

    void TmcmcEngine::bcast_runinfo() {
#ifdef _USE_TORC_
		MPI_Bcast(runinfo.SS[0], data.Nth*data.Nth, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(runinfo.p, data.MaxStages, MPI_DOUBLE, 0, MPI_COMM_WORLD);    /* just p[Gen]*/
		MPI_Bcast(&runinfo.Gen, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    }

    void TmcmcEngine::spmd_update_runinfo()    /* step*/
    {
#ifdef _USE_TORC_
        if (torc_num_nodes() == 1) return;
        for (int i = 0; i < torc_num_nodes(); i++) {
            torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())bcast_runinfo, 0);
        }
        torc_waitall();
#endif
    }

} //namespace tmcmc
