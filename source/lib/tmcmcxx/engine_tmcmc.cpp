#include <math.h>
#include <cstring>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <vector>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>

#include "fitfun.hpp"
#include "tmcmc_obj_fmin.hpp"
#include "engine_tmcmc.hpp"

using namespace priors;

namespace tmcmc
{

TmcmcEngine::TmcmcEngine(fitfun::IFitfun * ifitfun_ptr, Method method, const char * ftmcmcpar, const char * fpriorpar ) : 
    _method(method),
    data(data_t(ftmcmcpar)),
    nchains(data.Num[0]),
    leaders(new cgdbp_t[data.PopSize]),
    prior(fpriorpar),
    ifitfun_ptr_(ifitfun_ptr)
{

    gsl_error_handler = gsl_set_error_handler_off ();
    
        for (int i = 0; i< data.PopSize; ++i)
        leaders[i].point = new double[data.Nth];

    if (_method == Manifold) {
        for (int i = 0; i< data.PopSize; ++i) {
            leaders[i].gradient = new double[data.Nth];
            leaders[i].cov      = new double[data.Nth*data.Nth];
            leaders[i].evec     = new double[data.Nth*data.Nth];
            leaders[i].eval     = new double[data.Nth];
        }
    }

    curres_db.entries = 0;

    ifitfun_ptr_->initialize(0, NULL);
    init();
}

TmcmcEngine::~TmcmcEngine()
{
    for (int i = 0; i< data.PopSize; ++i) delete[] leaders[i].point;

    if (_method == Manifold) {
        for (int i = 0; i< data.PopSize; ++i) {
            delete[] leaders[i].gradient;
            delete[] leaders[i].cov;
            delete[] leaders[i].evec;
            delete[] leaders[i].eval;
        }
    }

    delete[] leaders;

    for(int i = 0; i< data.MaxStages; ++i) delete [] runinfo.meantheta[i];

    ifitfun_ptr_->finalize();

    //TODO: what else needs to be freed/deleted?? (DW)
    
    return;
    
}

void TmcmcEngine::run()
{

    if (runinfo.p[runinfo.Gen] == 1.0) {
        printf("p == 1 from previous run, nothing more to do\n");
        return;
    }

    prepare_newgen(nchains, leaders);
    spmd_update_runinfo();
    if (data.options.Display > 0) print_runinfo();

    while(runinfo.p[runinfo.Gen] < 1.0 && ++runinfo.Gen < data.MaxStages) {
        evalGen();
        if(data.options.Display > 0) print_runinfo();
    }

    print_matrix((char *)"runinfo.p", runinfo.p, runinfo.Gen+1);
    print_matrix((char *)"runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
    print_matrix_i((char *)"runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
    print_matrix((char *)"runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
    print_matrix((char *)"runinfo.logselection", runinfo.logselections, runinfo.Gen);

    print_runinfo();
    runinfo_t::save(runinfo, data.Nth, data.MaxStages);

    double logEvidence = compute_sum(runinfo.logselections, runinfo.Gen+1);
    printf("logEvidence = %f\n", logEvidence);

    FILE *fp;
    fp = fopen("log_evidence.txt", "w");
    fprintf(fp, "%lf\n", logEvidence);
    fclose(fp);

    if (data.icdump) {
        printf("lastgen = %d\n", runinfo.Gen);
        char cmd[256];
        sprintf(cmd, "cp curgen_db_%03d.txt final.txt", runinfo.Gen);
        system(cmd);
    }

    printf("total function calls = %d\n", get_tfc());

#if defined(_USE_TORC_)
    torc_finalize();
#endif

    return;
}


void TmcmcEngine::init()
{

    gt0 = t0 = torc_gettime();

    init_curgen_db();
    init_curres_db();
    init_full_db();

    runinfo_t::init(runinfo, data.Nth, data.MaxStages);

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

void TmcmcEngine::evalGen()
{

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


            for (int i = 0; i < nchains; ++i)
            {
                winfo[0] = runinfo.Gen;
                winfo[1] = i;
                winfo[2] = -1;    /* not used */
                winfo[3] = -1;    /* not used */

                for(int d = 0; d < data.Nth; ++d)
                    in_tparam[d] = leaders[i].point[d];

                nsteps = leaders[i].nsel;

                if (data.use_local_cov) {
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
                                data.beta2*runinfo.SS[d][e];

                    for (int d = 0; d < data.Nth; ++d) init_mean[d] = in_tparam[d];
                }

#ifdef _USE_TORC_
                if (_method == Standard )
                    torc_create(leaders[i].queue, (void (*)())chaintask, 7,
                                data.Nth, MPI_DOUBLE, CALL_BY_COP,
                                1, MPI_INT, CALL_BY_COP,
                                1, MPI_INT, CALL_BY_COP,
                                1, MPI_DOUBLE, CALL_BY_REF,
                                4, MPI_INT, CALL_BY_COP,
                                data.Nth, MPI_DOUBLE, CALL_BY_COP,
                                data.Nth*data.Nth, MPI_DOUBLE, CALL_BY_COP,
                                in_tparam, &data.Nth, &nsteps, leaders[i].F, winfo,
                                init_mean, chain_cov);
                else /* Manifold */ {
                    printf("Manifold method not yet implemented for TORC!!! aborting.. \n");
                    abort();
                }
#else
#ifdef _USE_OPENMP_
                #pragma omp task shared(data, leaders) firstprivate(i, nsteps, in_tparam, winfo, init_mean, chain_cov)
#endif
                if (_method == Standard ) 
                    chaintask( in_tparam, nsteps, leaders[i].F , winfo, init_mean, chain_cov );
                else /* Manifold */ 
                    manifold_chaintask(in_tparam, nsteps, leaders[i].F , 
                            leaders[i].error_flg, leaders[i].posdef, 
                            leaders[i].gradient, leaders[i].cov, 
                            leaders[i].evec, leaders[i].eval, winfo ); 
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
    printf("evalGen: Generation %d: total elapsed time = %lf secs, ",
           runinfo.Gen, gt1-t0);
    printf("generation elapsed time = %lf secs for function calls = %d\n",
           gt1-gt0, g_nfeval);

    reset_nfc();

    if (data.icdump) dump_curgen_db();
    if (data.ifdump) dump_full_db();

    runinfo_t::save(runinfo, data.Nth, data.MaxStages);

    if (data.restart) check_for_exit();

    curres_db.entries = 0;
    prepare_newgen(nchains, leaders);
    spmd_update_runinfo();

    return;
}

bool TmcmcEngine::load_data()
{
    bool res = runinfo_t::load(runinfo, data.Nth, data.MaxStages);
    if (res == 0) {
        load_curgen_db();
        printf("nchains = %d\n", data.Num[0]);
    }
    return res;
}

void TmcmcEngine::sample_from_prior()
{

#ifdef _USE_OPENMP_
    #pragma omp parallel
    {
        printf("Hello from thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
        #pragma omp for
#endif
        int winfo[4];
        for (int i = 0; i < nchains; ++i)
        {
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
                        4, MPI_INT, CALL_BY_COP,
                        in_tparam, &data.Nth, winfo);
#else
#ifdef _USE_OPENMP_
            //#pragma omp task shared(data, leaders) firstprivate(i, in_tparam, winfo)
            #pragma omp task firstprivate(i, winfo, in_tparam) shared(data, leaders)
#endif
            {
                initchaintask(in_tparam, winfo);
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

    printf("sample_from_prior: Generation %d: total elapsed time = %lf sec,\n",
           runinfo.Gen, gt1-t0);
    printf("generation elapsed time = %lf secs for function calls = %d\n",
           gt1-gt0, g_nfeval);

    reset_nfc();

    return;
}


void TmcmcEngine::init_full_db()
{
    pthread_mutex_init(&full_db.m, nullptr);
    full_db.entries = 0;
    full_db.entry   = new dbp_t[data.MaxStages*data.PopSize];
}


void TmcmcEngine::update_full_db(double point[], double F, double *G, int n, int surrogate)
{
    pthread_mutex_lock(&full_db.m);
    int pos = full_db.entries;
    full_db.entries++;
    pthread_mutex_unlock(&full_db.m);

    if (full_db.entry[pos].point == nullptr)
        full_db.entry[pos].point = new double[data.Nth] ;

    for (int i = 0; i < data.Nth; ++i) full_db.entry[pos].point[i] = point[i];
    full_db.entry[pos].F = F;
    full_db.entry[pos].nG = n;

    for (int i = 0; i < n; ++i) full_db.entry[pos].G[i] = G[i];
    full_db.entry[pos].surrogate = surrogate;
}


void TmcmcEngine::torc_update_full_db(double point[], double F, double *G, int n, int surrogate)
{
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


void TmcmcEngine::update_manifold_curgen_db(double point[], double F, double prior, int error_flg, bool posdef, double gradient[], double cov[], double evec[], double eval[])
{
	int i, pos;

	pthread_mutex_lock(&curgen_db.m);
	pos = curgen_db.entries;
	curgen_db.entries++;
	pthread_mutex_unlock(&curgen_db.m);

	if (curgen_db.entry[pos].point == nullptr) curgen_db.entry[pos].point = new double[data.Nth];

	for (i = 0; i < data.Nth; ++i) curgen_db.entry[pos].point[i] = point[i];
	curgen_db.entry[pos].F = F;
	curgen_db.entry[pos].prior = prior;

	curgen_db.entry[pos].error_flg = error_flg;
	curgen_db.entry[pos].posdef = posdef;

	if (curgen_db.entry[pos].gradient == nullptr) curgen_db.entry[pos].gradient = new double[data.Nth];
	for (i = 0; i < data.Nth; ++i) curgen_db.entry[pos].gradient[i] = gradient[i];

	if (curgen_db.entry[pos].cov == nullptr) curgen_db.entry[pos].cov = new double[data.Nth*data.Nth];
	for (i = 0; i < data.Nth*data.Nth; ++i) curgen_db.entry[pos].cov[i] = cov[i];

	if (curgen_db.entry[pos].evec == nullptr) curgen_db.entry[pos].evec = new double[data.Nth*data.Nth];
	for (i = 0; i < data.Nth*data.Nth; ++i) curgen_db.entry[pos].evec[i] = evec[i];

	if (curgen_db.entry[pos].eval == nullptr) curgen_db.entry[pos].eval = new double[data.Nth];
	for (i = 0; i < data.Nth; ++i) curgen_db.entry[pos].eval[i] = eval[i];

    return;
}


void TmcmcEngine::torc_update_manifold_curgen_db(double point[], double F, double prior, int error_flg, bool posdef, double gradient[], double cov[], double evec[], double eval[])
{
	if (torc_node_id() == 0) {
		update_manifold_curgen_db(point, F, prior, error_flg, posdef, gradient, cov, evec, eval);
		return;
	}

#ifdef _USE_TORC_
	torc_create_direct(0, torc_update_curgen_db_mala_task, 9,		/* message to the database manager (separate process?) or direct execution by server thread */
		data.Nth, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_INT, CALL_BY_COP,
		1, MPI_INT, CALL_BY_COP,
		data.Nth, MPI_DOUBLE, CALL_BY_COP,		// VAL...
		data.Nth*data.Nth, MPI_DOUBLE, CALL_BY_COP,
		data.Nth*data.Nth, MPI_DOUBLE, CALL_BY_COP,
		data.Nth, MPI_DOUBLE, CALL_BY_COP,
		point, &F, &prior, &error_flg, &posdef, gradient, cov, evec, eval);
	torc_waitall3();	/* wait without releasing the worker */
#endif
}




void TmcmcEngine::print_full_db()
{
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


void TmcmcEngine::dump_full_db()
{
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


void TmcmcEngine::init_curgen_db()
{
    pthread_mutex_init(&curgen_db.m, nullptr);
    curgen_db.entries = 0;
    curgen_db.entry   = new cgdbp_t[(data.MinChainLength+1)*data.PopSize];
}


void TmcmcEngine::update_curgen_db(double point[], double F, double prior)
{
    pthread_mutex_lock(&curgen_db.m);
    int pos = curgen_db.entries;
    curgen_db.entries++;
    pthread_mutex_unlock(&curgen_db.m);

    if (curgen_db.entry[pos].point == nullptr) curgen_db.entry[pos].point = new double[data.Nth];
    if (curgen_db.entry[pos].point == nullptr) curgen_db.entry[pos].point = new double[data.Nth];

    for (int i = 0; i < data.Nth; ++i) curgen_db.entry[pos].point[i] = point[i];
    curgen_db.entry[pos].F = F;
    curgen_db.entry[pos].prior = prior;
}


void TmcmcEngine::torc_update_curgen_db_task(double point[], double *pF, double *pprior)
{
    double F = *pF;
    double prior = *pprior;
    update_curgen_db(point, F, prior);
}


void TmcmcEngine::torc_update_curgen_db(double point[], double F, double prior)
{
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


void TmcmcEngine::dump_curgen_db()
{

    FILE *fp;
    char fname[256];

    sprintf(fname, "curgen_db_%03d.txt", runinfo.Gen);
    fp = fopen(fname, "w");
    for (size_t pos = 0; pos < curgen_db.entries; pos++) {

        for (int i = 0; i < data.Nth; ++i) {
            fprintf(fp, "%20.16lf ", curgen_db.entry[pos].point[i]);
        }
        fprintf(fp, "%20.16lf ", curgen_db.entry[pos].F);
        fprintf(fp, "%20.16lf ", curgen_db.entry[pos].prior);
        fprintf(fp,"\n");
    }
    fclose(fp);
}


int TmcmcEngine::load_curgen_db()
{

    FILE *fp;
    char fname[256];
    sprintf(fname, "curgen_db_%03d.txt", runinfo.Gen);
    fp = fopen(fname, "r");
    if (fp == nullptr) {
        printf("DB file: %s not found!!!\n", fname);
        exit(1);
        return 1;
    }

    curgen_db.entries = 0;
    char line[1024];
    while (fgets(line, 1024, fp) != nullptr)
        curgen_db.entries++;

    fclose(fp);
    fp = fopen(fname, "r");

    for (size_t pos = 0; pos < curgen_db.entries; ++pos) {
        for (int i = 0; i < data.Nth; ++i) {
            if (curgen_db.entry[pos].point == nullptr)
                curgen_db.entry[pos].point = new double[data.Nth];
            fscanf(fp, "%lf", &curgen_db.entry[pos].point[i]);
        }
        fscanf(fp, "%lf", &curgen_db.entry[pos].F);
        fscanf(fp, "%lf", &curgen_db.entry[pos].prior);
    }
    fclose(fp);

    return 0;
}


void TmcmcEngine::init_curres_db()
{
    pthread_mutex_init(&curres_db.m, nullptr);
    curres_db.entries = 0;
    curgen_db.entry   = new cgdbp_t[(data.MinChainLength+1)*data.PopSize];
}

void TmcmcEngine::evaluate_candidate(double point[], double *Fval, int worker_id,
                             int gen_id, int chain_id, int step_id, int ntasks)
{
    int winfo[4] = { gen_id, chain_id, step_id, 0 };

#if VERBOSE
    printf("running on worker %d\n", worker_id);
#endif

    inc_nfc();

    // TODO: this belongs somewher else for mTMCMC (DW)
    // also not very efficient
    fitfun::return_type *out = new fitfun::return_type;

    *Fval = ifitfun_ptr_->evaluate(point, data.Nth, out, winfo);
    
    delete out;

    return;
}


void TmcmcEngine::manifold_evaluate_candidate(double point[], double *Fval, int *perr, 
                                      bool *pposdef, double grad[], double cov[], 
                                      double evec[], double eval[], int worker_id, 
                                      int gen_id, int chain_id, int step_id, int ntasks)
{
    int Nth = data.Nth;
    int winfo[4] = { gen_id, chain_id, step_id, 0 };

#if VERBOSE
    printf("running on worker %d\n", worker_id);
#endif

    inc_nfc();

    // TODO: this belongs somewher else for mTMCMC (DW)
    fitfun::return_type *result = new fitfun::return_type;
    
    *Fval = ifitfun_ptr_->evaluateM(point, Nth, result, winfo);
    
    *perr    = result->error_flg;
    *pposdef = result->posdef;

    if (result->grad != nullptr) std::copy(result->grad->data, result->grad->data+Nth, grad);
    else std::memset(grad, 0, Nth*sizeof(double));

	if (result->eval != nullptr) std::copy(result->eval->data, result->eval->data+Nth, eval);
    else std::memset(eval, 0, Nth*sizeof(double));

	if (result->evec != nullptr) std::copy(result->evec->data, result->evec->data+Nth*Nth, evec);
	else std::memset(evec, 0, Nth*Nth*sizeof(double));
	
    if (result->cov != nullptr) std::copy(result->cov->data, result->cov->data+Nth*Nth, cov);
	else std::memset(cov, 0, Nth*Nth*sizeof(double));

    delete result;

    return;
}


void TmcmcEngine::initchaintask(double in_tparam[],  int winfo[4])
{
    int gen_id   = winfo[0];
    int chain_id = winfo[1];

    long   me = torc_worker_id();
    double leader[data.Nth], loglik_leader;

    for (int i = 0; i < data.Nth; ++i)
        leader[i] = in_tparam[i];

    double logprior = prior.eval_logpdf(leader);
   
    if( _method == Standard)
    {
        evaluate_candidate( leader, &loglik_leader, me, gen_id, chain_id, 0, 1);
        torc_update_curgen_db( leader, loglik_leader, logprior );

    } else /* Manifold */ {

        int c_err;
        bool c_posdef;
        double c_grad[data.Nth];
        double c_cov[data.Nth*data.Nth];
        double c_evec[data.Nth*data.Nth];
        double c_eval[data.Nth];

        manifold_evaluate_candidate(leader, &loglik_leader, &c_err, &c_posdef, 
                                    c_grad, c_cov, c_evec, c_eval, 
                                    me, gen_id, chain_id, 0, 1);
        torc_update_manifold_curgen_db(leader, loglik_leader, logprior, c_err, 
                                       c_posdef, c_grad, c_cov, c_evec, c_eval);

    }
    
    if (data.ifdump) torc_update_full_db(leader, loglik_leader, nullptr, 0, 0);

    return;
}

void TmcmcEngine::calculate_statistics(double flc[], int nselections,
                                       int gen, unsigned int sel[])
{

    int display = data.options.Display;

    double *coefVar       = runinfo.CoefVar;
    double *p             = runinfo.p;
    double *logselections = runinfo.logselections;

    double fmin = 0, xmin = 0;
    bool conv = 0;

#ifdef _USE_FMINCON_
    conv = fmincon(flc, curgen_db.entries, p[gen], data.TolCOV, &xmin, &fmin,
                   data.options);
    if (display) printf("calculate_statistics: \
                fmincon conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
#endif

#ifdef _USE_FMINSEARCH_
    if (!conv) {
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


    if ( conv && (xmin > p[gen]) ) {
        p[gen+1]       = xmin;
        coefVar[gen+1] = pow(fmin, 0.5) + data.TolCOV;
    } else {
        p[gen+1]       = p[gen] + data.MinStep;
        coefVar[gen+1] = pow(tmcmc_objlogp(p[gen+1], flc, curgen_db.entries,
                                       p[gen], data.TolCOV), 0.5) + data.TolCOV;
    }

    if (p[gen+1] > 1) {
        p[gen+1]       = 1;
        coefVar[gen+1] = pow(tmcmc_objlogp(p[gen+1], flc, curgen_db.entries,
                                       p[gen], data.TolCOV), 0.5) + data.TolCOV;
    }

    /* Compute weights and normalize*/
    unsigned int i;

    double *flcp = new double[curgen_db.entries];
    for (i = 0; i < curgen_db.entries; ++i) flcp[i] = flc[i]*(p[gen+1]-p[gen]);


    const double fjmax = gsl_stats_max(flcp, 1, curgen_db.entries);
    double *weight     = new double[curgen_db.entries];
    for (i = 0; i < curgen_db.entries; ++i)
        weight[i] = exp( flcp[i] - fjmax );
    if (display>2) print_matrix("weight", weight, curgen_db.entries);

    delete [] flcp;

    //TODO: logselections[gen+1] or ..[gen]? (DW)
    //logselection is log(Sj), what is logselections[0] if we switch to former?
    double sum_weight = std::accumulate(weight, weight+curgen_db.entries, 0.0);
    if (display) printf("sum_weight %f\n", sum_weight);
    
    logselections[gen] = log(sum_weight) + fjmax - log(curgen_db.entries);
    if (display) print_matrix("logselections", logselections, gen);

    double *q = new double[curgen_db.entries];
    for (i = 0; i < curgen_db.entries; ++i) q[i] = weight[i]/sum_weight;
    if (display>2) print_matrix("runinfo_q", q, curgen_db.entries);

    delete [] weight;

    if (display) print_matrix("CoefVar", coefVar, gen+1);

    for (i = 0; i < curgen_db.entries; ++i) sel[i] = 0;

    if (nselections == 0) nselections = curgen_db.entries;
    
    multinomialrand (curgen_db.entries, nselections, q, sel);

    if(display>2) {
        printf("\n s = [");
        int nonzeros = 0;
        for (i = 0; i < curgen_db.entries; ++i) { 
            printf("%d ", sel[i]); 
            if (sel[i] != 0) nonzeros++;
        }
        printf("] (total nonzeros %d )\n", nonzeros );
    }

    /* compute SS */
    unsigned int j;
    unsigned int PROBDIM = data.Nth;

    for (i = 0; i < PROBDIM; ++i) {
        runinfo.meantheta[gen][i] = 0;
        for (j = 0; j < curgen_db.entries; ++j) runinfo.meantheta[gen][i] += curgen_db.entry[j].point[i]*q[j];
    }

    if (display) print_matrix("runinfo.meantheta", runinfo.meantheta[gen], PROBDIM);

    for (i = 0; i < PROBDIM; ++i) {
        for (j = i; j < PROBDIM; ++j) {
            double s = 0;
            for (unsigned int k = 0; k < curgen_db.entries; ++k) {
                s += q[k]*(curgen_db.entry[k].point[i]-runinfo.meantheta[gen][i])*(curgen_db.entry[k].point[j]-runinfo.meantheta[gen][j]);
            }
            runinfo.SS[i][j] = runinfo.SS[j][i] = s;
        }
    }

    delete [] q;

#ifdef CHECK_POSDEF
    int fixed = make_posdef(runinfo.SS[0], PROBDIM, 2);
    if (fixed) printf("WARNING: runinfo.SS was forced to become positive definite\n");
#endif

    if (display) print_matrix_2d("runinfo.SS", runinfo.SS, PROBDIM, PROBDIM);

}


void TmcmcEngine::check_for_exit()
{
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
        exit(1);
    }

    FILE *fp;
    fp = fopen("exit.txt", "r");
    if (fp != nullptr) {
        printf("Found Exit File!!!\n");
        fclose(fp);
#ifdef _USE_TORC_
        torc_finalize();
#endif
        exit(1);
    }
}

//TODO: this func has not been checked (DW)
void TmcmcEngine::precompute_chain_covariances(const cgdbp_t* leader,double** init_mean, double** chain_cov, int newchains)
{
    bool display = data.options.Display;

    printf("Precomputing chain covariances for the current generation...\n");

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
        if (display) printf("Diameter %d: %.6lf\n", d, diam[d]);
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
                for (int j = 0; j < D; ++j) {
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
                chain_cov[pos][i*D+j] = data.beta2*runinfo.SS[i][j];
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


bool TmcmcEngine::compute_candidate(double candidate[], double chain_mean[])
{
    double bSS[data.Nth*data.Nth];

    for (int i = 0; i < data.Nth; ++i)
        for (int j = 0; j < data.Nth; ++j)
            bSS[i*data.Nth+j]= data.beta2*runinfo.SS[i][j];

    mvnrnd(chain_mean, (double *)bSS, candidate, data.Nth);

    int idx = 0;
    for (; idx < data.Nth; ++idx) {
        if (isnan(candidate[idx])) {
            printf("!!!!  compute_candidate: isnan in candidate point!\n");
            break;
        }
        if ((candidate[idx] < data.lowerbound[idx]) ||
                (candidate[idx] > data.upperbound[idx])) break;
    }

    if (idx < data.Nth) {
        runinfo.outside[runinfo.Gen] ++;
        return false;
    }
    
    return true;
}


bool TmcmcEngine::compute_candidate_cov(double candidate[], double chain_mean[],
                                       double chain_cov[])
{
    mvnrnd(chain_mean, (double *)chain_cov, candidate, data.Nth);
    for (int i = 0; i < data.Nth; ++i) {
        if (isnan(candidate[i])) {
            printf("!!!!  compute_candidate_cov: isnan in candidate point!\n");
            break;
        }
        if ((candidate[i] < data.lowerbound[i])||(candidate[i] > data.upperbound[i])) return false;
    }
    return true;
}


bool TmcmcEngine::compute_manifold_candidate(double candidate[], double leader[], double eps, double *SIG, double *grad, bool posdef)
{
	int i;
	double epsSIG[data.Nth*data.Nth];

	for (i = 0; i < data.Nth*data.Nth; ++i) epsSIG[i]= eps*SIG[i];

	double tmp[data.Nth];
	double theta[data.Nth];

	compute_mat_product_vect( SIG, grad, tmp, 0.5*eps, data.Nth);

	for(i = 0; i<data.Nth; ++i) theta[i] = leader[i] + tmp[i];

	mvnrnd( theta, epsSIG, candidate, data.Nth);

	for (i = 0; i < data.Nth; ++i) {
		if (isnan(candidate[i])) {
			printf("!!!!  compute_manifold_candidate: isnan in candidate point (%d) and posdef =%d!\n", i, posdef);
			break;
		}

		if ((candidate[i] < data.lowerbound[i])||(candidate[i] > data.upperbound[i])) {
            runinfo.outside[runinfo.Gen] ++;
            return false;	// out of bounds
        }
	}
	return true; // all good
}


double TmcmcEngine::accept_ratio( double lnfo_lik, double lnfo_pri, double lnfc_lik, double lnfc_pri) 
{
    double L = exp( runinfo.p[runinfo.Gen]*(lnfc_lik-lnfo_lik) + (lnfc_pri-lnfo_pri));
    return L;
}


double TmcmcEngine::manifold_accept_ratio( double lnfo_lik, double lnfo_pri, double lnfc_lik, double lnfc_pri, double eps, double* thetao, double* thetac, double* gradiento, double* gradientc, double* SIGo)
{
    int Nth = data.Nth;
   	double p = runinfo.p[runinfo.Gen];
  	double tmpc[Nth];
  	double tmpo[Nth];

    if (data.options.Display > 2) {
        printf("generation: %d\n", runinfo.Gen);
        print_matrix("thetao",thetao,Nth);
        print_matrix("thetac",thetac,Nth);
        print_matrix("gradiento",gradiento,Nth);
        print_matrix("gradientc",gradientc,Nth);
        print_matrix("SIGo",SIGo,Nth, Nth);
    }

    double o, c;
	int i, k;
	for (i=0; i<Nth; ++i)
	{
		o = 0;
		c = 0;
		for(k=0; k<Nth; ++k)
		{
			o += SIGo[i*Nth+k]*gradiento[k];
			c += SIGo[i*Nth+k]*gradientc[k];
		}

		tmpc[i] = thetac[i] - thetao[i] - 0.5*eps*o;
		tmpo[i] = thetao[i] - thetac[i] - 0.5*eps*c;
	}

	double inv_SIGo[Nth*Nth];

	inv_matrix(SIGo, inv_SIGo, Nth);

   	double tmp2c[Nth];
	double tmp2o[Nth];

	compute_mat_product_vect( inv_SIGo, tmpc, tmp2c, 1.0, Nth );
	compute_mat_product_vect( inv_SIGo, tmpo, tmp2o, 1.0, Nth );

	double qc = -0.5*compute_dot_product(tmpc,tmp2c,Nth) / eps;
	double qo = -0.5*compute_dot_product(tmpo,tmp2o,Nth) / eps;

    //printf("qo: %f, qc: %f, lnfc_lik: %f, lnfo_lik: %f\n", qo, qc, lnfc_lik, lnfo_lik);

	double L = exp( p*(lnfc_lik-lnfo_lik) + (lnfc_pri-lnfo_pri) + qo-qc );  
    
    return L;
}


void TmcmcEngine::manifold_calculate_grad(const double* grad, double* gradOut)
{
    for (int i = 0; i < data.Nth; ++i) gradOut[i] = runinfo.p[runinfo.Gen]*grad[i];
    return;
}


int TmcmcEngine::manifold_calculate_Sig(double *pSIGMA, bool posdef, double eval[], double evec[], const double* pos)
{
    int i, j, l;
	int Nth = data.Nth;
    double p = runinfo.p[runinfo.Gen];

	gsl_vector *gsl_eval = gsl_vector_alloc(Nth);
    std::copy(eval, eval+Nth, gsl_eval->data);

	// correct negative eigenvalues
	if (!posdef)
	{
		gsl_permutation *idx = gsl_permutation_alloc(data.Nth);
		gsl_sort_vector_index(idx, gsl_eval);

		/* TODO: move this part in main, calculating eval and evecs of sample cov must only be calculated once */
		gsl_matrix *SS = gsl_matrix_alloc(data.Nth, data.Nth);
		for (i = 0; i < Nth; ++i)
			for (j = 0; j < Nth; ++j)
				gsl_matrix_set(SS, i, j, runinfo.SS[i][j]);

		gsl_vector *evalSS = gsl_vector_alloc(Nth);
		gsl_matrix *evecSS = gsl_matrix_alloc(Nth, Nth);
		gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(Nth);

		//evaluate and sort eigenvalues of SS
		gsl_eigen_symmv(SS, evalSS, evecSS, w);
		gsl_eigen_symmv_sort(evalSS, evecSS, GSL_EIGEN_SORT_VAL_ASC);

		gsl_eigen_symmv_free(w);

		// replacing negative eigenvalues by evals from sample cov (ascending order)
		for (i=0; i<Nth; ++i)
		{
			if ( gsl_vector_get(gsl_eval, idx->data[i])<=0 )  
                //smallest eigenvalue of SS
				gsl_vector_set(gsl_eval, idx->data[i], data.beta2*evalSS->data[0]);
			else  
                //scale positive eigs
				gsl_vector_set( gsl_eval, idx->data[i], gsl_vector_get(gsl_eval, idx->data[i])/p );
		}
		gsl_vector_free(evalSS);
		gsl_matrix_free(evecSS);
		gsl_permutation_free(idx);
		gsl_matrix_free(SS);
	}
	else
	{
		for(int i = 0; i<Nth; ++i)
            //scale all positive eigs
            gsl_vector_set(gsl_eval, i, gsl_vector_get(gsl_eval, i)/p);
	}

	// make eval adaption to ebds (if necessary)
	// const double chi2 = gsl_cdf_chisq_Qinv(data.conf,data.Nth);	//TODO: check chi2 value if chisq_Pinv is correct choice (chi2 = 9.2704 for Nth = 8)
    
    gsl_matrix *out_lik_evec = gsl_matrix_alloc (Nth, Nth);
    std::copy(evec, evec+Nth*Nth, out_lik_evec->data);
    
    if (data.moptions.adapt)
    {

        gsl_vector *eigv = gsl_vector_alloc (Nth);

        for (l=0; l<Nth; ++l)
        {
            double scOrig, sc;
            scOrig = sqrt(gsl_vector_get(gsl_eval,l)*data.moptions.chi2);
            // correct eigenvectors to extended bounds
            gsl_matrix_get_col (eigv, out_lik_evec, l);
            sc = scale_to_box(pos,  scOrig, eigv->data, data.moptions.elbds, data.moptions.eubds, Nth);
            sc = scale_to_box(pos, -sc, eigv->data, data.moptions.elbds, data.moptions.eubds, Nth);
            gsl_vector_set(gsl_eval, l, sc*sc/data.moptions.chi2);
            if (scOrig != sc) runinfo.corrections[runinfo.Gen]++;
        }
        gsl_vector_free(eigv);
    }

	gsl_matrix *EVD_mat = gsl_matrix_alloc(Nth, Nth);
	for (i = 0; i < Nth; ++i)
    		for (j = 0; j < Nth; ++j)
      			gsl_matrix_set(EVD_mat, i, j, gsl_matrix_get(out_lik_evec,i,j)*gsl_vector_get(gsl_eval, j));

	gsl_matrix *SIGMA = gsl_matrix_alloc(Nth, Nth);
	gsl_matrix_set_zero(SIGMA);
	int err_gemm = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, EVD_mat, out_lik_evec, 1.0, SIGMA);
    if (err_gemm) printf("manifold_calculate_Sig: gsl_blas_dgemm failed (case untreated) !!!\n");


    // TEST AFTER ADAPTION
    int err_flg = 0;
    if (data.moptions.adapt)
    {
	    gsl_matrix *cSIGMA = gsl_matrix_alloc(Nth, Nth);
        std::copy(SIGMA->data, SIGMA->data+Nth*Nth,  cSIGMA->data); // XXX is this OK? vector-matrix
        
        err_flg = gsl_linalg_cholesky_decomp(cSIGMA);
        
        if (err_flg == GSL_EDOM) {
            printf("manifold_calculate_Sig: correcton failed, SIGMA not posdef!!!\n");
            runinfo.failedcorrections[runinfo.Gen]++;
        } else{
            std::copy(SIGMA->data, SIGMA->data+Nth*Nth,  pSIGMA); // XXX is this OK? vector-matrix
        }

	    gsl_matrix_free(cSIGMA);
    
    }
    else
    {
        std::copy(SIGMA->data, SIGMA->data+Nth*Nth,  pSIGMA); // XXX is this OK? vector-matrix
    }

	// clean memory
	gsl_matrix_free(EVD_mat);
	gsl_matrix_free(SIGMA);
	gsl_matrix_free(out_lik_evec);
	gsl_vector_free(gsl_eval);

	return err_flg;
}



void TmcmcEngine::chaintask(double in_tparam[], int pnsteps, double out_tparam, 
                            int winfo[4], double *init_mean, double *chain_cov)
{

    int gen_id   = winfo[0];
    int chain_id = winfo[1];

    long me = torc_worker_id();

    double leader[data.Nth], loglik_leader, logprior_leader;            /* old*/
    double candidate[data.Nth], loglik_candidate, logprior_candidate;   /* new*/

    // get initial leader and its value
    for (int i = 0; i < data.Nth; ++i) leader[i] = in_tparam[i];
    loglik_leader   = out_tparam;
    logprior_leader = prior.eval_logpdf(leader);

    double pj = runinfo.p[runinfo.Gen];

    int burn_in = data.burn_in;

    for (int step = 0; step < pnsteps + burn_in; ++step) {
        double chain_mean[data.Nth];
        if (step == 0)
            for (int i = 0; i < data.Nth; ++i) chain_mean[i] = init_mean[i];
        else
            for (int i = 0; i < data.Nth; ++i) chain_mean[i] = leader[i];

        bool candidate_inbds = compute_candidate(candidate, chain_mean); // I keep this for the moment, for performance reasons

        if (candidate_inbds) {

            evaluate_candidate(candidate, &loglik_candidate, me, gen_id, chain_id, step, 1);    // this can spawn many tasks

            if (data.ifdump && step >= burn_in) torc_update_full_db(candidate, loglik_candidate, nullptr, 0, 0);
            // last argument should be 1 if it is a surrogate

            logprior_candidate = prior.eval_logpdf(candidate);
            double L = exp((logprior_candidate-logprior_leader)+(loglik_candidate-loglik_leader)*pj);

            double P = 0;
            if (L > 1) L = 1;
            else P = uniformrand(0,1);
            if (P < L) {
                /* accept new leader */
                for (int i = 0; i < data.Nth; ++i) leader[i] = candidate[i];
                loglik_leader = loglik_candidate;
        
            }
        }
            
        /* increase counter or add the leader again in curgen_db */
        if (step >= burn_in) {
            logprior_leader = prior.eval_logpdf(leader);
            torc_update_curgen_db(leader, loglik_leader, logprior_leader);
        }
    }

    return;
}


//TODO: remove pointers where not needed
void TmcmcEngine::manifold_chaintask(double in_tparam[], int pnsteps, double out_tparam, 
                                    int t_err, bool t_posdef, double *t_grad, 
                                    double *t_cov, double *t_evec, double *t_eval, int winfo[4]) 
{
    int gen_id   = winfo[0];
    int chain_id = winfo[1];

    long me = torc_worker_id();

    double leader[data.Nth], loglik_leader;             /* old */
    double candidate[data.Nth], loglike_candidate;      /* new */

    int l_err;
    bool l_posdef;
    double l_grad[data.Nth], l_grad_scaled[data.Nth];
	double l_cov[data.Nth*data.Nth];
	double l_evec[data.Nth*data.Nth];
	double l_eval[data.Nth];
    
    int c_err;
    bool c_posdef;
    double c_grad[data.Nth], c_grad_scaled[data.Nth];
    double c_cov[data.Nth*data.Nth];
    double c_evec[data.Nth*data.Nth];
    double c_eval[data.Nth];

    double logprior_leader, logprior_candidate;

    /* get initial leader and its values */
    int i, step;
	for (i = 0; i < data.Nth; i++) leader[i] = in_tparam[i];
	loglik_leader = out_tparam;

	l_err    = t_err;
	l_posdef = t_posdef;
	for ( i=0; i < data.Nth; i++)			l_grad[i] = t_grad[i];
	for ( i=0; i < data.Nth*data.Nth; i++)	l_cov[i]  = t_cov[i];
	for ( i=0; i < data.Nth*data.Nth; i++)	l_evec[i] = t_evec[i];
	for ( i=0; i < data.Nth; i++)			l_eval[i] = t_eval[i];

    int burn_in = data.burn_in;

    /* sigma to be calculated */
    double l_SIG[data.Nth*data.Nth];
    
    for (step = 0; step < pnsteps + burn_in; ++step) {
        bool candidate_inbds = false;
        int s_err = 1;

        if (l_err == 0) {
                manifold_calculate_grad(l_grad, l_grad_scaled);
                s_err = manifold_calculate_Sig(l_SIG, l_posdef, l_eval, l_evec, leader);

                if (s_err == 0)
				    candidate_inbds = compute_manifold_candidate(candidate, leader, data.moptions.eps, l_SIG, l_grad_scaled, l_posdef);
			    else
				    candidate_inbds = compute_candidate(candidate, leader);

        } else {
        			candidate_inbds = compute_candidate(candidate, leader);
        }

        if (!candidate_inbds) continue; /* go to next step */

        /* candidate in bounds */
        manifold_evaluate_candidate(candidate, &loglike_candidate, &c_err, &c_posdef, c_grad, c_cov, c_evec, c_eval, me, gen_id, chain_id, step, 1);

        if (data.ifdump && step >= burn_in) torc_update_full_db(candidate, loglike_candidate, nullptr, 0, 0); //TODO: not the manifold update version? (DW)
    
    	logprior_candidate = prior.eval_logpdf(candidate);
		logprior_leader    = prior.eval_logpdf(leader);
    
        double L = 0;
        if( l_err > 0)
            L = accept_ratio(loglik_leader, logprior_leader, loglike_candidate, logprior_candidate);
        else 
        {
            if( c_err == 0 ){ 
                /* BASIS */
				if( s_err ==0 ){
                    /* mMALA */
					manifold_calculate_grad(c_grad, c_grad_scaled);
					L = manifold_accept_ratio(loglik_leader, logprior_leader, loglike_candidate, logprior_candidate, data.moptions.eps, leader, candidate, l_grad_scaled, c_grad_scaled, l_SIG);
				}
				else if( s_err > 0 )
					L = accept_ratio(loglik_leader, logprior_leader, loglike_candidate, logprior_candidate);
			}
			else if(c_err == 1){ 
                    /* reject (ODE failed && grad not available) */
                    L = 0.0;
			}
			else if(c_err == 2){ 
                    /* mMALA */
                    manifold_calculate_grad(c_grad, c_grad_scaled);
                    L = manifold_accept_ratio(loglik_leader, logprior_leader, loglike_candidate, logprior_candidate, data.moptions.eps, leader, candidate, l_grad_scaled, c_grad_scaled, l_SIG);
			}
        }

        double P = uniformrand(0,1);

        if (P < L) {
            /* accept new leader */
			for (i = 0; i < data.Nth; i++) leader[i] = candidate[i];

			/* update everything */
			loglik_leader = loglike_candidate;
			l_err = c_err;
			l_posdef = c_posdef;
			for (i = 0; i < data.Nth; i++) l_grad[i] = c_grad[i];
			for (i = 0; i < data.Nth; i++) l_grad_scaled[i] = c_grad_scaled[i];

			if( c_err == 0 ){
				for (i = 0; i < data.Nth*data.Nth; i++)	l_cov[i]  = c_cov[i];
				for (i = 0; i < data.Nth*data.Nth; i++)	l_evec[i] = c_evec[i];
				for (i = 0; i < data.Nth; i++)			l_eval[i] = c_eval[i];
			}

		}        
 
        /* increase counter or add the leader again in curgen_db*/
        if (step >= burn_in) {
            logprior_leader = prior.eval_logpdf(leader);
            torc_update_manifold_curgen_db(leader, loglik_leader, logprior_leader, l_err, l_posdef, l_grad, l_cov, l_evec, l_eval);
        }
    }
    
    return;

}


void TmcmcEngine::prepare_newgen(int nchains, cgdbp_t *leaders)
{
    /* process curgen_db -> calculate statitics */
    /* compute probs based on F values */
    /* draw new samples (nchains or user-specified) */
    /* find unique samples: fill the (new) leaders table */
    /* count how many times they appear -> nsteps */
    /* return the new sample size (number of chains) */

    int i, p;

    int n = curgen_db.entries;


    double **g_x = new double*[data.Nth];
    for (i = 0; i < data.Nth; ++i) g_x[i] = new double[n];

    
    /* calculate uniques & acceptance rate */
    {
        double **uniques = g_x;
        int un = 0,  j;
        bool unique_flag;

        /* copy first into uniques list */
        for( p = 0; p < data.Nth; ++p ) uniques[p][un] = curgen_db.entry[0].point[p];
        un++;

        for (i = 1; i < n; ++i) {
            double xi[data.Nth];
            for (p = 0; p < data.Nth; ++p) xi[p] = curgen_db.entry[i].point[p];

            for (j = 0; j < un; ++j) {  /* compare with  previous uniques */
                unique_flag = false;             /* is this point unique? */
                for (p = 0; p < data.Nth; ++p) {

                    /* do they differ in position? */
                    if (fabs(xi[p]-uniques[p][j]) > 1e-8) {
                        unique_flag = true; /* unique */
                        break;              /* check next */
                    }

                }
                if (unique_flag == false) break;
            }

            if (unique_flag) { /* unique, put it in the table */
                for (p = 0; p < data.Nth; ++p) uniques[p][un] = xi[p];
                un++;
            }
        }

        runinfo.currentuniques[runinfo.Gen] = un;

        runinfo.acceptance[runinfo.Gen] = (1.0*runinfo.currentuniques[runinfo.Gen])/data.Num[runinfo.Gen];


        printf("prepare_newgen: num uniques directly after :  %d\n", un);
    } 
    /* end block*/

    unsigned int *sel = new unsigned int[n];
    {
        /* calculate statistics and resample */
        double *fj = new double[n];
        double t0  = torc_gettime();
        for (i = 0; i < n; ++i) fj[i] = curgen_db.entry[i].F;    /* separate point from F ?*/
        double t1 = torc_gettime();
        calculate_statistics(fj, data.Num[runinfo.Gen], runinfo.Gen, sel);
        double t2 = torc_gettime();
        printf("prepare_newgen: init + calc stats : %lf + %lf = %lf seconds\n", t2-t1, t1-t0, t2-t0);
        delete[] fj;
    }

    
    /* calculate new chains */
    int totsel = 0;
    int newchains = 0;
    std::vector<sort_t> list(n);
    for (i = 0; i < n; ++i) {
        list[i].idx  = i;
        list[i].nsel = sel[i];
        if (sel[i] != 0) newchains++;
        totsel += list[i].nsel;
    }
    printf("prepare_newgen: newchains sampled :  %d\n", newchains);
    printf("prepare_newgen: total selections after sampling:  %d (sanity check) \n", totsel);

    std::sort(list.begin(), list.end(), compar_desc);

    if (data.MaxChainLength > 0) {
        /* UPPER THRESHOLD */
        /* splitting long chains */
        totsel = 0;
        int initial_newchains = newchains;
        for (i = 0; i < initial_newchains; ++i) {
            //printf("before list[i].nsel %d (totsel: %d) \n ", list[i].nsel, totsel);
            while (list[i].nsel > data.MaxChainLength) {
                list.push_back( { .idx = list[i].idx, .nsel = data.MaxChainLength } );
                totsel += data.MaxChainLength;
                newchains++;
                list[i].nsel = list[i].nsel - data.MaxChainLength;
            }
            totsel += list[i].nsel;
            //printf("after list[i].nsel %d (%d) \n", list[i].nsel, totsel);
        }
        printf("prepare_newgen: newchains after breaking long chains :  %d\n", newchains);
        printf("prepare_newgen: total selections after breaking long chains:  %d (sanity check) \n", totsel);
        std::sort(list.begin(), list.end(), compar_desc);
        //for(auto& el : list) printf("  (idx:%d, nsel:%d)  \n", el.idx, el.nsel);
    }

    if (data.MinChainLength > 0) {
        /* LOWER THRESHOLD */
        /* setting min chain length */
        //TODO: untested feature (DW)
        totsel = 0;
        int l_threshold = data.MinChainLength;
        for (i = 0; i < newchains; ++i) {
            if ((list[i].nsel > 0)&&(list[i].nsel < l_threshold)) {
                list[i].nsel = l_threshold;
            }
            totsel += list[i].nsel;
        }
        printf("prepare_newgen: newchains after upstepping short chains :  %d (sanity check) \n", newchains);
        printf("prepare_newgen: total selections after upstepping short chains:  %d\n", totsel);
        std::sort(list.begin(), list.end(), compar_desc);
    }

    /* TODO: do we need to copy this? (DW) esp. prior?
    int counter;    
    int surrogate;      
    double error;       
    */

    for (i = 0; i < newchains; ++i) {       /* newleader */
            int idx = list[i].idx;

            for (p = 0; p < data.Nth ; p++) leaders[i].point[p] = curgen_db.entry[idx].point[p];
    
            if( _method == Manifold ) {
                leaders[i].error_flg = curgen_db.entry[idx].error_flg;
                leaders[i].posdef = curgen_db.entry[idx].posdef;
                for (p = 0; p < data.Nth ; p++)          leaders[i].gradient[p] = curgen_db.entry[idx].gradient[p];
                for (p = 0; p < data.Nth ; p++)          leaders[i].eval[p] = curgen_db.entry[idx].eval[p];                
                for (p = 0; p < data.Nth*data.Nth ; p++) leaders[i].cov[p] = curgen_db.entry[idx].cov[p];
                for (p = 0; p < data.Nth*data.Nth ; p++) leaders[i].evec[p] = curgen_db.entry[idx].evec[p];
            }
        
            leaders[i].F     = curgen_db.entry[idx].F;
            leaders[i].prior = curgen_db.entry[idx].prior;
            leaders[i].nsel  = list[i].nsel;
    }

    printf("prepare_newgen: CURGEN DB (Gen: %d): [newchains=%d]\n", runinfo.Gen, newchains);

    if (data.use_local_cov)
        precompute_chain_covariances(leaders, data.init_mean, data.local_cov, newchains);

    curgen_db.entries = 0;

    for (i = 0; i < data.Nth; ++i) delete [] g_x[i];

    delete[] g_x;
    delete[] sel;

    this->nchains = newchains;
    return;
}

void TmcmcEngine::print_runinfo()
{
    printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    printf("runinfo.Gen = \n\n   %d\n\n", runinfo.Gen);
    print_matrix("runinfo.p", runinfo.p, runinfo.Gen+1);
    print_matrix("runinfo.mean_of_theta", runinfo.meantheta[runinfo.Gen], data.Nth);
    print_matrix_2d("runinfo.SS", runinfo.SS, data.Nth, data.Nth);
    print_matrixi("runinfo.outside", runinfo.outside, runinfo.Gen+1);
    if ( _method == Manifold ) {
        print_matrixi("runinfo.corrections", runinfo.corrections, runinfo.Gen+1);
        print_matrixi("runinfo.failedcorrections", runinfo.failedcorrections, runinfo.Gen+1);
    }
    printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");

    printf("Acceptance rate         :  %lf \n", runinfo.acceptance[runinfo.Gen]) ;
    printf("Annealing exponent      :  %lf \n", runinfo.p[runinfo.Gen]) ;
    printf("Coeficient of Variation :  %lf \n", runinfo.CoefVar[runinfo.Gen]) ;
    printf("----------------------------------------------------------------\n");

}

void TmcmcEngine::bcast_runinfo()
{
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


double TmcmcEngine::getLogEvidence() const
{
    return compute_sum(runinfo.logselections, runinfo.Gen+1);
}

double * TmcmcEngine::getNewMean() const
{
    double * newmean = new double[data.Nth];
    for(int i = 0; i < data.Nth; ++i) { newmean[i] =  runinfo.meantheta[runinfo.Gen][i]; }
    return newmean;
}

double** TmcmcEngine::getNewSampleCov() const
{
    double ** newSS = new double*[data.Nth];
    for(int i = 0; i < data.Nth; ++i) { 
        newSS[i] = new double[data.Nth];  
        for(int j = 0; j < data.Nth; ++j)
            newSS[i][j] = runinfo.SS[i][j]; 
    }
    return newSS;
}

} //namespace tmcmc
