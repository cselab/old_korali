/* -------------------------------------------------------------------------- */
/* --- File: engine_tmcmc.hpp--- Author: Daniel Waelchli -------------------- */
/* ---------------------- last modified: Jan 2019        -------------------- */
/* --------------------------------- by: Daniel Waelchli -------------------- */
/* -------------------------------------------------------------------------- */
/*
     Cpp-Class of Transitional Markov Chain Monte Carlo for sampling of
        posterior distributions of model parameter.

	 Implementation based on work of Panagiotis Hadjidoukas

*/

#ifndef ENGINE_TMCMC_HPP
#define ENGINE_TMCMC_HPP

#include "myrand.hpp"
#include "Ifitfun.hpp"
#include "tmcmc_utils.hpp"

namespace tmcmc
{

class TmcmcEngine
{

public:
    
    TmcmcEngine(fitfun::IFitfun * Ifitfun_ptr, 
                Method method = Standard, 
                const char * ftmcmcpar = "tmcmc.par",
                const char * fpriorpar = "priors.par");
    ~TmcmcEngine();
    void run();

private:
    
    Method _method;

    data_t    data;
    runinfo_t runinfo;

    int nchains;

    db_t    full_db;
    cgdb_t  curgen_db;
    resdb_t curres_db;

    double t0;
    double gt0;
    double gt1;

    double *out_tparam;
    cgdbp_t *leaders;

    priors::Prior prior;

    fitfun::IFitfun * ifitfun_ptr_;

    void init();
    void evalGen();

    bool load_data();
    void sample_from_prior();

    void init_full_db();
    void update_full_db(double point[], double F, double *G, int n, int surrogate);
    void torc_update_full_db(double point[], double F, double *G, int n, int surrogate);
    void print_full_db();
    void dump_full_db();

    void init_curgen_db();

    int load_curgen_db();
    void update_curgen_db(double point[], double F, double prior);
    void torc_update_curgen_db_task(double point[], double *pF, double *pprior);
    void torc_update_curgen_db(double point[], double F, double prior);
    void dump_curgen_db();

    void update_manifold_curgen_db(
            double point[], double F, double prior, 
            int error_flg, bool posdef, double gradient[], double cov[], 
            double evec[], double eval[]);
    void torc_update_manifold_curgen_db(
            double point[], double F, double prior, 
            int error_flg, bool posdef, double gradient[], double cov[], 
            double evec[], double eval[]);
    
    void init_curres_db();
    void dump_curres_db(int gen);

    void calculate_statistics(double flc[], int nselections, int gen,
                              unsigned int sel[]);


    void evaluate_candidate(double point[], double *Fval, int worker_id, int gen_id,
                    int chain_id, int step_id, int ntasks);
    
    void manifold_evaluate_candidate(double point[], double *Fval, int *perr, 
                             bool *pposdef, double grad[], double cov[], 
                             double evec[], double eval[], int worker_id, 
                             int gen_id, int chain_id, int step_id, int ntasks);
    
    void initchaintask(double in_tparam[], int winfo[4]);
    void check_for_exit();

    void precompute_chain_covariances(const cgdbp_t* leader,double** init_mean,
                                      double** chain_cov, int newchains);

    double accept_ratio( double lnfo_lik, double lnfo_pri, double lnfc_lik, double lnfc_pri);
    double manifold_accept_ratio( double lnfo_lik, double lnfo_pri, double lnfc_lik, double lnfc_pri, double eps, double* thetao, double* thetac, double* gradiento, double* gradientc, double* SIGo);
   
    void manifold_calculate_grad(const double* grad, double* gradOut);
    int manifold_calculate_Sig(double *pSIGMA, bool posdef, double eval[], double evec[], const double* pos);
    
    bool compute_candidate(double candidate[], double chain_mean[]);
    bool compute_candidate_cov(double candidate[], double chain_mean[], double chain_cov[]);
    bool compute_manifold_candidate(double candidate[], double leader[], double eps, double *SIG, double *grad, bool posdef);
    void chaintask(double in_tparam[], int pnsteps, double out_tparam, int winfo[4],
                   double *init_mean, double *chain_cov);
    void manifold_chaintask(double in_tparam[], int pnsteps, double out_tparam, 
                            int t_err, bool t_posdef, double *t_grad, 
                            double *t_cov, double *t_evec, double *t_eval, int winfo[4]);


    int prepare_newgen(int nchains, cgdbp_t *leaders);


    void print_runinfo();
    void bcast_runinfo();
    void spmd_update_runinfo();
};

} //namespace tmcmc

#endif //ENGINE_TMCMC_HPP
