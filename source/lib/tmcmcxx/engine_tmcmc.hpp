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

extern "C" {
#include "priors.h"
}
#include "tmcmc_types.hpp"
#include "tmcmc_utils.hpp"


namespace tmcmc {

    class TmcmcEngine {

    public:
        TmcmcEngine();
        ~TmcmcEngine();

        void run();
    private:
        data_t    data;
        runinfo_t runinfo;

        db_t    full_db;
        cgdb_t  curgen_db;
        resdb_t curres_db;

        Density *priors;

        int EXPERIMENTAL_RESULTS = 0; //TODO: where does this belong? (DW)

        void init(); // TODO: called data_init (reminder)

        void init_full_db();
        //void update_full_db();
        void update_full_db(double point[], double F, double *G, int n, int surrogate);
        void torc_update_full_db(double point[], double F, double *G, int n, int surrogate);
        void print_full_db();
        void dump_full_db(); 

        void init_curgen_db();      // (TODO: check if different funcs are needed for
                                   //   full, curgen and curres)

        //void update_curgen_db();
        void update_curgen_db(double point[], double F, double prior);
        void torc_update_curgen_db_task(double point[], double *pF, double *pprior);
        void torc_update_curgen_db(double point[], double F, double prior);
        //void print_curgen_db();  
        void dump_curgen_db(int gen);
        //void display_curgen_db(int gen);
        int load_curgen_db(int gen);

        void init_curres_db();
    //    void update_curres_db();
        void update_curres_db(double point[/* EXPERIMENTAL_RESULTS */], double F);
        void dump_curres_db(int gen);

        void calculate_statistics(double flc[], unsigned int n, int nselections, 
                                 int gen, unsigned int sel[]);


        void taskfun(const double *x, int *pN, double *res, int winfo[4]);
        void evaluate_F(double point[], double *Fval, int worker_id, int gen_id, 
                            int chain_id, int step_id, int ntasks);
        void initchaintask(double in_tparam[], int *pdim, double *out_tparam, int winfo[4]);
        void check_for_exit();

        void precompute_chain_covariances(const cgdbp_t* leader,double** init_mean, 
                                            double** chain_cov, int newchains); // TODO: chain_cov, mean etc can be rmvd from args?

        int compute_candidate(double candidate[], double chain_mean[], double var);
        int compute_candidate_cov(double candidate[], double chain_mean[], double chain_cov[]);
        void chaintask(double in_tparam[], int *pdim, int *pnsteps, double *out_tparam, int winfo[4],
                        double *init_mean, double *chain_cov);
        int prepare_newgen(int nchains, cgdbp_t *leaders);
    };

} //namespace tmcmc

#endif //ENGINE_TMCMC_HPP
