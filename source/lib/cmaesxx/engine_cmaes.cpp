#include <stdexcept>

#include "engine_cmaes.hpp"
#include "engine_cmaes_utils.hpp"

CmaesEngine::CmaesEngine(double (*fun) (double*, int, void*, int*), 
	std::string workdir, std::string cmaes_par, 
	std::string cmaes_bounds_par, std::string priors_par, 
	int restart) 
		: workdir_(workdir), 
		  cmaes_par_(cmaes_par),
		  cmaes_bounds_par_(cmaes_bounds_par),
		  priors_par_(priors_par),
		  restart_(restart),
		  step_(0), 
		  stt_(0.0) {

		if (!GetCurrentDir(exeDir_, sizeof(exeDir_))) {
			throw std::runtime_error("Current directory could not be localized.");
		}

		if (ChangeDir(workdir_.c_str()) != 0) {
			throw std::runtime_error("Working directory '" + workdir_ + "' does not exist.");
		}

		if ( !cmaes_utils_file_exists(cmaes_par.c_str()) ) { 
			throw std::runtime_error("Cmaes param file '" + cmaes_par_ + "' does not exist.");
		}

		if ( !cmaes_utils_file_exists(cmaes_bounds_par.c_str()) ) {
			throw std::runtime_error("Cmaes bounds param file '" + cmaes_bounds_par_ + "' does not exist.");
		}
		if ( !cmaes_utils_file_exists(priors_par.c_str()) ) { 
			throw std::runtime_error("Priors param file '" + priors_par_ + "' does not exist.");
		}

		gt0_ = get_time();
		
		CmaesEngine::fitfun_ = fun;
		
        arFunvals_ = cmaes_init(&evo_, 0, NULL, NULL, 0, 0,
						cmaes_par_.c_str());
		
		printf("%s\n", cmaes_SayHello(&evo_));
		cmaes_ReadSignals(&evo_, cmaes_par_.c_str()); 

		dim_    = cmaes_Get(&evo_, "dim");
		lambda_ = cmaes_Get(&evo_, "lambda");

		cmaes_utils_read_bounds(VERBOSE, cmaes_bounds_par_.c_str(), 
			&lower_bound_, &upper_bound_, dim_);

		int Nprior;
		read_priors( priors_par_.c_str(), &priors_, &Nprior );
	
		if( Nprior != dim_ ){
			throw std::runtime_error("The dimension of the prior is different from the dimension of the fitfun.");
		}

#if defined(_USE_TORC_)
	    torc_register_task(CmaesEngine::taskfun_);
        int argc = 0; //dummy (TODO)
        char **argv;  //dummy (TODO)
		torc_init(argc, argv, MODE_MS); 
#endif

		if (ChangeDir(exeDir_) != 0) {
			std::runtime_error("Could not return to exe dir.");
		}

}

CmaesEngine::~CmaesEngine(){
    cmaes_exit(&evo_); /* release memory */

/* (TODO: do we still need this?)
#if defined(_USE_TORC_)
    torc_finalize();
#endif
*/
}

cmaes_t* CmaesEngine::getEvo() {
	return &evo_;
}

double* CmaesEngine::getBestEver() {
	return evo_.rgxbestever;
}

double CmaesEngine::getBestFunVal() {
	return cmaes_Get(&evo_,"fbestever");
}

double CmaesEngine::evaluate_population( cmaes_t *evo, double *arFunvals, 
	double *const* pop, Density *d, int step ) {

    int info[4];
    double tt0, tt1 ;
    	
    tt0 = get_time();
	
    for( int i = 0; i < lambda_; ++i){
        info[0] = 0; info[1] = 0; info[2] = step; info[3] = i;     /* gen, chain, step, task */
		
#if defined(_USE_TORC_)

        torc_create( -1, CmaesEngine::taskfun_, 4,
                     dim_, MPI_DOUBLE, CALL_BY_VAL,
                     1, MPI_INT, CALL_BY_COP,
                     1, MPI_DOUBLE, CALL_BY_RES,
                     4, MPI_INT, CALL_BY_COP,
                     pop_[i], &dim_, &arFunvals_[i], info);

#else
        CmaesEngine::taskfun_(pop_[i], &dim_, &arFunvals_[i], info);
#endif
    }
	
#if defined(_USE_TORC_)

#if defined(_STEALING_)
    torc_enable_stealing();
#endif

    torc_waitall();

#if defined(_STEALING_)
    torc_disable_stealing();
#endif

#endif

    // subtract the log-prior from the log-likelohood
    for( int i=0; i<lambda_; i++){
        arFunvals_[i] -= prior_log_pdf(priors_, dim_, pop_[i]);
    }

    tt1 = get_time();
  
    return tt1-tt0;
};


double CmaesEngine::run() {

	gt1_ = get_time();

	if (ChangeDir(workdir_.c_str()) != 0) {
		printf("Working directory '%s' does not exist. \
			Exit with exit(1)...\n", workdir_.c_str());
		exit(1);
	}

	double dt;

	while( !cmaes_TestForTermination(&evo_) ){

        pop_ = cmaes_SamplePopulation(&evo_); 
		
        if (restart_) {
            dt = cmaes_utils_load_pop_from_file(VERBOSE, step_, pop_, 
					arFunvals_, dim_, lambda_, &restart_);
    	} else {
			cmaes_utils_make_all_points_feasible( &evo_, pop_
				, lower_bound_, upper_bound_);
            dt = evaluate_population( &evo_, arFunvals_, pop_, priors_, step_);
        }
        stt_ += dt;
	
        cmaes_UpdateDistribution(1, &evo_, arFunvals_);

        cmaes_ReadSignals(&evo_, cmaes_par_.c_str()); fflush(stdout);

        if (VERBOSE) cmaes_utils_print_the_best(evo_, step_);
		
       	if (_IODUMP_ && !restart_){
            cmaes_utils_write_pop_to_file(evo_, arFunvals_, pop_, step_);
        }

#if defined(_RESTART_)
        cmaes_WriteToFile(&evo_, "resume", "allresumes.dat");
#endif

        if( ! cmaes_utils_is_there_enough_time( JOBMAXTIME, gt0_, dt ) ){
            evo_.sp.stopMaxIter=step_+1;
            break;
        }
        
        step_++;
    }

    gt2_ = get_time();

    printf("Stop:\n %s \n",  cmaes_TestForTermination(&evo_)); /* print termination reason */
    cmaes_WriteToFile( &evo_, "all", "allcmaes.dat" );         /* write final results */

    gt3_ = get_time();
    
    printf("Total elapsed time      = %.3lf  seconds\n", gt3_-gt0_);
    printf("Initialization time     = %.3lf  seconds\n", gt1_-gt0_);
    printf("Processing time         = %.3lf  seconds\n", gt2_-gt1_);
    printf("Funtion Evaluation time = %.3lf  seconds\n", stt_);
    printf("Finalization time       = %.3lf  seconds\n", gt3_-gt2_);

    if (ChangeDir(exeDir_) != 0) {
		printf("Could not return to exe dir '%s'. \
			Exit with exit(1)...\n", exeDir_);
		exit(1);
	}
	return 0.0;
}

double (*CmaesEngine::fitfun_) (double*, int, void*, int*);

void CmaesEngine::taskfun_(double *x, int *n, double *res, 
	int /* deprecated */ *info) {
	
    (*res) = - CmaesEngine::fitfun_(x, *n, (void*)NULL, info);    // minus for minimization
	
	return;
}
