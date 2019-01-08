/* --------------------------------------------------------- */
/* --- File: engine_cmaes.cpp--- Author: Daniel Waelchli --- */
/* ---------------------- last modified: Jan 2019        --- */
/* --------------------------------- by: Daniel Waelchli --- */
/* --------------------------------------------------------- */
/*   
     Cpp Wrapper for parallel CMA-ES for non-linear function minimization. 
	
	 Implementation based on work of Panagiotis Hadjidoukas (see engine_cmaes.c) 
		and Nikolaus Hansen (see cmaes.c)
      
*/

#include <string>

extern "C" {

#include <cmaes.h>
#include <cmaes_utils.h>
#include <priors.h>

}


#if defined(_USE_TORC_)

#include <mpi.h>
#include <torc.h>

#endif

#define VERBOSE 1
#define JOBMAXTIME 0
#define _IODUMP_ 1

class CmaesEngine {

public:
	CmaesEngine(double (*fitfun) (double*, int), 
		std::string cmaes_par, std::string cmaes_bounds_par, 
		std::string prios_par, int restart = 0); 
	
	double run();
	double evaluate_population( cmaes_t *evo, double *arFunvals, double * const* pop, Density *d, int step );

private:

	std::string cmaes_par_, cmaes_bounds_par_, priors_par_;

	int restart_;
	cmaes_t evo_;
	int lambda_;
	int step_;
	
	double gt0_, gt1_, gt2_, gt3_;
	double stt_;
	
	int dim_;
	double *lower_bound_, *upper_bound_;

	Density *priors_;

	double *const*pop_;
    double *arFunvals_; 

	double (*fitfun_) (double*, int);
	void taskfun_(double *x, int *no, double* res, int *info);

};



