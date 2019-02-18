/* -------------------------------------------------------------------------- */
/* --- File: engine_cmaes.cpp--- Author: Daniel Waelchli -------------------- */
/* ---------------------- last modified: Jan 2019        -------------------- */
/* --------------------------------- by: Daniel Waelchli -------------------- */
/* -------------------------------------------------------------------------- */
/*   
     Cpp Wrapper for parallel CMA-ES for gradient free non-linear 
     	function minimization. 
	
	 Implementation based on work of Panagiotis Hadjidoukas (see engine_cmaes.c) 
		and Nikolaus Hansen (see cmaes.c)
      
*/

#ifndef ENGINE_CMAES_HPP
#define ENGINE_CMAES_HPP

#include <string>
#include <stdio.h>

extern "C" {

#include <cmaes.h>
#include <cmaes_utils.h>
#include <priors.h>

}


#if defined(_USE_TORC_)

extern "C" {
    #include <torc.h>
}

#endif

#define VERBOSE 0
#define JOBMAXTIME 0
#define _IODUMP_ 1

class CmaesEngine {

public:
	CmaesEngine(double (*fun) (double*, int, void*, int*), 
		std::string workdir = ".", 
		std::string cmaes_par = "cmaes_initials.par", 
		std::string cmaes_bounds_par = "cmaes_bounds.par", 
		std::string prios_par = "priors.par", 
		int restart = 0); 

	~CmaesEngine();

	double run();

	cmaes_t* getEvo();
	double   getBestFunVal();
	double*  getBestEver();

private:

	char exeDir_[FILENAME_MAX];
	std::string workdir_, cmaes_par_, cmaes_bounds_par_, priors_par_;

	cmaes_t evo_;
	int restart_;
	int lambda_;
	int step_;
	
	double gt0_, gt1_, gt2_, gt3_;
	double stt_;
	
	int dim_;
	double *lower_bound_, *upper_bound_;

	Density *priors_;

	double *const*pop_;
    double *arFunvals_; 

	static double (*fitfun_) (double*, int, void*, int*);
	static void taskfun_(double *x, int *no, double* res, int *info);
    double evaluate_population( cmaes_t *evo, double *arFunvals, double * const* pop, Density *d, int step );
};

#endif //ENGINE_CMAES_HPP


