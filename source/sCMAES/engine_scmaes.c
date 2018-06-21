#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "surrogate.h"
#include "surrogate_samples.h"

#include <cmaes.h>
#include <cmaes_utils.h>
#include <fitfun.h>
#include <priors.h>

#if defined(_USE_TORC_)
#include <mpi.h>
#include <torc.h>
#endif


enum {
    NSTEPS_ARCHIVE = 3
};

#define VERBOSE 1
#define _STEALING_
#define _IODUMP_ 1
#define JOBMAXTIME    0
#define _RESTART_

void taskfun(double *x, int *pn, double *res, int *info);
double evaluate_population( cmaes_t *evo, double *arFunvals, double * const* pop, Density *d, int step );

void add_population_to_surrogate(int lambda, double *const *pop, double *arFunvals, Surrogate *s);
double evaluate_population_surrogate( cmaes_t *evo, double *arFunvals, double * const* pop, Density *d, int step, Surrogate *s );

int main(int argn, char **args) {
    cmaes_t evo; 
    double *arFunvals, *const*pop;
    int lambda, dim;
    double gt0, gt1, gt2, gt3;
    double stt = 0.0, dt;
    char dim_str[12];
    int step = 0;

    Surrogate *surrogate;
    Surrogate_pop *surrogate_pop;
    Archive *archive;
    
    static int checkpoint_restart = 0;

    double *lower_bound, *upper_bound;


#if defined(_USE_TORC_)
    torc_register_task(taskfun);
    torc_init(argn, args, MODE_MS);
#endif


    gt0 = get_time();

        
    if ( argn==2  &&  !strcmp(args[1], "-cr") )
	checkpoint_restart = 1;

    
    arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "cmaes_initials.par");
    printf("%s\n", cmaes_SayHello(&evo));
    cmaes_ReadSignals(&evo, "cmaes_signals.par");  


    dim    = cmaes_Get(&evo, "dim");
    lambda = cmaes_Get(&evo, "lambda");
    cmaes_utils_read_bounds(VERBOSE, "cmaes_bounds.par", &lower_bound, &upper_bound, dim );

    surrogate_ini(dim, &surrogate);
    surrogate_pop_ini(dim, lambda, &surrogate_pop);
    archive_ini(dim, NSTEPS_ARCHIVE * lambda, &archive);
    

    // Initialize prior distributions
    Density *priors;
    int Nprior;
    read_priors( "priors.par", &priors, &Nprior );

    if( Nprior != dim ){
        printf("The dimension of the prior is different from the dimension of the problem. Exit...");
        exit(1);
    }


    // Initialize log-likelihood
    sprintf(dim_str, "%d", dim );
    fitfun_initialize( dim_str );


    gt1 = get_time();
        
    while ( !cmaes_TestForTermination(&evo) ) {
        
        pop = cmaes_SamplePopulation(&evo); 

        if( checkpoint_restart ){
            dt = cmaes_utils_load_pop_from_file( VERBOSE, step, pop, arFunvals, dim, lambda, &checkpoint_restart);
    	} else {
            cmaes_utils_make_all_points_feasible( &evo, pop, lower_bound, upper_bound );
            dt = evaluate_population( &evo, arFunvals, pop, priors, step );
        }
        stt += dt;

        archive_add(archive, lambda, pop, arFunvals);
        archive_mark_candidates(archive, 8.0, &evo);
        surrogate_pop_select_from_archive(surrogate_pop, archive);

        surrogate_reset(surrogate);
        add_population_to_surrogate(
            surrogate_pop_get_n(surrogate_pop),
            surrogate_pop_get_pop(surrogate_pop),
            surrogate_pop_get_funvals(surrogate_pop),
            surrogate);
        surrogate_optimize(surrogate);
        
        cmaes_UpdateDistribution(1, &evo, arFunvals);

        cmaes_ReadSignals(&evo, "cmaes_signals.par"); fflush(stdout);

        if (VERBOSE) cmaes_utils_print_the_best( evo, step );
		
       	if (!checkpoint_restart){
            cmaes_utils_write_pop_to_file( evo, arFunvals, pop, step );
        }

#if defined(_RESTART_)
        cmaes_WriteToFile(&evo, "resume", "allresumes.dat");
#endif

        if( ! cmaes_utils_is_there_enough_time( JOBMAXTIME, gt0, dt ) ){
            evo.sp.stopMaxIter=step+1;
            break;
        }
        
        step++;
    }
    gt2 = get_time();

    printf("Stop:\n %s \n",  cmaes_TestForTermination(&evo)); /* print termination reason */
    cmaes_WriteToFile( &evo, "all", "allcmaes.dat" );         /* write final results */
    cmaes_exit(&evo); /* release memory */

    gt3 = get_time();
    
    printf("Total elapsed time      = %.3lf  seconds\n", gt3-gt0);
    printf("Initialization time     = %.3lf  seconds\n", gt1-gt0);
    printf("Processing time         = %.3lf  seconds\n", gt2-gt1);
    printf("Funtion Evaluation time = %.3lf  seconds\n", stt);
    printf("Finalization time       = %.3lf  seconds\n", gt3-gt2);

    surrogate_fin(surrogate);
    surrogate_pop_fin(surrogate_pop);
    archive_fin(archive);


#if defined(_USE_TORC_)
    torc_finalize();
#endif
    
		
    return 0;
}


/*
  Assumptions: the feasible domain is convex, the optimum is
  not on (or very close to) the domain boundary, initialX is
  feasible and initialStandardDeviations are sufficiently small
  to prevent quasi-infinite looping. 
*/

// the function to be minimized
void taskfun(double *x, int *n, double *res, int *info)
{
    (*res) = - fitfun(x, *n, (void *)NULL, info);    // minus for minimization
}

double evaluate_population( cmaes_t *evo, double *arFunvals, double * const* pop, Density *d, int step ) {
    int lambda = cmaes_Get( evo, "lambda");
    int dim    = cmaes_Get( evo, "dim");
    int info[4];
    double tt0, tt1 ;
    	
    tt0 = get_time();
	
    for( int i = 0; i < lambda; ++i){
        info[0] = 0; info[1] = 0; info[2] = step; info[3] = i;     /* gen, chain, step, task */
		
#if defined(_USE_TORC_)
        torc_create( -1, taskfun, 4,
                     dim, MPI_DOUBLE, CALL_BY_VAL,
                     1, MPI_INT, CALL_BY_COP,
                     1, MPI_DOUBLE, CALL_BY_RES,
                     4, MPI_INT, CALL_BY_COP,
                     pop[i], &dim, &arFunvals[i], info);
#else
        taskfun(pop[i], &dim, &arFunvals[i], info);
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
    for( int i=0; i<lambda; i++){
        arFunvals[i] -= prior_log_pdf(d, dim, pop[i]);
    }


    tt1 = get_time();
  
    return tt1-tt0;
}

double evaluate_population_surrogate( cmaes_t *evo, double *arFunvals, double * const* pop, Density *d, int step, Surrogate *s ){
    int i, dim, lambda;
    double tt0, tt1 ;

    lambda = cmaes_Get( evo, "lambda");
    dim    = cmaes_Get( evo, "dim");
    
    tt0 = get_time();
	
    for (i = 0; i < lambda; ++i)
        arFunvals[i] = surrogate_eval(pop[i], s);
    
    // subtract the log-prior from the log-likelohood
    for (i = 0; i < lambda; ++i)
        arFunvals[i] -= prior_log_pdf(d, dim, pop[i]);

    tt1 = get_time();
  
    return tt1-tt0;
}

void add_population_to_surrogate(int lambda, double *const *pop, double *arFunvals, Surrogate *s) {
    int i;
    double *x, y;
    for (i = 0; i < lambda; ++i) {
        x = pop[i];
        y = arFunvals[i];
        surrogate_add_point(x, y, s);
    }
}
