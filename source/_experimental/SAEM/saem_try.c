#define LOGISTIC	// comment out to enable NORMAL
#define DEBUG 0

#include <stdio.h>
#include "global.h"
#include "auxil.c"
#include "functions/make_transf_z.c"
#include "functions/chol_eval.c"

#ifdef LOGISTIC
	#include "models/logistic/my_model.c"
	#include "load_logistic.c"
#else // NORMAL
	#include "models/normal/my_model.c"
	#include "load_normal.c"
#endif

#include "functions/evaluate_all_models_p.c"
#include "functions/init_saem.c"
#include "functions/logp_eta_eval.c"
#include "functions/draw_z_samples_p.c"
#include "functions/update_parameters.c"
#include "functions/estimate_loglike_p.c"


model_t model;
data_t 	data;
par_t 	par;

stats_t stats;

double **z;

MCMC_t	MCMC;
accum_t accum;

int main(int argc, char *argv[])
{
	gsl_rand_init(0);

	load_logistic( &model, &data, &par );

	int Ngen = 100;

	init_saem( &par, &data, &model, &z, &stats );

	//variables for visualization
	double *par_transf_vec            = calloc(1, sizeof(double) * Ngen * par.N );
	double *omega_par_transf_vec      = calloc(1, sizeof(double) * Ngen * par.N );
	double *err_par_transf_vec        = calloc(1, sizeof(double) * Ngen);
	double *omega_err_par_transf_vec  = calloc(1, sizeof(double) * Ngen);


	for( int i = 0; i<Ngen; i++ ){
		MCMC.steps[0] = 2; MCMC.steps[1] = 2; MCMC.steps[2] = 2;
		MCMC.outer_steps = 1;
		
		draw_z_samples_p( z, &par , &stats,  &data, &model, &MCMC );

		update_parameters(&par, i, z, &accum, &stats);

		//visualization
		printf("==================== Generation %d ====================\n", i);

		double sqrt_omega_par_transf_vec_i[par.N];
		double sqrt_omega_err_par_transf_vec_i;

		double omega_tmp[par.N];
		for( int j=0; j<par.N; j++ ) omega_tmp[j] = par.omega[j*par.N+j];

		int deriv = 0;

		model.par_tranf( par.beta, par.Nmp, deriv, &par_transf_vec[i*par.N] );

		for (int j=0; j<par.N; j++ ) omega_par_transf_vec[i*par.N+j] = omega_tmp[j];

		if (strcmp(model.error,"ind")==0){
			model.err_par_tranf( &par.beta[par.N], 1, deriv, &err_par_transf_vec[i]);
			omega_err_par_transf_vec[i] = omega_tmp[par.Nmp+1];
		}
		else{
			model.err_par_tranf(&par.alpha, 1, deriv, &err_par_transf_vec[i]);
		}

		// get the square roots
		for (int j = 0; j < par.N; j++)
			sqrt_omega_par_transf_vec_i[j] = sqrt(omega_par_transf_vec[i*par.N+j]);

		sqrt_omega_err_par_transf_vec_i = sqrt(omega_err_par_transf_vec[i]);


		print_matrix(" par_transf_vec_i", &par_transf_vec[i*par.N], par.N);
		print_matrix("sqrt(omega_par_transf_vec_i)", sqrt_omega_par_transf_vec_i, par.N);
		print_matrix("err_par_transf_vec_i", &err_par_transf_vec[i], 1);
		print_matrix("sqrt(omega_err_par_transf_vec_i)", &sqrt_omega_err_par_transf_vec_i, 1);

#if 1
		MCMC.outer_steps = 50;
		double LL = estimate_loglike_p( z, &par, &stats, &data, &model, &MCMC);
		printf("\n -2*loglikelihood = %e \n\n", -2*LL );
#endif

	}


	free(par_transf_vec);
	free(omega_par_transf_vec);
	free(err_par_transf_vec);
	free(omega_err_par_transf_vec);

	return 0;

}
