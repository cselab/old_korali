#include <stdio.h>
#include "global.h"
#include "auxil.c"

#include "models/normal/my_model.c"

#include "functions/make_transf_z.c"
#include "functions/chol_eval.c"
#include "load_logistic.c"
#include "functions/evaluate_all_models_p.c"
#include "functions/init_saem.c"
#include "functions/logp_eta_eval.c"
#include "functions/draw_z_samples_p.c"
#include "functions/update_parameters.c"
#include "functions/estimate_loglike_p.c"
#include "functions/estimate_FIM.c"

//% estimate likelihood and FIM in one script

//addpath('./functions');
//clear; clc


model_t model;
data_t data;
par_t par;

stats_t stats;
double **z;

MCMC_t MCMC;
accum_t accum;


int main(int argc, char *argv[])
{
	memset(&model, 0, sizeof(model));
	memset(&data, 0, sizeof(data));
	memset(&par, 0, sizeof(par));
	memset(&stats, 0, sizeof(stats));
	memset(&MCMC, 0, sizeof(MCMC));
	memset(&accum, 0, sizeof(accum));

	gsl_rand_init(0);

	//%% load model + data + parameters

	//[model, data, par] = load_logistic();
	load_logistic(&model,&data,&par);

//	par.beta   = [ 201.4979   42.3234    1.0007 ]';
	par.beta[0] = 201.4979;
	par.beta[1] = 42.3234;
	par.beta[2] = 1.007;
	par.alpha   = 4.6279;

	par.L = 1;
	//	par.omega  = [ 17.8767    9.2859     0.1381 ].^2;
	double par_omega_data[3] = {17.8767, 9.2859, 0.1381};
	for (int i = 0; i < 3*3; i++) par.omega[i] = 0.0; // memset
	//par.omega  = diag( par.omega );
	for (int i = 0; i < 3; i++) par.omega[i*3+i] = pow(par_omega_data[i],2);
//	par.omega_chol = chol_eval( par.omega );
	chol_eval(par.omega, par.omega_chol, par.Nmp);

	//%% initialize variables

	//[ z, stats ] = init_saem( par, data, model );
	init_saem(&par,&data,&model,&z,&stats);

	print_stats(&stats,&par,&MCMC);

	//return 0;
	stats.save_samples = 1;

	//MCMC.steps = [1 1 1];	% [2 2 2]
	//MCMC.outer_steps = 10; % [1000]
	MCMC.steps[0] = 2; MCMC.steps[1] = 2; MCMC.steps[2] = 2; // MCMC.steps = [ 2 2 2 ];
	MCMC.outer_steps = 100;
	par.M_ll = 2e1; // 1e4


	//%%Estimate log-likelihood and save the z-samples to estimate the FIM

	//if (par.vectorized)
	//{
	//	[ LL, stats ] = estimate_loglike( z, par, stats, data, model,  MCMC );
	//}
	//else
	//{
	//	[ LL, stats ] = estimate_loglike_p( z, par, stats, data, model,  MCMC )

	double LL = estimate_loglike_p( z, &par, &stats, &data, &model, &MCMC);

	//}

	printf("\n -2*loglikelihood = %e \n\n", -2*LL );

	print_stats(&stats,&par,&MCMC);

	return 0;

	//%% Estimate and print FIM
	printf("\n\n");

	//[ FIM, stats ] = estimate_FIM( z, par, stats, model );

#if DEBUG
	{
		int Np = par.N;
		int Ni = par.Ni;
		int L = par.L;
		int Ns = (MCMC.steps[0]*1 + MCMC.steps[1]*Np + MCMC.steps[2]*Np) * MCMC.outer_steps;
		for (int i=0; i<Ns; i++)
		{
			int l = 0;
			for (int j = 0; j < Ni; j++)
			{
				for (int k = 0; k < Np; k++)
				  printf("%f ", stats.samples[(i)*L*Ni*Np+l*Ni*Np+j*Np+k]);
				printf("\n");
			}
  	}
	}
#endif


	printf("before estimate_FIM()\n");
	int Ns, M;
	double **FIM = estimate_FIM(NULL, &par, &stats, &model, &MCMC, &Ns, &M);
	printf("after  estimate_FIM()\n");

	//F = squeeze( FIM(end,:,:));
	double F[M][M];
	int end = Ns-1;

	printf("Ns = %d, M = %d\n", Ns, M);
	print_matrix_2d_linear("FIM[0]", (double *)FIM[0], M, M);
	print_matrix_2d_linear("FIM[end]", (double *)FIM[end], M, M);

	for (int j = 0; j < M; j++)
	for (int k = 0; k < M; k++)
		F[j][k] = FIM[end][j*M+k];


#if DEBUG
	double FIM_data[] = {
			0.031036,  -0.00014113,     0.084524,   -0.0073608,    0.0021019,     0.069913,      0.22411,
   -0.00014113,      0.10884,     0.069067,   0.00044366,     -0.06969,     -0.19031,      0.11918,
      0.084524,     0.069067,       484.91,     0.067533,     -0.70762,       477.39,      -90.555,
    -0.0073608,   0.00044366,     0.067533,     0.067597,    0.0023276,     0.077757,      0.18188,
     0.0021019,     -0.06969,     -0.70762,    0.0023276,      0.12415,     -0.89276,      -1.9108,
      0.069913,     -0.19031,       477.39,     0.077757,     -0.89276,       2536.1,       -72.42,
       0.22411,      0.11918,      -90.555,      0.18188,      -1.9108,       -72.42,      -163.59,
	};

	for (int j = 0; j < M; j++)
	for (int k = 0; k < M; k++)
		F[j][k] = FIM_data[j*M+k];
#endif

	//F = inv(F);
	double invF[M][M];

  print_matrix_2d_linear("F", (double *)F, M, M);

	inverse((double *)F, (double *)invF, M);

	memcpy((double *)F, (double *)invF, M*M*sizeof(double));

	print_matrix_2d_linear("invF", (double *)F, M, M);

	//F1 = F( 1:par.N, 1:par.N );
	//F2 = F( par.N+1:2*par.N, par.N+1:2*par.N );
	double F1[par.N][par.N];
	double F2[par.N][par.N];

	for (int j = 0; j < par.N; j++)
	for (int k = 0; k < par.N; k++)
		F1[j][k] = F[j][k];

	for (int j = 0; j < par.N; j++)
	for (int k = 0; k < par.N; k++)
		F2[j][k] = F[j+par.N][k+par.N];

	print_matrix_2d_linear("F1", (double *)F1, par.N, par.N);
	print_matrix_2d_linear("F2", (double *)F2, par.N, par.N);

	// peh: todo
	//T =  make_transf_z( par.transf , 1 );
	//J = diag( T(par.beta) );

	double J[par.N];

//	for (int j = 0; j < par.N; j++)
//		J[j] = 1.0; //par.beta[i];

//	transf_z(par.beta, par.N, 1, J);
	transf_z_ex(par.beta, par.transf, par.N, 1, J);


	print_matrix("J", J, par.N);
//	exit(1);

	//se1 = sqrt( diag( J'*F1*J ) );    %'
  double se1[par.N];

	double se1tmp[par.N][par.N];
	for (int j = 0; j < par.N; j++)
	for (int k = 0; k < par.N; k++)
		se1tmp[j][k] = F1[j][k];

	for (int j = 0; j < par.N; j++)
	for (int k = 0; k < par.N; k++)
		se1tmp[j][k] = J[j] * se1tmp[j][k];

	for (int j = 0; j < par.N; j++)
	for (int k = 0; k < par.N; k++)
		se1tmp[j][k] = se1tmp[j][k] * J[k];

	print_matrix_2d_linear("se1tmp", (double *)se1tmp, par.N, par.N);

	for (int j = 0; j < par.N; j++)
		se1[j] = sqrt(se1tmp[j][j]);

	//se2 = sqrt(diag(F2));
	double se2[par.N];
	for (int i = 0; i < par.N; i++)
		se2[i] = sqrt(F2[i][i]);

	//beta_pop = model.par_transf( par.beta );
	double beta_pop[par.N];
	for (int i = 0; i < par.N; i++)
		beta_pop[i] = par.beta[i];


	printf("======== beta_pop ============\n");
	for (int i=0; i <par.N; i++)
	{
		printf(" %f    %f    %2.1f \n",beta_pop[i], se1[i], 100*se1[i]/beta_pop[i] );
	}

	printf("======== omega    ============\n");
	for (int i=0; i <par.N; i++)
	{
		printf(" %f    %f    %2.1f \n",sqrt(par.omega[i*par.N+i]), se2[i], 100*se2[i]/sqrt(par.omega[i*par.N+i]) );
	}


	if (strcmp(model.error, "common") == 0)
	{
		printf("======== error    ============\n");
		int end = 2*par.N;
		double se3;
		if (F[end][end] < 0)
			se3 = 0;
		else
			se3	= sqrt(F[end][end]);
		printf(" %f    %f    %2.1f \n", par.alpha, se3, 100*se3/par.alpha );
	}
	printf("\n\n");

	return 0;
}
