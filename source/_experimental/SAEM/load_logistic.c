#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern void my_model(double *t, double *theta, int flag, double *x);

//function [model, data, par] = load_logistic()
void load_logistic(model_t *model, data_t *data, par_t *par)
{
	//%% LOAD MODEL
	strcpy(model->path, "models/logistic/");
	//model.path = '../models/logistic/';
	//addpath(model.path);

	//model.fun    = @(theta,t) my_model(t,theta,1);
#if 1
	model->fun = my_model;
#endif


	//%% LOAD DATA
	//filename = '../data/logistic/all_data.txt';
	//delimiterIn = '\t'; headerlinesIn = 1;
	//A = importdata(filename,delimiterIn,headerlinesIn);
	//data.id =  A.data(:,1);            % subject ids
	//data.x =  A.data(:,2);             % independent variable
	//data.y =  A.data(:,3);             % dependent variable
	//data.uid = unique(data.id);        % subject's unique ids
	//data.Ni = length(data.uid);        % number of individuals
	//data.Nd = length(data.y);          % total number of data
	//for k = 1:data.Ni                  % number of data for i-th individual
	//	data.Ndi(k) = sum( data.id==k );
	//end
	FILE *fp;

	fp = fopen("data/logistic/all_data.txt", "r");
	char line[256];
	fgets(line, 256, fp);	// skip header
	int n = 0;
	while (fgets(line, 256, fp)!= NULL)
	{
		n++;
	}
	printf("nlines = %d\n", n);
	data->id = malloc(n*sizeof(int));
	data->x= malloc(n*sizeof(double));
	data->y= malloc(n*sizeof(double));

	rewind(fp);
	fgets(line, 256, fp);	// skip header
	for (int i = 0; i < n; i++)
	{
		fgets(line, 256, fp);	// skip header
		double id;
		sscanf(line, "%lf %lf %lf", &id, &data->x[i], &data->y[i]);
		data->id[i] = (int)id;
	}
	fclose(fp);

	//print_matrix("data->x", data->x, n);


	data->uid= malloc(n*sizeof(double));	// max unique ids

	data->uid[0] = data->id[0];
	int u_k = 0;
	for (int i = 1; i < n; i++)
	{
		if (data->uid[u_k] != data->id[i])
		{
			u_k++;
			data->uid[u_k] = data->id[i];
		}
	}

	//print_matrix_i("data->uid", data->uid, u_k+1);

	data->Ni = u_k+1;
	data->Nd = n;

	data->Ndi = calloc(1,data->Ni*sizeof(int));
	for (int i = 0; i < data->Ni; i++)	//% number of data for i-th individual
	{
		//data.Ndi(k) = sum( data.id==k );
		for (int k = 0; k < data->Nd; k++)
			if (data->id[k] == i) data->Ndi[i]++;
	}

	//%%
	strcpy(model->error, "common"); //% ind or common
	strcpy(model->error_model, "constant"); //% constant or proportional

	int *par_transf = malloc(3*sizeof(int)); //[3];
	int error_par_transf;
	int length_par_beta;
	int length_omega;

	double *omega;

	if (strcmp(model->error,"ind")==0)
	{
		par->beta = malloc(4*sizeof(double));
		for (int i = 0; i < 4; i++) par->beta[i] = 1;	// par.beta   = [ 1; 1; 1; 1 ];
		length_par_beta = 4;

		omega = malloc(4*sizeof(double));
		for (int i = 0; i < 4; i++) omega[i] = 100*1;	//par.omega  = 100*[ 1; 1; 1; 1 ].^2;
		length_omega = 4;

		par->Nmp    = 3;	//length( par.beta  ) - 1;

		for (int i = 0; i < 3; i++) par_transf[i] = 0; //par_transf = [0 0 0];
		error_par_transf = 1 ;
	}
	else
	{
		par->beta = malloc(3*sizeof(double));
		par->beta[0] = 100; par->beta[1] = 1; par->beta[2] = 1;	//par.beta   = [ 100; 1; 1 ];
		length_par_beta = 3;

		omega = malloc(3*sizeof(double));
		for (int i = 0; i < 3; i++) omega[i] = 2*1; // par.omega  = 2*[ 1; 1; 1 ].^2;
		length_omega = 3;

		par->alpha   = 1;

		par->Nmp     = 3; //length( par.beta  ) ;

		for (int i = 0; i < 3; i++) par_transf[i] = 0; //par_transf = [0 0 0];
		error_par_transf = 0;
	}

	//%%
	par->vectorized = 0; //false;

	par->transf = par_transf;			// peh: copy of pointer to allocated memory
	par->error_transf = error_par_transf;		// peh: added this, to store it somewhere
	par->length_omega = length_omega;		// peh: added this

	//par.omega  = diag( par.omega );
	par->omega = calloc(1, length_omega*length_omega*sizeof(double));
	for (int i = 0; i < length_omega; i++)
		par->omega[i*length_omega+i] = omega[i];

	par->N = length_par_beta; // length( par.beta  );
	par->Nd = data->Nd;

#if 1
//	model->par_transf     = make_transf_z( par_transf );
//	model->err_par_transf = make_transf_z( error_par_transf );

	model->par_tranf = transf_z; 
	model->err_par_tranf = error_transf_z;
#endif

	//par->omega_chol = chol_eval( par.omega );
	par->omega_chol = calloc(1, length_omega*length_omega*sizeof(double));
	for (int i = 0; i < length_omega; i++)
	{
		par->omega_chol[i*length_omega+i] = sqrt(par->omega[i*length_omega+i]);
	}

	par->Ni = data->Ni;

	//% number of simulated chains per individual
	par->L = 5;

	//% number of samples for log-likelihood estimation
	par->M_ll = 5000;

	double sigma0 = 1;
	if (par->vectorized)
	{
		exit(1);
#if 0
		par->sigma = sigma0*eye(par->N,par->N);
#endif
	}
	else
	{
		//par->sigma = zeros(par->L,par->Ni,par->N,par->N);
		//for l = 1 : par.L
			//for i = 1 : par.Ni
				//par.sigma(l,i,:,:) = sigma0*eye(par.N,par.N);
			//end
		//end

		par->sigma = malloc(par->L*par->Ni*sizeof(double *));
		for (int i = 0; i < par->L*par->Ni; i++)
		{
			par->sigma[i] = calloc(1, par->N*par->N*sizeof(double));
			for (int si = 0; si < par->N; si++) par->sigma[i][si*par->N+si] = sigma0;
		}

		//equivalent, more detailed, initialization code
		//for (int l = 0; l < par->L; l++)
		//for (int i = 0; i < par->Ni; i++)
		//{
		//	for (int si = 0; si < par->N; si++)
		//	par->sigma[l*par->Ni+i][si*par->N+si] = sigma0;
		//}
	}

	free(omega);
}
