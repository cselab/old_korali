#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern void my_model(double *t, double *theta, int flag, double *x);

void load_logistic(model_t *model, data_t *data, par_t *par)
{
	// LOAD MODEL
	strcpy(model->path, "models/normal/");

	model->fun = my_model;


	// LOAD DATA
	FILE *fp;

	fp = fopen("data/normal/all_data.txt", "r");
	char line[256];
	fgets(line, 256, fp);	// skip header
	int n = 0;
	while (fgets(line, 256, fp)!= NULL){
		n++;
	}
	printf("nlines = %d\n", n);
	data->id = malloc(n*sizeof(int));
	data->x	 = malloc(n*sizeof(double));
	data->y  = malloc(n*sizeof(double));

	rewind(fp);
	fgets(line, 256, fp);	// skip header
	for (int i = 0; i < n; i++){
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


	strcpy(model->error, "common"); //% ind or common
	strcpy(model->error_model, "constant"); //% constant or proportional

	int *par_transf = malloc(3*sizeof(int)); //[3];
	int error_par_transf;
	int length_par_beta;
	int length_omega;

	double *omega;

	if (strcmp(model->error,"ind")==0){
		par->beta = malloc(2*sizeof(double));
		for (int i = 0; i < 2; i++) par->beta[i] = 2;	// par.beta   = [ 2; 2 ];
		length_par_beta = 2;

		omega = malloc(2*sizeof(double));
		for (int i = 0; i < 2; i++) omega[i] = 1*1;	//par.omega  = [ 1; 1; ].^2;
		length_omega = 2;

		par->Nmp    = 1;	//length( par.beta  ) - 1;

		for (int i = 0; i < 1; i++) par_transf[i] = 0; //par_transf = [0 0 0];
		error_par_transf = 1 ;
	}
	else{
		par->beta = malloc(1*sizeof(double));
		par->beta[0] = 2; ;	//par.beta   = 2;
		length_par_beta = 1;

		omega = malloc(1*sizeof(double));
		for (int i = 0; i < 1; i++) omega[i] = 1; // par.omega  = [ 1 ].^2;
		length_omega = 1;

		par->alpha   = 1;

		par->Nmp     = 1; //length( par.beta  ) ;

		for (int i = 0; i < 1; i++) par_transf[i] = 0; //par_transf = [0 0 0];
		error_par_transf = 0;
	}

	par->vectorized = 0; //false;

	par->transf 	  = par_transf;			// peh: copy of pointer to allocated memory
	par->error_transf = error_par_transf;		// peh: added this, to store it somewhere
	par->length_omega = length_omega;		// peh: added this

	par->omega = calloc(1, length_omega*length_omega*sizeof(double));
	for (int i = 0; i < length_omega; i++)
		par->omega[i*length_omega+i] = omega[i];

	par->N  = length_par_beta; 
	par->Nd = data->Nd;

	model->par_tranf = transf_z;
	model->err_par_tranf = error_transf_z;

	par->omega_chol = calloc(1, length_omega*length_omega*sizeof(double));
	for (int i = 0; i < length_omega; i++){
		par->omega_chol[i*length_omega+i] = sqrt(par->omega[i*length_omega+i]);
	}

	par->Ni = data->Ni;

	// number of simulated chains per individual
	par->L = 5;

	// number of samples for log-likelihood estimation
	par->M_ll = 5000;

	double sigma0 = 1;
	if (par->vectorized)
	{
		exit(1);
	}
	else
	{
		par->sigma = malloc( par->L * par->Ni * sizeof(double *) );
		for(int i = 0; i < par->L * par->Ni; i++){
			
			par->sigma[i] = calloc(1, par->N * par->N * sizeof(double) );
			
			for (int si = 0; si < par->N; si++) 
				par->sigma[i][si*par->N+si] = sigma0;
		}

		//equivalent, more detailed, initialization code
		//for (int l = 0; l < par->L; l++)
		//for (int i = 0; i < par->Ni; i++)
		//{
		//	for (int si = 0; si < par->N; si++)
		//	par->sigma[l*par->Ni+i][si*par->N+si] = sigma0;
		//}
	}

}
