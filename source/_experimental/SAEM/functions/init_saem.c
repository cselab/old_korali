//function [ z, stats ] = init_saem( par, data, model )
void init_saem(par_t *par, data_t *data, model_t *model, double ***zz, stats_t *stats)
{
	int Np = par->N;
	int Ni = par->Ni;
	int L  = par->L;

	//%% initialize z
	double **z;
	z = malloc(par->L*sizeof(double *));
	for (int i = 0; i < par->L; i++)
		z[i] = calloc(1, par->N * data->Ni * sizeof(double));

	*zz = z; 	// write to the output

	//for i=1:par.L
	//	z(:,:,i) = mvnrnd( par.beta', par.omega, data.Ni )';
	//end

	for (int i = 0; i < L; i++)
	{
		//z(:,:,i) = mvnrnd( par.beta', par.omega, data.Ni )';
		//z[i][] = mvnrnd(...);
		double rnd[Np];
		for (int j = 0; j < Ni; j++)
		{
			mvnrnd(par->beta, par->omega, rnd, Np);
			for (int k = 0; k < Np; k++)
				z[i][k*Ni+j] = rnd[k];
		}
	}

	//%% initialize statss
	if (par->vectorized)
	{
		//stats.num   = zeros(1,Np);
		//stats.denom = zeros(1,Np);
		//stats.ar    = zeros(1,Np);
	}
	else
	{
		stats->num   = calloc(1, L*Ni*Np*sizeof(double));	// zeros(L,Ni,Np);
		stats->denom = calloc(1, L*Ni*Np*sizeof(double));	// zeros(L,Ni,Np);
		stats->ar    = calloc(1, L*Ni*Np*sizeof(double));	// zeros(L,Ni,Np);
	}

	stats->m_cond = 0;
	stats->v_cond = 0;

	stats->samples = NULL; //[];

	//% compute mean and variance
	stats->compute_mv = 1;

	//% save all samples
	stats->save_samples = 1;


	//% initialize likelihoods and standard errors
	//% this loop can be parallelized
	stats->ll    = calloc(1, L*Ni*sizeof(double));	//zeros(L,Ni);
	stats->sserr = calloc(1, L*Ni*sizeof(double));	//zeros(L,Ni);

	if (par->vectorized)
	{
		printf("par->vectorized = %d: this should not happen!", par->vectorized);
		abort();
//		for l = 1 : L
//			[ stats.ll(l,:), stats.sserr(l,:) ] = evaluate_all_models( z(:,:,l), par, data, model );
//		end
        }
	else
	{
		//for l = 1 : L
		//	for i=1:Ni
		//		[ stats.ll(l,i), stats.sserr(l,i) ] = evaluate_all_models_p( i, z(:,i,l), par, data, model);
		double zc[Np];

		for (int l = 0; l < L; l++)
		{
			for (int i = 0; i < Ni; i++)
			{
				for (int k = 0; k < Np; k++) zc[k] = z[l][k*Ni+i];

#if VERBOSE
				print_matrix("zc", zc, Np);
#endif
				//[ stats.ll(l,i), stats.sserr(l,i) ] = evaluate_all_models_p( i, z(:,i,l), par, data, model);
				//zc = z(:,i,l);
				evaluate_all_models_p(i, zc, par, data, model, &stats->ll[l*Ni+i], &stats->sserr[l*Ni+i]);
			}
		}
	}

	return;
}
