/*
typedef struct stats_s
{
        double *num;
        double *denom;
        double *ar;
        double *m_cond;
        double *v_cond;
        double *samples;
        int compute_mv;
        int save_samples;
        double *ll;
        double *sserr;
        double *sserr_tot;
        int N_data_tot;
} stats_t;
*/


//function [ z, stats ] = init_saem( par, data, model )

void copy_stats(par_t *par, MCMC_t *MCMC, stats_t *stats, stats_t *stats_new)
{
	int Np = par->N;
	int Ni = par->Ni;
	int L  = par->L;
	int Nsteps = (MCMC->steps[0]*1 + MCMC->steps[1]*Np + MCMC->steps[2]*Np) * MCMC->outer_steps;

#if DEBUG
	printf("copy_stats: Nsteps = %d, stats_Nsteps = %d, new_stats_Nsteps = %d\n", Nsteps, stats->Nsteps, stats_new->Nsteps);
#endif
	if (stats_new->Nsteps == 0) stats_new->Nsteps = Nsteps;
	int min_Nsteps = 0;
	if (stats_new->Nsteps < stats->Nsteps)
		min_Nsteps = stats_new->Nsteps;
	else
		min_Nsteps = stats->Nsteps;
#if DEBUG
	printf("min_Nsteps = %d\n", min_Nsteps);
#endif

	if (stats_new->num == NULL)
	stats_new->num =       (double *)calloc(1,        L*Ni*Np*sizeof(double));	// zeros(L,Ni,Np);
	if (stats_new->denom == NULL)
	stats_new->denom =     (double *)calloc(1,        L*Ni*Np*sizeof(double));	// zeros(L,Ni,Np);
	if (stats_new->ar == NULL)
	stats_new->ar =        (double *)calloc(1,        L*Ni*Np*sizeof(double));	// zeros(L,Ni,Np);

	if (stats_new->m_cond == NULL)
	stats_new->m_cond =    (double *)calloc(1,          Ni*Np*sizeof(double));
	if (stats_new->v_cond == NULL)
	stats_new->v_cond =    (double *)calloc(1,          Ni*Np*sizeof(double));
	if (stats_new->samples == NULL)
	stats_new->samples =   (double *)calloc(1, Nsteps*L*Ni*Np*sizeof(double));;

	if (stats_new->ll == NULL)
	stats_new->ll =        (double *)calloc(1,        L*Ni*sizeof(double));	//zeros(L,Ni);
	if (stats_new->sserr == NULL)
	stats_new->sserr =     (double *)calloc(1,        L*Ni*sizeof(double));	//zeros(L,Ni);
	if (stats_new->sserr_tot == NULL)
	stats_new->sserr_tot = (double *)calloc(1, Nsteps*L*Ni*sizeof(double));


	for (int i = 0; i < L*Ni*Np; i++) stats_new->num[i] = stats->num[i];
	for (int i = 0; i < L*Ni*Np; i++) stats_new->denom[i] = stats->denom[i];
	for (int i = 0; i < L*Ni*Np; i++) stats_new->ar[i] = stats->ar[i];


	if (stats->m_cond != NULL) for (int i = 0; i < Ni*Np; i++) stats_new->m_cond[i] = stats->m_cond[i];
	if (stats->v_cond != NULL) for (int i = 0; i < Ni*Np; i++) stats_new->v_cond[i] = stats->v_cond[i];
	if (stats->samples != NULL) for (int i = 0; i < min_Nsteps*L*Ni*Np; i++) stats_new->samples[i] = stats->samples[i];

	if (stats->ll != NULL) for (int i = 0; i < L*Ni; i++) stats_new->ll[i] = stats->ll[i];
	if (stats->sserr != NULL) for (int i = 0; i < L*Ni; i++) stats_new->sserr[i] = stats->sserr[i];
	if (stats->sserr_tot != NULL) for (int i = 0; i < min_Nsteps*L*Ni; i++) stats_new->sserr_tot[i] = stats->sserr_tot[i];

	stats_new->compute_mv = stats->compute_mv;
	stats_new->save_samples = stats->save_samples;
	stats_new->N_data_tot = stats->N_data_tot;
}


void copy_stats2(par_t *par, MCMC_t *MCMC, stats_t *stats, stats_t *stats_new)
{
}

void print_stats(stats_t *stats, par_t *par, MCMC_t *MCMC)
{
	int Ni = par->Ni;
	int L  = par->L;
#if 0
	int Np = par->N;
	int Nsteps = (MCMC->steps[0]*1 + MCMC->steps[1]*Np + MCMC->steps[2]*Np) * MCMC->outer_steps;

	print_matrix("stats->num", stats->new, L*Ni*Np);
	print_matrix("stats->denom", stats->denom, L*Ni*Np);
	print_matrix("stats->ar", stats->ar, L*Ni*Np);
#endif

#if 0
	print_matrix("stats->m_cond", stats->m_cond, Ni*Np);
	print_matrix("stats->v_cond", stats->v_cond, Ni*Np);
	print_matrix("stats->samples", stats->samples, Nsteps*L*Ni*Np);
#endif

	print_matrix("stats->ll", stats->ll, L*Ni);
	print_matrix("stats->sserr", stats->sserr, L*Ni);
#if 0
	print_matrix("stats->sserr_tot", stats->sserr_tot, Nsteps*L*Ni*Np);
#endif
	print_matrix_i("stats->compute_mv", &stats->compute_mv, 1);
	print_matrix_i("stats->save_samples", &stats->save_samples, 1);
	print_matrix_i("stats->N_data_tot", &stats->N_data_tot, 1);
}


void free_stats(stats_t *stats)
{
	if (stats->num) free(stats->num);
	if (stats->denom) free(stats->denom);
	if (stats->ar) free(stats->ar);
	if (stats->m_cond) free(stats->m_cond);
	if (stats->v_cond) free(stats->v_cond);
	if (stats->samples) free(stats->samples);
        //int compute_mv;
        //int save_samples;
	if (stats->ll) free(stats->ll);
	if (stats->sserr) free(stats->sserr);
	if (stats->sserr_tot) free(stats->sserr_tot);
	// int N_data_tot;
	// int Nsteps;
}
