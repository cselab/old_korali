

//function [ z, par, stats ] = draw_z_samples_p( z, par, stats, data, model, MCMC )
void draw_z_samples_p( double **z, par_t *par, stats_t *stats, data_t *data, model_t *model, MCMC_t *MCMC)
{
	double a_star = 0.4;
	double delta  = 0.4;

	int Np = par->N;
	int Ni = data->Ni;
	int L  = par->L;

	if (Np==1)
	{
		MCMC->steps[2] = 0;	//% force zero steps in the 3rd MCMC
	}

	//index_diag = ( eye(Np,Np)~=0 );
  //3x3 logical array
  // 1   0   0
  // 0   1   0
  // 0   0   1
	int index_diag[Np*Np];
	memset(index_diag, 0, Np*Np*sizeof(int));
	for (int i = 0; i < Np; i++) index_diag[i*Np+i] = 1;

	int Nsteps = (MCMC->steps[0]*1 + MCMC->steps[1]*Np + MCMC->steps[2]*Np) * MCMC->outer_steps;	//sum( MCMC.steps.*[1,Np,Np] ) * MCMC.outer_steps;

	//printf("Nsteps = %d\n", Nsteps);

	//stats.samples   = zeros(Nsteps,Np,Ni,L);
	//stats.sserr_tot = zeros(Nsteps,Ni,L);
	if (Nsteps > stats->Nsteps) {
		free(stats->samples);
		free(stats->sserr_tot);
		stats->samples = NULL;
		stats->sserr_tot = NULL;
		stats->Nsteps = Nsteps;
	}

	if (stats->samples == NULL)
		stats->samples = calloc(1, Nsteps*L*Ni*Np*sizeof(double));	// Nsteps*Np*Ni*L
	else
		memset(stats->samples, 0, Nsteps*L*Ni*Np*sizeof(double));

	if (stats->sserr_tot == NULL)
		stats->sserr_tot = calloc(1, Nsteps*L*Ni*sizeof(double));
	else
		memset(stats->sserr_tot, 0, Nsteps*L*Ni*sizeof(double));


	//zero_mean = zeros(Np,1);
	double zero_mean[Np];
	memset(zero_mean, 0, Np*sizeof(double));


	double **eta;
	eta = malloc(par->L*sizeof(double *));
	for (int i = 0; i < par->L; i++)
		eta[i] = calloc(1, par->N * data->Ni * sizeof(double));


	//eta = bsxfun( @minus, z, par.beta);
	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < Ni; j++)
		{
			for (int k = 0; k < Np; k++)
				eta[i][k*Ni+j] = z[i][k*Ni+j] - par->beta[k];
		}
	}


	//if( isempty( stats.ll ) )
	//else
	//	ll    = stats.ll;
	//	sserr = stats.sserr;
	//end
	double *ll    = calloc(1, L*Ni*sizeof(double));
	double *sserr = calloc(1, L*Ni*sizeof(double));

	if (stats->ll != NULL)
		for (int i = 0; i < L*Ni; i++) ll[i] = stats->ll[i];

	if (stats->sserr != NULL)
		for (int i = 0; i < L*Ni; i++) sserr[i] = stats->sserr[i];


	double *cov_annealed = par->omega;		// peh: this should be ok
	double *cov_annealed_chol = par->omega_chol;

	int cnt_s = 0;	// moved here

	//% main loop
	for (int l = 0; l < L; l++)	//% number of chains for each individual (l=1:L)
	{
		for (int i = 0; i < Ni; i++)	//% loop over the individuals
		{
			double eta_loc[Np];
			for (int c = 0; c < Np; c++) eta_loc[c] = eta[l][c*Ni+i];				//eta_loc   = squeeze( eta(:,i,l) );

			double ll_loc    = ll[l*Ni+i];
			double sserr_loc = sserr[l*Ni+i];

			double sigma_loc[Np*Np];
			double num_loc[Np];
			double denom_loc[Np];
			double ar_loc[Np];

			//sigma_loc = squeeze( par->sigma(l,i,:,:) );
			for (int ci = 0; ci < Np; ci++)
			for (int cj = 0; cj < Np; cj++)
				sigma_loc[ci*Np+cj] = par->sigma[l*Ni+i][ci*Np+cj];

#if VERBOSE
			print_matrix_2d_linear("sigma_loc loaded", sigma_loc, Np, Np);
#endif

			for (int ci = 0; ci < Np; ci++)
				num_loc[ci] = stats->num[l*Ni*Np+i*Np+ci];		//num_loc   = squeeze( stats->num(l,i,:) );

			for (int ci = 0; ci < Np; ci++)
				denom_loc[ci] = stats->denom[l*Ni*Np+i*Np+ci];		//denom_loc = squeeze( stats->denom(l,i,:) );

			for (int ci = 0; ci < Np; ci++)
			{
				ar_loc[ci]    = num_loc[ci] / denom_loc[ci];
				if (num_loc[ci] == 0 && denom_loc[ci] == 0) ar_loc[ci] = 0.0;
			}

			//int cnt_s = 0;	// was 1, moved outside the loops
			cnt_s = 0;

			double z_new[Np];
			double z_loc[Np];

			for (int k = 0; k <MCMC->outer_steps; k++)
			{
				for (int m = 0; m < MCMC->steps[0]; m++)
				{
					double eta_new[Np];
					mvnrnd( zero_mean, cov_annealed, eta_new, Np);

#if VERBOSE
					print_matrix("eta_new", eta_new, Np);
#endif
					//xxx -> double z_new[Np];
					//z_new = bsxfun( @plus, par.beta, eta_new);
					for (int c = 0; c < Np; c++)
						z_new[c] = par->beta[c] + eta_new[c];

					double ll_new, sserr_new;
					evaluate_all_models_p( i, z_new, par, data, model, &ll_new, &sserr_new);

					//printf("ll_new = %lf\n", ll_new);
					//printf("sserr_new = %lf\n", sserr_new);

					double dif = ll_new - ll_loc;

					int accept = log(uniformrand(0,1)) < dif;	// peh: for testing xxx
					if (accept)
					{
						for (int c = 0; c < Np; c++) eta_loc[c] = eta_new[c];	// eta_loc = eta_new;
						ll_loc    =  ll_new;
						sserr_loc =  sserr_new;
					}

					//xxx->double z_loc[Np];
					for (int c = 0; c < Np; c++) z_loc[c] = par->beta[c] + eta_loc[c];

					//% always save samples
					//stats->samples(cnt_s,:,i,l) = z_loc;
					for (int c = 0; c < Np; c++)
						stats->samples[cnt_s*L*Ni*Np+l*Ni*Np+i*Np+c] = z_loc[c];	// Nsteps*L*Ni*Np

					stats->sserr_tot[cnt_s*L*Ni+l*Ni+i] = sserr_loc;
					cnt_s = cnt_s + 1;
				} // % end of MCMC(1) loop

				double eta_new[Np];
				for (int c = 0; c < Np; c++) eta_new[c] = eta_loc[c];	// eta_new = eta_loc;

				double logp_eta = logp_eta_eval( eta_loc, cov_annealed_chol, Np);
#if DEBUG
				printf("logp_eta = %lf\n", logp_eta);
#endif
				for (int m = 0; m < MCMC->steps[1]; m++)
				{
					for (int j=0; j < Np; j++)
					{
						eta_new[j] = normalrand(eta_loc[j], sqrt(sigma_loc[j*Np+j]));
						z_new[j]   = par->beta[j] + eta_new[j];

						double ll_new, sserr_new;
						evaluate_all_models_p( i, z_new, par, data, model, &ll_new, &sserr_new);

						double logp_eta_new = logp_eta_eval( eta_new, cov_annealed_chol, Np);

						double dif = ll_new - ll_loc + logp_eta_new - logp_eta;

						int accept = log(uniformrand(0,1)) < dif;	// peh: xxx testing
						if (accept)
						{
							for (int c = 0; c < Np; c++) eta_loc[c] = eta_new[c];	// eta_loc = eta_new;
							ll_loc    =  ll_new;
							sserr_loc =  sserr_new;
							logp_eta  =  logp_eta_new;
						}

						for (int c = 0; c < Np; c++) z_loc[c] = par->beta[c] + eta_loc[c];
						for (int c = 0; c < Np; c++) z_new[c] = z_loc[c];
						for (int c = 0; c < Np; c++) eta_new[c] = eta_loc[c];

						num_loc[j]   = num_loc[j] + accept;
						denom_loc[j] = denom_loc[j] + 1;
						ar_loc[j]    = num_loc[j] / denom_loc[j];

						//% always save samples
						//stats->samples(cnt_s,:,i,l) = z_loc;
						//stats->sserr_tot(cnt_s,i,l) = sserr_loc;
						for (int c = 0; c < Np; c++)
							stats->samples[cnt_s*L*Ni*Np+l*Ni*Np+i*Np+c] = z_loc[c];	// Nsteps*L*Ni*Np
						stats->sserr_tot[cnt_s*L*Ni+l*Ni+i] = sserr_loc;
						cnt_s = cnt_s + 1;
					}

				}

				//...sigma_loc(index_diag) = sigma_loc(index_diag) .* ( 1 + delta*( ar_loc-a_star ) );
				for (int c = 0; c < Np; c++)
					sigma_loc[c*Np+c] = sigma_loc[c*Np+c] * (1 + delta*(ar_loc[c]-a_star));

#if VERBOSE
				print_matrix_2d_linear("sigma_loc", sigma_loc, Np, Np);
#endif

				for (int c = 0; c < Np; c++) eta_new[c] = eta_loc[c];	// eta_new = eta_loc;

				logp_eta = logp_eta_eval( eta_loc, cov_annealed_chol, Np);

				for (int m = 0; m < MCMC->steps[2]; m++)
				{
					for (int j=0; j < Np; j++)
					{ // peh: check this
						//Neta = randi(Np-1);
						//rset = sort( randsample(2:Np,Neta) );

						//long Neta = 2;
						long Neta = randi(Np);						// 1 + gsl_rng_uniform_int(r[0], Np-1);
#if DEBUG
						printf("Neta = %ld\n", Neta);
#endif
						double rset[Neta];
						randsample_sorted(2, Np, rset, Neta);

						//double rset[Neta], b[Np-1];
						//for (int i = 0; i < Np-1; i++)
						//{
						//	b[i] = 2 + (double) i;
						//}
						//gsl_ran_choose (r[0], rset, Neta, b, Np-1, sizeof (double));

						//...eta_new(rset) = mvnrnd( eta_loc(rset), sigma_loc(rset,rset) );
						double eta_new_rset[Neta];
						double eta_loc_rset[Neta];
						double sigma_loc_rset[Neta*Neta];

#if VERBOSE
						print_matrix("rset", rset, Neta);
#endif

						for (int r = 0; r < Neta; r++)
						{
							long idx = (long)rset[r]-1;
							eta_loc_rset[r] = eta_loc[idx];
						}

						for (int ri = 0; ri < Neta; ri++)
						for (int rj = 0; rj < Neta; rj++)
						{
							int ri_idx = (long)rset[ri]-1;
							int rj_idx = (long)rset[rj]-1;
							sigma_loc_rset[ri*Neta+rj] = sigma_loc[ri_idx*Np+rj_idx];
						}

#if VERBOSE
						print_matrix_2d_linear("sigma_loc_rset", sigma_loc_rset, Neta, Neta);
#endif

						mvnrnd(eta_loc_rset, sigma_loc_rset, eta_new_rset, Neta);

						for (int r = 0; r < Neta; r++)
						{
							long idx = (long)rset[r]-1;
							eta_new[idx] = eta_new_rset[r];
						}

						//z_new(rset,:) = bsxfun( @plus, par.beta(rset), eta_new(rset,:) );
						for (int r = 0; r < Neta; r++)
						{
							long idx = rset[r]-1;
							z_new[idx] = par->beta[idx] + eta_new[idx];
						}

						double ll_new, sserr_new;
						evaluate_all_models_p( i, z_new, par, data, model, &ll_new, &sserr_new);

						double logp_eta_new = logp_eta_eval( eta_new, cov_annealed_chol, Np );

						double dif = ll_new - ll_loc + logp_eta_new - logp_eta;

						int accept =  log(uniformrand(0,1)) < dif ;	// peh: xxx testing
						if (accept)
						{
							for (int c = 0; c < Np; c++) eta_loc[c] = eta_new[c];	// eta_loc = eta_new;
							ll_loc    =  ll_new;
							sserr_loc =  sserr_new;
							logp_eta  =  logp_eta_new;
						}

						//z_loc = bsxfun( @plus, par.beta, eta_loc);
						for (int c = 0; c < Np; c++)
							z_loc[c] = par->beta[c] + eta_loc[c];

						for (int c = 0; c < Np; c++) z_new[c] = z_loc[c];	//z_new = z_loc;
						for (int c = 0; c < Np; c++) eta_new[c] = eta_loc[c];	//eta_new = eta_loc;

						for (int r = 0; r < Neta; r++)
						{
							long idx = rset[r];

							//num_loc(rset)   = num_loc(rset) + accept;
							//denom_loc(rset) = denom_loc(rset) + 1;
							//ar_loc(rset)    = num_loc(rset) / denom_loc(rset);
							num_loc[idx] = num_loc[idx] + accept;
							denom_loc[idx] = denom_loc[idx] + 1;
							ar_loc[idx] = num_loc[idx] / denom_loc[idx];
						}

						//% always save samples
						//stats->samples(cnt_s,:,i,l) = z_loc;
						//stats->sserr_tot(cnt_s,i,l) = sserr_loc;
						for (int c = 0; c < Np; c++)
							stats->samples[cnt_s*L*Ni*Np+l*Ni*Np+i*Np+c] = z_loc[c];	// Nsteps*L*Ni*Np
						stats->sserr_tot[cnt_s*L*Ni+l*Ni+i] = sserr_loc;

						cnt_s = cnt_s + 1;
					}
				}

				if (MCMC->steps[2] > 0)
				{
					//sigma_loc(index_diag) = sigma_loc(index_diag) .* ( 1 + delta*( ar_loc-a_star ) );
					for (int c = 0; c < Np; c++)
						sigma_loc[c*Np+c] = sigma_loc[c*Np+c] * (1 + delta*(ar_loc[c]-a_star));
				}
			}	// % end of outer MCMC loop

			//eta(:,i,l) = eta_loc;
			for (int c = 0; c < Np; c++)
				eta[l][c*Ni+i] = eta_loc[c];

			ll[l*Ni+i]    = ll_loc;
			sserr[l*Ni+i] = sserr_loc;

			for (int ci = 0; ci < Np; ci++)
			for (int cj = 0; cj < Np; cj++)
				par->sigma[l*Ni+i][ci*Np+cj] = sigma_loc[ci*Np+cj];

			for (int ci = 0; ci < Np; ci++)
				stats->num[l*Ni*Np+i*Np+ci] = num_loc[ci];

			for (int ci = 0; ci < Np; ci++)
				stats->denom[l*Ni*Np+i*Np+ci] = denom_loc[ci];

			for (int ci = 0; ci < Np; ci++)
				stats->ar[l*Ni*Np+i*Np+ci] = num_loc[ci] / denom_loc[ci];
		} //% end of Ni individuals loop
	} //% end of L chains loop

	//% estimate conditional mean and variances
	//% it can simplified in a for loop
	if (stats->m_cond == NULL)
		stats->m_cond = calloc(1, Ni*Np*sizeof(double));
	else
		memset(stats->m_cond, 0, Ni*Np*sizeof(double));

	if (stats->v_cond == NULL)
		stats->v_cond = calloc(1, Ni*Np*sizeof(double));
	else
		memset(stats->v_cond, 0, Ni*Np*sizeof(double));

#if DEBUG
	printf("cnt_s = %d\n", cnt_s);
#endif
	for (int i = 0; i < Ni; i++)
	{
		//tmp = squeeze( stats->samples(:,:,i,:) );
		//tmp = reshape( shiftdim(tmp,2), [], Np);
		//stats->m_cond(:,i) = mean(tmp,2);
		//stats->v_cond(:,i) = var(tmp,[],2);

		for (int c = 0; c < Np; c++)
		{
			double sum;
			sum = 0.0;
			for (int cnt = 0; cnt < cnt_s; cnt++)
			{
				for (int l = 0; l < L; l++)
				{
						double val = stats->samples[cnt*L*Ni*Np+l*Ni*Np+i*Np+c]; 	// Nsteps*L*Ni*Np
						sum += val;
				} // l
			} // cnt
			double mean_v = sum/(cnt_s*L);

			double sum2;
			sum2 = 0.0;
			for (int cnt = 0; cnt < cnt_s; cnt++)
			{
				for (int l = 0; l < L; l++)
				{
					double val = stats->samples[cnt*L*Ni*Np+l*Ni*Np+i*Np+c]; 	// Nsteps*L*Ni*Np
					sum2 += (val-mean_v)*(val-mean_v);
				} // l
			} // cnt
			double var_v = sum2/(cnt_s*L-1);

			stats->m_cond[i*Np+c] = mean_v;
			stats->v_cond[i*Np+c] = var_v;

		} // c	(parameter)
	} // i (individual)

	//%	stats.m_cond
	//%	stats.v_cond
	print_matrix_2d_linear("stats->m_cond", stats->m_cond, Ni, Np);
	print_matrix_2d_linear("stats->v_cond", stats->v_cond, Ni, Np);

	// z = bsxfun( @plus, par.beta, eta);
	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < Ni; j++)
		{
			for (int k = 0; k < Np; k++)
				z[i][k*Ni+j] = eta[i][k*Ni+j] + par->beta[k];
		}
	}


	for (int i = 0; i < L*Ni; i++) stats->ll[i] = ll[i];
	for (int i = 0; i < L*Ni; i++) stats->sserr[i] = sserr[i];

	stats->N_data_tot = data->Nd;    //% total number of data

	free(ll);
	free(sserr);
	for (int i = 0; i < par->L; i++)
		free(eta[i]);
	free(eta);
}
