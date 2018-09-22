#include "stats.c"
#include "log_tpdf.c"

//function [ LL, stats_loc ] = estimate_loglike_p( z, par, stats, data, model, MCMC )
//double LL = estimate_loglike_p( z, &par, &stats, &data, &model, &MCMC);


#define COPY_STATS
double estimate_loglike_p(double **z, par_t *par, stats_t *stats, data_t *data, model_t *model, MCMC_t *MCMC)
{
  int nu = 5; //% degrees of freedom for t-student distribution

  int Np = par->N;  // = size(z,1);
  int Ni = par->Ni; // = size(z,2);
  int M_ll = par->M_ll;


  //printf("Np = %d, Ni = %d\n", Np, Ni);
  par->L = 1;

  double *LL = (double *)calloc(1, M_ll*Ni*sizeof(double));

  //[ ~, ~, stats_loc ] = draw_z_samples_p( z, par, stats, data, model, MCMC );
  //  z, par, stats_loc

#ifdef COPY_STATS
  stats_t stats_loc;
  memset(&stats_loc, 0, sizeof(stats_loc));
  // copy stats to stats_loc
  copy_stats(par, MCMC, stats, &stats_loc);

  draw_z_samples_p( z, par, &stats_loc, data, model, MCMC );
#else
  draw_z_samples_p( z, par, stats, data, model, MCMC );
#endif

  //m_cond = stats_loc.m_cond;
  //v_cond = stats_loc.v_cond;
  double m_cond[Ni*Np];
  double v_cond[Ni*Np];

#ifdef COPY_STATS
  for (int i = 0; i < Ni*Np; i++) m_cond[i] = stats_loc.m_cond[i];
  for (int i = 0; i < Ni*Np; i++) v_cond[i] = stats_loc.v_cond[i];
#else
  for (int i = 0; i < Ni*Np; i++) m_cond[i] = stats->m_cond[i];
  for (int i = 0; i < Ni*Np; i++) v_cond[i] = stats->v_cond[i];
#endif

  double s_cond[Ni*Np];
  for (int i = 0; i < Ni*Np; i++) s_cond[i] = sqrt(v_cond[i]);

  double *cov_annealed_chol = par->omega_chol;

  for (int i = 0; i < Ni; i++)
  {
    double m_cond_loc[Np];
    double s_cond_loc[Np];

    for (int c = 0; c < Np; c++)
      m_cond_loc[c] = m_cond[Np*i+c];

    for (int c = 0; c < Np; c++)
      s_cond_loc[c] = s_cond[Np*i+c];

    for (int k = 0; k < par->M_ll; k++)
    {
      double r[Np];
        //r = trnd(nu,Np,1);

      for (int c = 0; c < Np; c++) r[c] = trnd(nu);

      double zp[Np];

      //z = m_cond_loc + s_cond_loc.*r; % Np x 1
      for (int c = 0; c < Np; c++)
        zp[c] = m_cond_loc[c] + s_cond_loc[c]*r[c]; //% Np x 1


      //[ e1, ~ ] = evaluate_all_models_p( i, z, par, data, model );
      double e1, sserr_tmp;
      evaluate_all_models_p(i, zp, par, data, model, &e1, &sserr_tmp);

      //void evaluate_all_models_p(int idx, double *z, par_t *par, data_t *data, model_t *model, double *ll, double *sserr)

      //dif = bsxfun( @minus, z, par.beta );
      double dif[Np];
      for (int c = 0; c < Np; c++)
        dif[c] = zp[c] - par->beta[c];

      double e2 = logp_eta_eval( dif, cov_annealed_chol, Np );

      //e3 = sum( log_tpdf(z,nu,m_cond_loc,s_cond_loc), 1 );  % 1 x Ni
      double e3tmp[Np];
      log_tpdf(zp, Np, nu, m_cond_loc, s_cond_loc, e3tmp);

      double e3 = 0;
      for (int c = 0; c < Np; c++) e3 += e3tmp[c];

      LL[k*Ni+i] = e1 + e2 - e3;

    }

  }


#if DEBUG
  print_matrix_2d_linear("LL", LL, M_ll, Ni);
#endif
  //LL = bsxfun( @rdivide, cumsum(exp(LL),1), (1:size(LL,1))' ); //'
  //LL = sum(log(LL),2);
  //stats_loc.LL = LL;

  for (int k = 0; k < M_ll; k++)
    for (int i = 0; i < Ni; i++)
      LL[k*Ni+i] = exp(LL[k*Ni+i]);

  for (int k = 1; k < M_ll; k++)  // peh: check memory layout of LL
    for (int i = 0; i < Ni; i++)
      LL[k*Ni+i] += LL[(k-1)*Ni+i];

#if DEBUG
  print_matrix_2d_linear("LL_cs", LL, M_ll, Ni);
#endif

  for (int k = 0; k < M_ll; k++)
    for (int i = 0; i < Ni; i++)
      LL[k*Ni+i] /= (k+1);

#if DEBUG
  print_matrix_2d_linear("LL_cs_dv", LL, M_ll, Ni);
#endif

  for (int k = 0; k < M_ll; k++)
    for (int i = 0; i < Ni; i++)
      LL[k*Ni+i] = log(LL[k*Ni+i]);

#if DEBUG
  print_matrix_2d_linear("log(LL_cs_dv)", LL, M_ll, Ni);
#endif

  double LLout[M_ll];
  memset(LLout, 0, M_ll*sizeof(double));

  for (int k = 0; k < M_ll; k++)
  {
    LLout[k] = 0;
    for (int i = 0; i < Ni; i++)
      LLout[k] += (LL[k*Ni+i]);

#if DEBUG
    printf("LLout[%d]= %lf\n", k, LLout[k]);
#endif
  }

  //LLv = LL(end);
  double LLv = LLout[M_ll-1]; // = LL(end);

#ifdef COPY_STATS
  // copy stats_loc to stats
  copy_stats(par, MCMC, &stats_loc, stats);
  free_stats(&stats_loc);
#endif

  free(LL);
  return LLv;
}
