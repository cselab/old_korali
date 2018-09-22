//function [ par, s ] = update_parameters( par, step, z, s, stats )
void update_parameters(par_t *par, int step0, double **z, accum_t *s, stats_t *stats )
{
  int step = step0+1; // adjustment to fit range 1..Ngen
  int K1 = 200;
  int Ka = 200;
  double taf = 0.95;

  int L = par->L;
  int Ni = par->Ni; //size(z,2);
  int Np = par->N;

  //tmp1 = sum( sum( z, 2) , 3 );
  double tmp1[Np];

  for (int k = 0; k < Np; k++) tmp1[k] = 0;

  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < Ni; j++)
    {
      for (int k = 0; k < Np; k++)
        tmp1[k] += z[i][k*Ni+j];
    }
  }

#if VERBOSE
  print_matrix("tmp1", tmp1, Np);
#endif

  double *tmp2 = calloc(1, par->N*par->N*sizeof(double));

  for (int i = 0; i < L; i++)
  {
    double zi[Np][Ni], ziT[Ni][Np];
    for (int j = 0; j < Ni; j++)
      for (int k = 0; k < Np; k++)
        zi[k][j] = z[i][k*Ni+j];

#if VERBOSE
    print_matrix_2d_linear("zi", (double *)zi, Np, Ni);
#endif

    for (int j = 0; j < Ni; j++)
      for (int k = 0; k < Np; k++)
        ziT[j][k] = zi[k][j];

#if VERBOSE
    print_matrix_2d_linear("ziT", (double *)ziT, Ni, Np);
#endif

    double zz[Np][Np];
    for (int ki = 0; ki < Np; ki++)
    for (int kj = 0; kj < Np; kj++)
    {
      zz[ki][kj] = 0.0;
      for (int kk = 0; kk < Ni; kk++)
        zz[ki][kj] += zi[ki][kk]*ziT[kk][kj];
    }

#if VERBOSE
    print_matrix_2d_linear("zz", (double *)zz, Np, Np);
#endif

    for (int ki = 0; ki < Np; ki++)
    for (int kj = 0; kj < Np; kj++)
      tmp2[ki*Np+kj] += zz[ki][kj];

  }

#if VERBOSE
  print_matrix_2d_linear("tmp2", tmp2, Np, Np);

  print_matrix_2d_linear("update_parameters/stats->serr", stats->sserr, L, Ni);
#endif

  //tmp3 = sum( sum( stats->sserr ) );
  double tmp3 = 0.0;
  for (int i = 0; i < L; i++)
  for (int j = 0; j < Ni; j++)
    tmp3 += stats->sserr[i*Ni+j];

  if (s->s1 == NULL) s->s1 = calloc(1, Np*sizeof(double));
  if (s->s2 == NULL) s->s2 = calloc(1, Np*Np*sizeof(double));

  if( step < K1+1 )
  {
    //s{1} = tmp1/L;
    //s{2} = tmp2/L;
    //s{3} = tmp3/L;
    for (int i = 0 ; i < Np; i++) s->s1[i] = tmp1[i]/L;
    for (int i = 0 ; i < Np*Np; i++) s->s2[i] = tmp2[i]/L;
    s->s3 = tmp3/L;
  }
  else
  {
    double g = 1.0/(step-K1);
    //s{1} = s{1} + g*( tmp1/L - s{1} );
    //s{2} = s{2} + g*( tmp2/L - s{2} );
    //s{3} = s{3} + g*( tmp3/L - s{3} );

    for (int i = 0 ; i < Np; i++) s->s1[i] = s->s1[i] + g*(tmp1[i]/L - s->s1[i]);
    for (int i = 0 ; i < Np*Np; i++) s->s2[i] = s->s2[i] + g*(tmp2[i]/L - s->s2[i]);
    s->s3 = s->s3 + g*(tmp3/L - s->s3);
  }

  double beta[Np];

  //beta  = s->s1 / Ni;
  for (int i = 0 ; i < Np; i++) beta[i] = s->s1[i]/Ni;

  double omega[Np][Np];
  double beta2[Np][Np];

  for (int i = 0; i < Np; i++)
  for (int j = 0; j < Np; j++)
  {
    beta2[i][j] = beta[i]*beta[j];
  }

  for (int i = 0; i < Np; i++)
  for (int j = 0; j < Np; j++)
    omega[i][j] = s->s2[i*Np+j] / Ni - beta2[i][j];

//  omega = s->s2 / Ni - beta*beta;
  //omega = diag(diag(omega));
  for (int i = 0; i < Np; i++)
  for (int j = 0; j < Np; j++)
    if (i != j) omega[i][j] = 0.0;

  double alpha = sqrt( s->s3/stats->N_data_tot);

  //par->beta  = beta;
  for (int i = 0; i < Np; i++)
  {
      par->beta[i] = beta[i];
  }

  if(step < Ka)
  {
    //par->omega =  max( taf*par->omega,  omega );
    //par->alpha =  max( taf*par->alpha,  alpha );

    for (int i = 0; i < Np; i++)
    for (int j = 0; j < Np; j++)
      par->omega[i*Np+j] = max(taf*par->omega[i*Np+j], omega[i][j]);

    par->alpha = max(taf*par->alpha, alpha);
  }
  else
  {
    //par->omega = omega;
    //par->alpha = alpha;
    for (int i = 0; i < Np; i++)
    for (int j = 0; j < Np; j++)
      par->omega[i*Np+j] = omega[i][j];

    par->alpha = alpha;
  }

  chol_eval(par->omega, par->omega_chol, Np); //par->omega_chol = chol_eval( par->omega );

  free(tmp2);
}
