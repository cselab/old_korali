#include "d_loglike.c"
//function [  FIM, stats ] = estimate_FIM( z, par, stats, model )

double **estimate_FIM(double **z0, par_t *par, stats_t *stats, model_t *model, MCMC_t *MCMC, int *outNs, int *outM)
{
  //int Np = size(z,1);
  int Np = par->N;
  int Ni = par->Ni;
  int L  = par->L;
  int Nsteps = (MCMC->steps[0]*1 + MCMC->steps[1]*Np + MCMC->steps[2]*Np) * MCMC->outer_steps;

  printf("Nsteps = %d\n", Nsteps);
  //exit(1);

  int M;
  if (strcmp(model->error, "ind") == 0)
    M = 2*Np;
  else
    M = 2*Np + 1;

  //int Ns0 = size(stats->samples,1);
  int Ns = Nsteps; //*L*Ni*Np;
  int Nb = ceil(Ns/2);
  Ns = Ns - Nb;

  printf("---> Ns = %d\n", Ns);
  *outNs = Ns;
  *outM = M;
  //D = zeros(Ns,M);
  //H = zeros(Ns,M,M);
  //G = zeros(Ns,M,M);

  double **D = (double **)calloc(1, Ns*sizeof(double *));
  for (int i = 0; i < Ns; i++) D[i] = (double *)calloc(1, M*sizeof(double));

  double **H = (double **)calloc(1, Ns*sizeof(double *));
  for (int i = 0; i < Ns; i++) H[i] = (double *)calloc(1, M*M*sizeof(double));

  double **G = (double **)calloc(1, Ns*sizeof(double *));
  for (int i = 0; i < Ns; i++) G[i] = (double *)calloc(1, M*M*sizeof(double));

  double **FIM = (double **)calloc(1, Ns*sizeof(double *));
  for (int i = 0; i < Ns; i++) FIM[i] = (double *)calloc(1, M*M*sizeof(double));

  for (int i=0; i<Ns; i++)
  {
    /*
    // info, L=1!
    //stats_new->samples =   (double *)calloc(1, Nsteps*L*Ni*Np*sizeof(double));;
    //stats->samples(cnt_s,:,i,l) = z_loc;
    //for (int c = 0; c < Np; c++)
    //    stats->samples[cnt_s*L*Ni*Np+l*Ni*Np+i*Np+c] = z_loc[c];        // Nsteps*L*Ni*Np
    */

    //z = squeeze( stats->samples(Nb+i,:,:,1) );
    //if(Np==1)
    //    z=z';
    //end
    int l = 0;
    double z[Ni*Np];
    for (int j = 0; j < Ni; j++)
    for (int k = 0; k < Np; k++)
      z[j*Np+k] = stats->samples[(Nb+i)*L*Ni*Np+l*Ni*Np+j*Np+k];


    //sserr = stats->sserr_tot(Nb+i,:);
    double sserr[Ni];

    for (int j = 0; j < Ni; j++)
        sserr[j] = stats->sserr_tot[(Nb+i)*L*Ni+l*Ni+j];


    //[ d, h ] = d_loglike(  z, sserr, par, model );
    double d[M];
    double h[M*M];

    d_loglike(z, sserr, par, model, d, h, M);

    //D(i,:)   =  d;
    //H(i,:,:) =  h;
    //G(i,:,:) =  d'*d;  // '
    for (int j = 0; j < M  ; j++) D[i][j] = d[j];
    for (int j = 0; j < M*M; j++) H[i][j] = h[j];

    for (int j = 0; j < M; j++)
    for (int k = 0; k < M; k++)
      G[i][j*M+k] = d[j]*d[k];

    double tmp1[M*M];
    double tmp2[M*M];
    double tmp3_v[M], tmp3[M*M];

    memset(tmp1  , 0, M*M*sizeof(double));
    memset(tmp2  , 0, M*M*sizeof(double));
    memset(tmp3_v, 0, M*sizeof(double));
    memset(tmp3  , 0, M*M*sizeof(double));

    //tmp1 = squeeze( mean( H(1:i,:,:) , 1 ) );
    //tmp2 = squeeze( mean( G(1:i,:,:) , 1 ) );
    //tmp3 = mean( D(1:i,:) , 1 );
    //tmp3 = tmp3'*tmp3;  // '

    for (int p = 0; p <= i; p++)
    {
      for (int j = 0; j < M; j++)
      for (int k = 0; k < M; k++)
      	tmp1[j*M+k] += H[p][j*M+k];
    }

    for (int j = 0; j < M; j++)
    for (int k = 0; k < M; k++)
      tmp1[j*M+k] = tmp1[j*M+k] / (i+1);

    for (int p = 0; p <= i; p++)
    {
      for (int j = 0; j < M; j++)
      for (int k = 0; k < M; k++)
        tmp2[j*M+k] += G[p][j*M+k];
    }

    for (int j = 0; j < M; j++)
    for (int k = 0; k < M; k++)
      tmp2[j*M+k] = tmp2[j*M+k] / (i+1);


    for (int p = 0; p <= i; p++)
    {
  	for (int j = 0; j < M; j++)
  		tmp3_v[j] += D[p][j];
    }

    for (int j = 0; j < M; j++)
      tmp3_v[j] = tmp3[j] / (i+1);

    for (int j = 0; j < M; j++)
    for (int k = 0; k < M; k++)
      tmp3[j*M+k] = tmp3_v[j]*tmp3_v[k];

    //FIM(i,:,:) =  - ( tmp1 + tmp2  - tmp3 );
    for (int j = 0; j < M; j++)
    for (int k = 0; k < M; k++)
      FIM[i][j*M+k] = - (tmp1[j*M+k] + tmp2[j*M+k] - tmp3[j*M+k]);

    char msg[64];
    sprintf(msg, "FIM[%d]", i);
    print_matrix_2d_linear(msg, FIM[i], M, M);
  }

  return FIM;
}
