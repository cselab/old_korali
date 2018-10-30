//function [ D, H ] = d_loglike(  z, sserr, par, model )

/*
z =

       191.91       178.73       231.17       208.99       216.91       174.05       194.41       219.37       180.58       210.66
       48.319       36.795       56.386       30.507       45.893        43.88       24.993        33.65       34.641        48.74
       1.0866       1.3358       1.0586      0.89236       0.9229       1.2026      0.96497      0.93139       1.1507      0.93681


sserr =

        314.5       534.72       584.27       572.58       420.36       665.86       882.49       560.68       343.12       343.83


H =

    -0.031291            0            0    0.0028673            0            0            0
            0     -0.11597            0            0     0.048535            0            0
            0            0      -524.34            0            0      -361.22            0
    0.0028673            0            0    -0.069447            0            0            0
            0     0.048535            0            0     -0.24508            0            0
            0            0      -361.22            0            0      -1242.4            0
            0            0            0            0            0            0       -24.35


D =

    -0.025629     -0.22534       24.943     0.040905     0.040656        8.919       7.3119
*/



void d_loglike(double *z, double *sserr, par_t *par, model_t *model, double *D_, double *H_, int M)
{
  //double z[Ni*Np];  we have Ni z's of size Np
  //double sserr[Ni];

  //%
  //%   D: gradient
  //%   H: Hessian

  int N = par->N;   // N is Np in estimate_FIM()
  int Np = N;
  int Nd = par->Nd;
  int Ni = par->Ni;

  double *beta  = par->beta;

  beta[0] = 201.5;
  beta[1] = 42.323;
  beta[2] = 1.0007;

  double alpha = par->alpha;

  double omega[N], omega2[N], omega3[N], omega4[N];

#if DEBUG
  print_matrix_2d_linear("d_loglike/z", z, Ni, Np);
  print_matrix("d_loglike/sserr", sserr, Ni);
#endif
  //omega = diag(par.omega);
  //omega  = sqrt(omega);
  //omega2 = omega.^2;
  //omega3 = omega.^3;
  //omega4 = omega.^4;
  for (int i = 0; i < N; i++)
  {
    omega[i] = sqrt(par->omega[i*N+i]);
    omega2[i] = pow(omega[i], 2.0);
    omega3[i] = pow(omega[i], 3.0);
    omega4[i] = pow(omega[i], 4.0);
  }

  //%% ========================================================================

  /*
  dif  = bsxfun( @minus, z, beta );
  sum1 = sum( dif,    2 );
  sum2 = sum( dif.^2, 2 );
  */

  //dif = bsxfun( @minus, z, par.beta );
  double dif[Ni*Np];
  for (int i = 0; i < Ni; i++)
    for (int j = 0; j < Np; j++)
       dif[i*Np+j] = z[i*Np+j] - beta[j];

#if DEBUG
  print_matrix_2d_linear("d_loglike/diff", dif, Ni, Np);
#endif

  double sum1[Np], sum2[Np];
  for (int j = 0; j < Np; j++)
  {
    sum1[j] = 0;
    sum2[j] = 0;
    for (int i = 0; i < Ni; i++)
    {
      sum1[j] += dif[i*Np+j];
      sum2[j] += dif[i*Np+j]*dif[i*Np+j];
    }
  }

  //% gradient
  // D = zeros(2*N,1);
  //D(1:N)     = sum1./omega2;
  //D(N+1:2*N) = -Ni./omega + sum2./omega3;
  double *D = (double *)calloc(1, M*sizeof(double)); // = zeros(2*N,1);
  for (int i = 0; i < N; i++) D[i] = sum1[i]/omega2[i];
  for (int i = 0; i < N; i++) D[i+N] = -Ni/omega[i] + sum2[i]/omega3[i];

  //% hessian
  double tmp1[N], tmp2[N], tmp3[N];

  //tmp1 = -Ni./omega2;
  //tmp1 = diag(tmp1);
  for (int i = 0; i < N; i++) tmp1[i] = -Ni/omega2[i];

  //tmp2 = -2*sum1./omega3;
  //tmp2 = diag(tmp2);
  for (int i = 0; i < N; i++) tmp2[i] = -2*sum1[i]/omega3[i];

  //tmp3 = Ni./omega2 - 3*sum2./omega4;
  //tmp3 = diag(tmp3);
  for (int i = 0; i < N; i++) tmp3[i] = Ni/omega2[i] -3*sum2[i]/omega4[i];

  double *H = (double *)calloc(1, M*M*sizeof(double));
  //H = [ tmp1 , tmp2 ; tmp2 , tmp3 ];
  for (int i = 0; i < N; i++) {
    H[i*M+i] = tmp1[i];
    H[(i+N)*M+i] = tmp2[i];
    H[i*M+(N+i)] = tmp2[i];
    H[(i+N)*M+(N+i)] = tmp3[i];
  }


  //% handle fixed parameters
  //% reallocate arrays
  double sum_sserr = compute_sum(sserr, Ni);
  if (strcmp(model->error, "common") == 0) {
    //D(end+1) = -Nd/alpha + sum(sserr)/(alpha^3);
    D[M-1] = -Nd/alpha + sum_sserr/pow(alpha, 3);

    //H(end+1,end+1) = Nd/(alpha^2) - 3*sum(sserr)/(alpha^4);
    H[(M-1)*M+(M-1)] = Nd/pow(alpha,2) - 3*sum_sserr/pow(alpha,4);
  }

  //D = D';

#if DEBUG
  print_matrix_2d_linear("d_loglike/H", H, M, M);
  print_matrix("d_loglike/D", D, M);
#endif

  memcpy(H_, H, M*M*sizeof(double));
  memcpy(D_, D, M*sizeof(double));

  free(H);
  free(D);
}
