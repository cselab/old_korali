#include <gsl/gsl_sf_gamma.h>

#ifdef TESTING
#include "../auxil.c"
/*
x = [209.5203 42.2562 0.8914];
v = 5;
mu = [209.5203 42.1791 0.9062];
sigma = [0.0001 0.4041 0.0169];
*/
#endif

void log_tpdf(double *xp, int Np, double v, double *mu, double *sigma, double *y)
{
  //%
  //% logarithm of location-scale t-student distribution
  //% adapted  from tpdf.m

#if VERBOSE
  //x,v,mu,sigma
  print_matrix("xp", xp, Np);
  print_matrix("mu", mu, Np);
  print_matrix("sigma", sigma, Np);
#endif

  //x = (x-mu)./sigma;
  double x[Np];
  for (int i = 0; i < Np; i++) x[i] = (xp[i]-mu[i])/sigma[i];

  //% Use gammaln function to avoid overflows.
  double term1[Np], term2[Np], term3[Np];

  for (int i = 0; i < Np; i++)
  {
    //term1 = gammaln((v + 1) / 2) - gammaln(v/2);
    //term2 = -0.5*log(v*pi)-log(sigma);
    //term3 = -0.5*(v + 1) .* log( (1 + (x.^2)./v) );

    term1[i] = gsl_sf_lngamma((v + 1) / 2) - gsl_sf_lngamma(v/2);
    term2[i] = -0.5*log(v*M_PI)-log(sigma[i]);
    term3[i] = -0.5*(v + 1) * log( (1 + (x[i]*x[i])/v) );
  }

#if VERBOSE
  print_matrix("term1", term1, Np);
  print_matrix("term2", term2, Np);
  print_matrix("term3", term3, Np);
#endif

  for (int i = 0; i < Np; i++)
  {
    //y =  term1 + term2 + term3
    y[i] =  term1[i] + term2[i] + term3[i];
  }

#if VERBOSE
  print_matrix("y", y, Np);
#endif
}


#ifdef TESTING
int main()
{
	#define Np 3
	double x[Np] = {209.5203, 42.2562, 0.8914};
	double v = 5;
	double mu[Np] = {209.5203, 42.1791, 0.9062};
	double sigma[Np] = {0.0001, 0.4041, 0.0169};
	double yref[Np] = {8.2417, -0.0843, 2.6837};


	double y[Np];

	log_tpdf(x, Np, v, mu, sigma, y);

	for (int i = 0; i < Np; i++) printf("y[%d] = %lf\n", i, y[i]);

	for (int i = 0; i < Np; i++) printf("y[%d] = %lf\n", i, yref[i]);

	return 0;
}
#endif
