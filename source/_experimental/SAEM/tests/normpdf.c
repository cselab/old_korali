#include <stdio.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>


double normpdf(double x)
{
  double pdf = gsl_ran_ugaussian_pdf(x);
  return pdf;
}

double normcdf(double x)
{
  double cdf = gsl_cdf_ugaussian_P (x);
  return cdf;
}


int
main (void)
{
  double x = 0.45;

  double pdf = normpdf(x);
  printf ("pdf(%f) = %f\n", x, pdf);

  double cdf = normcdf(x);
  printf ("cdf(%f) = %f\n", x, cdf);

  return 0;
}
