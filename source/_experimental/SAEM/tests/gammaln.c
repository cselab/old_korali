#include <gsl/gsl_sf_gamma.h>
#include <stdio.h>


int main()
{
	printf("%lf\n", gsl_sf_lngamma(4.0));
	printf("%lf\n", gsl_sf_lngamma(5.0));

	return 0;
}
