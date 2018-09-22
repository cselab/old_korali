#include <stdio.h>
#include <gsl/gsl_linalg.h>
//function z =  logp_eta_eval( x, sigma_chol )
double logp_eta_eval(double *x, double *sigma_chol, int N)
{
	//% log probability of zero mean normal distribution
	//%
	//% columns of x must be equal to columns of sigma_chol

	int d = N;

	double log_det_omega = 0;

	for (int i = 0; i < d; i++)
		log_det_omega += log(sigma_chol[i*d+i]);

	double c = -0.5*d*log(2*M_PI) - log_det_omega;

//	printf("c = %lf\n", c);

#if 0
	one_d = size(sigma_chol,1) < 2;
	utri  = istriu(sigma_chol);

	if( one_d || utri )
	{
		q = x / sigma_chol ;
	}
	else
	{
		q =  sigma_chol \ x'; //'
		q = q';	//'
	}
#else

	double q[d];

#if 0
	for (int i = 0; i < d; i++)
		q[i] = x[i]/sigma_chol[i*d+i];

#else
	gsl_matrix_view m = gsl_matrix_view_array (sigma_chol, d, d);

	gsl_vector_view b = gsl_vector_view_array (x, d);

	gsl_vector *sol = gsl_vector_alloc (d);

	int s;

	gsl_permutation * p = gsl_permutation_alloc (d);
	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, sol);

#if VERBOSE
	printf ("sol = \n");
	gsl_vector_fprintf (stdout, sol, "%g");
#endif

	for (int i = 0; i < d; i++)
		q[i] = sol->data[i];

	gsl_permutation_free (p);
	gsl_vector_free (sol);
#endif

#endif

#if VERBOSE
	for (int i = 0; i < d; i++)
		printf("q[%d]=%lf\n", i, q[i]);
#endif

	//z = c - 0.5*sum(q.^2,2);
	double sumq2 = 0;
	for (int i = 0; i < d; i++) sumq2 += q[i]*q[i];

	double z = c - 0.5*sumq2;

	return z;
}
