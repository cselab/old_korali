#ifndef MYRAND_HPP
#define MYRAND_HPP

#include <cstddef>

namespace priors {

	//void gsl_rand_init(int seed);
    void spmd_gsl_rand_init(int seed);

	double normal_pdf(double x, double *p);
	double normal_log_pdf(double x, double *p);
	double normal_rnd( double *p );

	double uniform_pdf(double x, double *p);
	double uniform_log_pdf(double x, double *p);
	double uniform_rnd( double *p );
	double uniformrand(double a, double b);
	
    double exp_pdf(double x, double *p);
	double exp_log_pdf(double x, double *p);
	double exp_rnd( double *p );

	double gamma_pdf(double x, double *p);
	double gamma_log_pdf(double x, double *p);
	double gamma_rnd( double *p );

	double normalrand(double mu, double sigma);

    //TODO: prefer this in rand file (DW)
    int mvnrnd(double *mean, double *sigma, double *out, int N);
    void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[]);

}

#endif//MYRAND_HPP

