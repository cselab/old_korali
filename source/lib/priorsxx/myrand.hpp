#ifndef MYRAND_HPP
#define MYRAND_HPP

namespace priors {

	void gsl_rand_init(int seed);

	double normal_pdf(double x, double *p);
	double normal_log_pdf(double x, double *p);
	double normal_rnd( double *p );

	double uniform_pdf(double x, double *p);
	double uniform_log_pdf(double x, double *p);
	double uniform_rnd( double *p );

	double exp_pdf(double x, double *p);
	double exp_log_pdf(double x, double *p);
	double exp_rnd( double *p );

	double gamma_pdf(double x, double *p);
	double gamma_log_pdf(double x, double *p);
	double gamma_rnd( double *p );

	double normalrand(double mu, double sigma);
	double uniformrand(double a, double b);

}

#endif//MYRAND_HPP

