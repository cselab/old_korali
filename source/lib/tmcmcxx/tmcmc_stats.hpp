#ifndef TMCMC_STATS_HPP
#define TMCMC_STATS_HPP

typedef struct fparam_s {
	const double *fj;
	int           fn;
	double        pj;
	double        tol;
} fparam_t;

void calculate_statistics(double flc[], unsigned int n, int nselections, 
                          int gen, unsigned int sel[]);

#endif // TMCMC_STATS_HPP
