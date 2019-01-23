#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics.h>

#include "tmcmc_utils.hpp"
#include "tmcmc_types.hpp"
#include "tmcmc_stats.hpp"

#define CHECK_POSDEF
#define LARGE_SCALE_POPS

double tmcmc_objlogp(double x, const double *fj, int fn, double pj, double tol)
{
    const double fjmax = gsl_stats_max(fj, 1, fn);
   
#ifdef LARGE_SCALE_POPS
    double *weight = new double[fn];
#else
    double weight[fn];
#endif

    for(int i = 0; i <fn; ++i)
        weight[i] = exp((fj[i]-fjmax)*(x-pj));

    double sum_weight = std::accumulate(weight, weight+fn, 0.0);

#ifdef LARGE_SCALE_POPS
    double *q = new double[fn];
#else
    double q[fn];
#endif

    for(int i = 0; i < fn; ++i) {
        q[i] = weight[i]/sum_weight;
    }

    double mean_q = gsl_stats_mean(q, 1, fn);
    double std_q  = gsl_stats_sd_m(q, 1, fn, mean_q);
    double coef   = pow(std_q/mean_q-tol, 2); // (TODO: is this true and not cvar? (DW))

#ifdef LARGE_SCALE_POPS
    delete[] weight;
    delete[] q;
#endif

    return coef;
}


double tmcmc_objlogp_gsl(double x, void *param) {
	fparam_t *fp = (fparam_t *) param;
	return tmcmc_objlogp(x, fp->fj, fp->fn, fp->pj, fp->tol);
}


double tmcmc_objlogp_gsl2(const gsl_vector *v, void *param) {
	double x = gsl_vector_get(v, 0);
	return tmcmc_objlogp_gsl(x, param);
}


/* simple min search, check stepwise for x < opt.Tol; if unsuccessfull increase
    range and refine steps (DW: could be optimized) */
int fzerofind(const optim_options& opt, double const *fj, int fn, double pj, 
              double objTol, double *xmin, double *fmin) {

    bool display = opt.Display;
    double tol   = opt.Tol;
    double step  = opt.Step;
    
    double x_lo = 0.0;
    double x_hi = 4.0;

    int first_try = 1;
    bool dump = 0;

    FILE *fp = NULL;
    char fname[64];

    static int counter = 0;
    
    size_t iter;
    size_t niters;


    double m, fm;
    bool converged = false;
    
    counter++;
    while (converged == false && 1e-16 < step) {
    
    if(display) printf("fminzero: x_lo %e x_hi %ei step %e\n", x_lo, x_hi, step);
    niters = (size_t) ((x_hi-x_lo)/step);

    first_try++;
    if (first_try) dump = true;

    if (dump) {
        sprintf(fname, "fzero_%03d.txt", counter);
        fp = fopen(fname, "w");
    }

    m  = 0;
    fm = std::numeric_limits<double>::max();
    double t0  = torc_gettime();

#if !defined(_USE_OPENMP_)
    for (iter = 0; iter < niters; ++iter) {
        
        double x  = x_lo + iter*step;
        double fx = tmcmc_objlogp(x, fj, fn, pj, tol);
		
        if (dump) fprintf(fp, "%.16f %.16f\n", x, fx);

		if (fx < fm) {
			fm = fx;
			m  = x;
		}
		if (fabs(fx) <= tol) {
			converged = true;
			break;
		}
    }
#else
	#pragma omp parallel 
    {
        double lm = 0;
        double lfm = std::numeric_limits<double>::max();
        #pragma omp for
        for (iter = 0; iter < niters; ++iter)
        {
            double x, fx;

            if (converged == false)
            {
                x  = x_lo + iter*Step;
                fx = tmcmc_objlogp_gsl(x, fj, fn, pj, tol);
                if (fx < lfm)
                {
                    lfm = fx;
                    lm  = x;
                }
                if (fabs(fx) <= tol) {
                    converged = true;
                    #pragma omp flush(converged)
                }
            } /* (PH: task cancellation?) */
        }
        
        #pragma omp critical
        {
            if (lfm < fm) {
                fm = lfm;
                m  = lm;
            }
        }
    }
#endif
    double t1 = torc_gettime();

    if (converged) {
        if (display) printf("fzerofind: m=%.16f fm=%.16f iter=%ld, time=%lf s\n",
                             m, fm, niters, t1-t0);
    } else {
        x_lo = m - 10*step;
        if (x_lo < 0) x_lo = 0.0;
        x_hi = m + 10*step;
        if (x_hi > 4) x_hi = 4.0;
        step *= 0.1;
        if (display) printf("fzerofind (%e): m=%.16f fm=%.16f iter=%ld, time=%lf s\n", 
                             step, m, fm, niters, t1-t0);
    }

    }

	*xmin = m; // (TODO: check if we should set this if converged == false (DW))
	*fmin = fm;

    if (dump) fclose(fp);

    return converged;
}


int fminsearch(const optim_options& opt, double const *fj, int fn, double pj,
               double objTol, double *xmin, double *fmin)
{
	const gsl_multimin_fminimizer_type *T;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	bool converged = 0;

	size_t iter     = 0; 
    size_t max_iter = opt.MaxIter;
    double tol      = opt.Tol;
	bool display    = opt.Display;
	double step     = opt.Step;
	
    int status;
	double size;

	fparam_t fp;
	fp.fj = fj; fp.fn = fn; fp.pj = pj; fp.tol = objTol;

	x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, pj);

	ss = gsl_vector_alloc (1);
	gsl_vector_set_all (ss, step); 

	minex_func.n      = 1;
	minex_func.f      = tmcmc_objlogp_gsl2;
	minex_func.params = &fp;

/*	T = gsl_multimin_fminimizer_nmsimplex;*/
	T = gsl_multimin_fminimizer_nmsimplex2;
/*	T = gsl_multimin_fminimizer_nmsimplex2rand;*/
	s = gsl_multimin_fminimizer_alloc (T, 1);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	if (display) printf ("using %s method\n", gsl_multimin_fminimizer_name (s));

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status) break;

		size   = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, tol);

		if (status == GSL_SUCCESS) {
			converged = true;
			if (display) printf ("converged to minimum at\n");
		}
		else if (fabs(s->fval) <= tol) {
			converged = true;
			status = GSL_SUCCESS;
			if (display) printf ("found minimum at\n");
		}
		if (display) printf ("%3ld x =  %.16lf f() = %.16f size = %.16f\n",
				              iter, gsl_vector_get (s->x, 0), s->fval, size);

	} while (status == GSL_CONTINUE && iter < max_iter);

	/* double-check */
	if ((converged == 1) && (fabs(s->fval) > tol)) {
		converged = false;
		if (display) printf ("fminsearch: converged but not found minimum.\n");
	}

	if (converged) {
		*fmin = s->fval;
		*xmin = gsl_vector_get(s->x, 0);
	} else {
		*fmin = 0;
		*xmin = 0.0;
	}

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return converged;
}

int fmincon(const optim_options& opt, const double *fj, int fn, double pj, 
            double objTol, double *xmin, double *fmin) {
	
    int status;
	size_t iter     = 0; 
    size_t max_iter = opt.MaxIter;
    double tol      = opt.Tol;
	bool display    = opt.Display;

	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	gsl_function F;
	gsl_vector *x;
    
    double x_lo = 0.0, x_hi = 4.0;    /* input */
	double m = 0.5, fm = 0.0;
	bool converged = 0;
	int i;

	x = gsl_vector_alloc (1);

	fparam_t fp;
	fp.fj = fj; fp.fn = fn; fp.pj = pj; fp.tol = objTol;

	F.function = tmcmc_objlogp_gsl;
	F.params = &fp;

	T = gsl_min_fminimizer_brent;
/*	T = gsl_min_fminimizer_goldensection;*/
/*	T = gsl_min_fminimizer_quad_golden;*/
	s = gsl_min_fminimizer_alloc (T);

	double f_lo = tmcmc_objlogp_gsl(x_lo, &fp);
	double f_hi = tmcmc_objlogp_gsl(x_hi, &fp);
	if (f_lo < f_hi) {
		m  = x_lo;
		fm = f_lo;
	} else {
		m  = x_hi;
		fm = f_hi;
	}

	for (i = 0; i < max_iter; i++) {
		double x = x_lo + i*(x_hi-x_lo)/max_iter;
		double fx = tmcmc_objlogp_gsl(x, &fp);
		if (fx < fm) {
			m = x;
			fm = fx;
		}
	}

	if (fabs(fm) <= tol) {
		converged = true;
		gsl_vector_free(x);
		gsl_min_fminimizer_free (s);
		if (display) printf("fmincon: early return with m = %.16f fm = %.16f\n", m, fm);
		return converged;
	}

	if ((fm < f_lo) && (fm < f_hi)) {
		if (display) printf("fmincon: initialized with %d tries and m = %f (fm = %f)\n", i, m, fm);
	} else {
		if (display) printf("failed to initialize fmincon (%.16f, %.16f)!\n", f_lo, f_hi);
		return 0;
	}

	gsl_min_fminimizer_set (s, &F, m, x_lo, x_hi);

	if (display) {
		printf ("using %s method\n", gsl_min_fminimizer_name (s));
		printf ("%5s [%18s, %18s] %18s %18s %18s\n", "iter", "lower", "upper", "min", "fmin", "err(est)");
		printf ("%5zu [%.16f, %.16f] %.16f %.16f %.16f\n", iter, x_lo, x_hi, m, fm, x_hi - x_lo);
	}

	do {
		iter++;
		status = gsl_min_fminimizer_iterate (s);

		m    = gsl_min_fminimizer_x_minimum (s);
		x_lo = gsl_min_fminimizer_x_lower (s);
		x_hi = gsl_min_fminimizer_x_upper (s);

		status = gsl_min_test_interval (x_lo, x_hi, tol, tol);
		if (status == GSL_SUCCESS) {
			if (display) printf ("Converged:\n");
			converged = true;
		}
		else if (fabs(gsl_min_fminimizer_f_minimum(s)) <= tol) {
			converged = true;
			status = GSL_SUCCESS;
			if (display) printf ("found minimum at\n");
		}

		if (display)
			printf ("%5zu [%.16f, %.16f] %.16f f()=%.16f %.16f\n",
				iter, x_lo, x_hi, m, gsl_min_fminimizer_f_minimum(s), x_hi - x_lo);

	} while (status == GSL_CONTINUE && iter < max_iter);

	/* double-check */
	if ((converged == true) && (fabs(gsl_min_fminimizer_f_minimum(s)) > tol)) {
		converged = false;
		if (display) printf ("converged but not found minimum.\n");
	}

	if (converged) {
		converged = 1;
		gsl_vector_set (x, 0, m);
		*fmin = tmcmc_objlogp_gsl(m, &fp);
		*xmin = m;
	} else {
		*fmin = 0;
		*xmin = 0.0;
	}

	gsl_vector_free(x);
	gsl_min_fminimizer_free (s);

	return converged;
}
