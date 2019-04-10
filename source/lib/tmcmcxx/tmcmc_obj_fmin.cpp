#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics.h>

#include "tmcmc_utils.hpp"
#include "tmcmc_types.hpp"
#include "tmcmc_obj_fmin.hpp"

#define LARGE_SCALE_POPS


namespace tmcmc
{

double tmcmc_objlogp(double x, const double *fj, int fn, double pj, double zero)
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
    double cov2   = pow(std_q/mean_q-zero, 2);

#ifdef LARGE_SCALE_POPS
    delete[] weight;
    delete[] q;
#endif

    return cov2;
}


double tmcmc_objlogp_gsl(double x, void *param)
{
    fparam_t *fp = (fparam_t *) param;
    return tmcmc_objlogp(x, fp->fj, fp->fn, fp->pj, fp->tol);
}


double tmcmc_objlogp_gsl2(const gsl_vector *v, void *param)
{
    double x = gsl_vector_get(v, 0);
    return tmcmc_objlogp_gsl(x, param);
}


int fmincon(const double *fj, int fn, double pj, double objTol,
            double *xmin, double *fmin, const optim_options& opt)
{

    int status;
    size_t iter     = 0;
    size_t max_iter = opt.MaxIter;
    double tol      = opt.Tol;
    bool display    = opt.Display;
    double x_lo     = opt.LowerBound;
    double x_hi     = opt.UpperBound;

    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function F;
    gsl_vector *x;

    double m  = 0.5*(x_hi-x_lo);
    double fm = 0.0;

    bool converged = false;

    x = gsl_vector_alloc (1);

    fparam_t fp;
    fp.fj = fj;
    fp.fn = fn;
    fp.pj = pj;
    fp.tol = objTol;

    F.function = tmcmc_objlogp_gsl;
    F.params = &fp;

    // SELECT ONE MINIMIZER STRATEGY
    T = gsl_min_fminimizer_brent;
    /*	T = gsl_min_fminimizer_goldensection;*/
    /*	T = gsl_min_fminimizer_quad_golden;; */
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

    size_t i;
    for (i = 1; i < max_iter; ++i) {
        double x = x_lo + i*(x_hi-x_lo)/max_iter;
        double fx = tmcmc_objlogp_gsl(x, &fp);
        if (fx < fm) {
            m  = x;
            fm = fx;
        }
    }

    if (fm <= tol) {
        converged = true;
        gsl_vector_free(x);
        gsl_min_fminimizer_free (s);
        if (display)
            printf("fmincon: Early return with m = %.16f fm = %.16f\n (%dtries).", m, fm, i);
        return converged;
    }

    if ((fm < f_lo) && (fm < f_hi)) {
        if (display)
            printf("fmincon: Initialized with %d tries and m = %f (fm = %f)\n", i, m, fm);
    } else {
        if (display)
            printf("fmincon: Failed to initialize (%.16f, %.16f) (%d tries)!\n", f_lo, f_hi, i);
        return 0;
    }

    gsl_min_fminimizer_set (s, &F, m, x_lo, x_hi);

    if (display) {
        printf ("fmincon: Using %s method\n", gsl_min_fminimizer_name (s));
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
            if (gsl_min_fminimizer_f_minimum(s) <= tol) {
                if (display) printf ("fmincon: Converged:\n");
                converged = true;
            } else {
                if (display) printf ("fmincon: Converged but did not find minimum:\n");
            }
        } else if (gsl_min_fminimizer_f_minimum(s) <= tol) {
            converged = true;
            status = GSL_SUCCESS;
            if (display) printf ("fmincon: Did not converge but found minimum at\n");
        }

        if (display)
            printf ("%5zu [%.16f, %.16f] %.16f f()=%.16f %.16f\n",
                    iter, x_lo, x_hi, m, gsl_min_fminimizer_f_minimum(s), x_hi - x_lo);

    } while (status == GSL_CONTINUE && iter < max_iter);

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


int fminsearch(double const *fj, int fn, double pj, double objTol,
               double *xmin, double *fmin, const optim_options& opt)
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
    fp.fj = fj;
    fp.fn = fn;
    fp.pj = pj;
    fp.tol = objTol;

    x = gsl_vector_alloc (1);
    gsl_vector_set (x, 0, pj);

    ss = gsl_vector_alloc (1);
    gsl_vector_set_all (ss, step);

    minex_func.n      = 1;
    minex_func.f      = tmcmc_objlogp_gsl2;
    minex_func.params = &fp;

    // SELECT ONE MINIMIZER STRATEGY
    /*	T = gsl_multimin_fminimizer_nmsimplex;*/
    T = gsl_multimin_fminimizer_nmsimplex2;
    /*	T = gsl_multimin_fminimizer_nmsimplex2rand;*/
    s = gsl_multimin_fminimizer_alloc (T, 1);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    if (display) printf ("fminsearch: Using %s method\n", gsl_multimin_fminimizer_name (s));

    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status) break;

        size   = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, tol);

        if (status == GSL_SUCCESS) {
            if(s->fval <= tol) {
                converged = true;
                if (display) printf ("fminsearch: Converged to minimum at\n");
            } else {
                converged = false;
                if (display) printf ("fminsearch: Converged but did not find minimum.\n");
            }
        } else if (s->fval <= tol) {
            converged = true;
            status = GSL_SUCCESS;
            if (display)
                printf ("fminsearch: NOT converged but found minimum at\n");
        } else {
            converged = false;
            if (display)
                printf("fminsearch: NOT converged and did not find minimum.\n");
        }
        if (display) printf ("%3ld x =  %.16lf f() = %.16f size = %.16f\n",
                                 iter, gsl_vector_get (s->x, 0), s->fval, size);

    } while (status == GSL_CONTINUE && iter < max_iter);

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

/* simple min search, check stepwise for x < opt.Tol; if unsuccessfull
 * adjust range and refine steps (DW: could be optimized) */
int fzerofind(double const *fj, int fn, double pj, double objTol,
              double *xmin, double *fmin, const optim_options& opt)
{

    bool display = opt.Display;
    bool dump    = opt.Zdump;
    double tol   = opt.Tol;
    double step  = opt.Step;
    double x_lo  = opt.LowerBound;
    double x_hi  = opt.UpperBound;

    int first_try = true;

    FILE *fp = NULL;
    char fname[64];

    static int counter = 0;

    size_t iter;
    size_t niters;

    double min  = 0;
    double fm = std::numeric_limits<double>::max();

    bool converged = false;

    counter++;
    while (converged == false && 1e-16 < step) {

        if(display) printf("fzerofind: x_lo %e x_hi %ei step %e\n", x_lo, x_hi, step);
        niters = (size_t) ((x_hi-x_lo)/step);

        if (dump && first_try) {
            first_try = false;
            sprintf(fname, "fzero_%03d.txt", counter);
            fp = fopen(fname, "w");
            fprintf(fp, "%-19s%-19s\n","x","fx=(cov-TolCOV)^2");
        }

        double t0  = torc_gettime();

#ifndef _USE_OPENMP_
        for (iter = 0; iter < niters; ++iter) {
            double x  = x_lo + iter*step;
            double fx = tmcmc_objlogp(x, fj, fn, pj, objTol);

            if (dump) fprintf(fp, "%.16f %.16f\n", x, fx);

            if (fx < fm) {
                fm  = fx;
                min = x;
            }
            if (fx <= tol) {
                converged = true;
                break;
            }
        }
#else
        #pragma omp parallel
        {
            double lmin  = 0;
            double lfmin = std::numeric_limits<double>::max();
            #pragma omp for
            for (iter = 0; iter < niters; ++iter) {
                double x, fx;

                if (converged == false) {
                    x  = x_lo + iter*Step;
                    fx = tmcmc_objlogp_gsl(x, fj, fn, pj, tol);
                    if (fx < lfmin) {
                        lfmin = fx;
                        lmin  = x;
                    }
                    if (fx <= tol) {
                        converged = true;
                        #pragma omp flush(converged)
                    }
                } /* (PH: task cancellation?) */
            }

            #pragma omp critical
            {
                if (lfmin < fm) {
                    fm  = lfmin;
                    min = lmin;
                }
            }
        }
#endif
        double t1 = torc_gettime();

        if (converged) {
            if (display) printf("fzerofind converged: m=%.16f fm=%.16f iter=%ld, time=%lf s\n",
                                    min, fm, niters, t1-t0);
        } else {
            x_lo = min - 10*step;
            if (x_lo < 0) x_lo = opt.LowerBound;
            x_hi = min + 10*step;
            if (x_hi > 4) x_hi = opt.UpperBound;
            step *= 0.1;
            if (display) printf("fzerofind (%e): m=%.16f fm=%.16f iter=%ld, time=%lf s\n",
                                    step, min, fm, niters, t1-t0);
        }

    }

    *xmin = min;
    *fmin = fm;

    if (dump) fclose(fp);

    return converged;
}

} //namespace tmcmc
