#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


static gsl_vector *mu;
static gsl_matrix *sigma;

static double f_Ackley(double *x, int N);
static double f_Dixon_Price(double *x, int N);
static double f_Griewank(double *x, int N);
static double f_Levy(double *x, int N);
static double f_Perm(double *x, int N);
static double f_Perm0(double *x, int N);
static double f_Rastrigin(double *x, int N); 
static double f_Rosenbrock(double *x, int N);
static double f_Rotated_Hyper_Ellipsoid(double *x, int N);
static double f_Schwefel(double *x, int N);
static double f_Sphere(double *x, int N);
static double f_Styblinski_Tang(double *x, int N);
static double f_Sum_Of_Power(double *x, int N);
static double f_Sum_Of_Squares(double *x, int N);
static double f_Zakharov(double *x, int N);

void fitfun_initialize(char *s) {
    int i, N = atoi(s);

    mu = gsl_vector_calloc(N);
    for( i = 0; i < N; ++i)
        gsl_vector_set( mu, i, (double)i );

    sigma = gsl_matrix_calloc( N,N );

    double D0;
    for( i = 0; i < N; ++i) {
        D0 = (i%2)*0.5 + 0.02;
        gsl_matrix_set( sigma, i, i, D0 );
    }

    for( i = 1; i < N; ++i) {
        D0 = (2*(i%2) - 1) * 0.05; //change the sign
        gsl_matrix_set( sigma, i, i-1, D0 );
    }

    gsl_linalg_cholesky_decomp( sigma );
}

void fitfun_finalize() {
	gsl_vector_free(mu);
	gsl_matrix_free(sigma);
}

double fitfun(double *x, int N, void *output, int *info) {
	gsl_vector *xg 	= gsl_vector_calloc(N);
	gsl_vector *work = gsl_vector_calloc(N);

	double result;

	for( int i=0; i<N; i++)
		gsl_vector_set( xg, i, x[i] );

	gsl_ran_multivariate_gaussian_log_pdf( xg, mu, sigma, &result, work);

	gsl_vector_free(xg);
	gsl_vector_free(work);

	return result;

}


static double f_Ackley(double *x, int N) {
    int i;
    double a = 20, b = .2, c = 2.*M_PI, s1 = 0., s2 = 0.;
    for (i = 0; i < N; ++i) {
        s1 += x[i]*x[i];
        s2 += cos(c*x[i]);
    }
    return -a*exp(-b*sqrt(s1/N)) - exp(s2/N) + a + exp(1.);
}

static double f_Dixon_Price(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 1; i < N; ++i)
        s += (i+1.)*pow(2*x[i]*x[i]-x[i-1], 2);
    return pow(x[0]-1., 2) + s;
}

static double f_Griewank(double *x, int N) {
     int i;
     double s = 0., p = 1.;
     for (i = 0; i < N; ++i) {
         s += x[i]*x[i];
         p *= cos(x[i]/sqrt(1.+i));
     }
     return s/4000. - p + 1.;
 }

static double f_Levy(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 0; i < N-1; ++i) {
        s += .0625*pow(x[i]-1, 2)
            * (1.+10.*pow(sin(M_PI*.25*(3+x[i]) +1), 2));
    }
    return pow(sin(M_PI*.25*(3.+x[0])), 2) + s
        + .0625*pow(x[N-1]-1, 2)*(1+pow(sin(M_PI*.5*(3+x[N-1])), 2));
}

static double f_Perm(double *x, int N) {
    int i;
    double beta = .5;
    double s2 = 0.;
    int j;
    for (i = 0; i < N; ++i) {
        double s1 = 0.;
        for (j = 0; j < N; ++j)
            s1 += (pow(j+1, i+1)+beta)*(pow(x[j]/(j+1.), i+1) - 1.);
        s2 += s1*s1;
    }
    return s2;
}

static double f_Perm0(double *x, int N) {
    double beta = 10.;
    double s2 = 0.;
    int i, j;
    for (i = 0; i < N; ++i) {
        double s1 = 0.;
        for (j = 0; j < N; ++j)
            s1 += (j+1.+beta)*(pow(x[j], i+1) - 1./pow(j+1, i+1));
        s2 += s1*s1;
    }
    return s2;
}

static double f_Rastrigin(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 0; i < N; ++i)
        s += x[i]*x[i] - 10.*cos(2.*M_PI*x[i]);
    return 10.*N+s;
}
 
static double f_Rosenbrock(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 0; i < N-1; ++i)
        s += 100.*pow(x[i+1]-x[i]*x[i], 2) + pow(x[i]-1., 2);
    return s;
}

static double f_Rotated_Hyper_Ellipsoid(double *x, int N) {
    int i, j;
    double s = 0.;
    for (i = 0; i < N; ++i)
        for (j = 0; j <= i; ++j)
            s += x[j]*x[j];
    return s;
}

static double f_Schwefel(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 0; i < N; ++i)
        s += x[i]*sin(sqrt(fabs(x[i])));
    return 418.9829*N-s;
}

static double f_Sphere(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 0; i < N; ++i)
        s += x[i]*x[i];
    return s;
}

static double f_Styblinski_Tang(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 0; i < N; ++i)
        s += pow(x[i], 4) - 16.*x[i]*x[i] + 5.*x[i];
    return 39.16599*N + .5*s;
}

static double f_Sum_Of_Power(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 0; i < N; ++i)
        s += pow(fabs(x[i]), i+2);
    return s;
}

static double f_Sum_Of_Squares(double *x, int N) {
    int i;
    double s = 0.;
    for (i = 0; i < N; ++i)
        s += (i+1.)*x[i]*x[i];
    return s;
}

static double f_Zakharov(double *x, int N) {
    int i;
    double s1 = 0., s2 = 0.;
    for (i = 0; i < N; ++i) {
        s1 += x[i]*x[i];
        s2 += .5*(i+1.)*x[i];
    }
    return s1 + pow(s2, 2) + pow(s2, 4);
}
