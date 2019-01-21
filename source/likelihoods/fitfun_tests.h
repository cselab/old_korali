#ifndef FITFUN_TESTS_H
#define FITFUN_TESTS_H

typedef double (*fitfun_t)(double*, int);

static fitfun_t my_fitfun;

double fitfun(double *x, int N, void *output, int *info);

void fitfun_initialize_simple(const char *func);

static double f_multivariate_gaussian(double *x, int N);
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

#endif // FITFUN_TESTS_H 
