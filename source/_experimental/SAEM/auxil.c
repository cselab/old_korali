/*
 *  auxil.c
 *
 *  Created by CSE Lab, D-MAVT. 
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */


#include <stdio.h>
#if defined(_USE_TORC_)
#include <torc.h>
#else
int torc_i_worker_id() { return 0; }
int torc_i_num_workers() { return 1; }
int torc_node_id() { return 0; }
#endif

#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

/**********************************************/
/* Timer routine */
/**********************************************/

#include <sys/time.h>

double gettime()
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}


/**********************************************/
/* Helper routines */
/**********************************************/

void print_matrix(char *title, double *v, int n)
{
	int i;

	printf("\n%s =\n\n", title);
	for (i = 0; i < n; i++) {
		printf("   %12.6lf\n", v[i]);
	}
	printf("\n");
}

void print_matrix_i(char *title, int *v, int n)
{
	int i;

	printf("\n%s =\n\n", title);
	for (i = 0; i < n; i++) {
		printf("  %8d\n", v[i]);
	}
	printf("\n");
}

void print_matrix_2d(char *title, double **v, int n1, int n2)
{
	int i, j;

/*	if (!display) return;*/

	printf("\n%s =\n\n", title);
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			printf("   %12.6lf", v[i][j]);
//			printf("   %20.15lf", v[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_matrix_2d_linear(char *title, double *v, int n1, int n2)
{
	int i, j;

/*	if (!display) return;*/

	printf("\n%s =\n\n", title);
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			printf("   %12.6lf", v[i*n2+j]);
//			printf("   %20.15lf", v[i*n2+j]);
		}
		printf("\n");
	}
	printf("\n");
}

void fprint_matrix_1d(FILE *fp, char *title, double *v, int n)
{
	int i;

	if (fp == stdout)
		fprintf(fp, "\n%s =\n\n", title);
	for (i = 0; i < n; i++) {
		fprintf(fp, "%12.4lf ", v[i]);
	}
	fprintf(fp, "\n");
}

void fprint_matrix_2d(FILE *fp, char *title, double **v, int n1, int n2)
{
	int i, j;

	if (fp == stdout)
		fprintf(fp, "\n%s =\n\n", title);
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			fprintf(fp, "   %20.15lf", v[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}

double compute_max(double *v, int n)
{
	int i;
	double vmax = v[0];
	for (i = 1; i < n; i++)
		if (v[i] > vmax) vmax = v[i];

	return vmax;
}

double compute_min(double *v, int n)
{
	int i;
	double vmin = v[0];
	for (i = 1; i < n; i++)
		if (v[i] < vmin) vmin = v[i];

	return vmin;
}

int compute_min_idx_i(int *v, int n)
{
	int i;
	double vmin = v[0];
	int idx = 0;

	for (i = 1; i < n; i++)
		if (v[i] < vmin) {
			vmin = v[i];
			idx = i;
		}

	return idx;
}

double compute_sum(double *v, int n)
{
	int i;
	double s = 0;
	for (i = 0; i < n; i++) s += v[i];

	return s;
}

double compute_mean(double *v, int n)
{
	int i;
	double s = 0;
	for (i = 0; i < n; i++) s += v[i];

	return s/n;
}

double compute_std(double *v, int n, double mean)
{
	int i;
	double s = 0;
	for (i = 0; i < n; i++) s += pow(v[i]-mean,2);

	return sqrt(s/(n-1));
}


/**********************************************/
/* Random number generators */
/**********************************************/

const gsl_rng_type *T;
gsl_rng **r;

void gsl_rand_init(int seed)
{
	int i, local_workers = torc_i_num_workers();
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = (gsl_rng **)malloc(local_workers*sizeof(gsl_rng *));
	for (i = 0; i < local_workers; i++) {
		r[i] = gsl_rng_alloc (T);
	}

	if (seed == 0) seed = time(0);
	for (i = 0; i < local_workers; i++) {
		gsl_rng_set(r[i], seed+i+local_workers*torc_node_id());
	}
}


void gsl_rand_finalize()
{
	int i, local_workers = torc_i_num_workers();

	for (i = 0; i < local_workers; i++) {
		gsl_rng_free(r[i]);
	}
	free(r);
}

/* normal distribution random number N(mu,rho^2)*/
double normalrand(double mu, double var)
{
	double res;

	int me = torc_i_worker_id();
	res = mu + gsl_ran_gaussian(r[me], var);

	return res;

/*	return mu + gsl_ran_gaussian(r, var);*/
}

/* uniform (flat) distribution random number between a and b */
double uniformrand(double a, double b)
{
	double res;

	int me = torc_i_worker_id();
	res = gsl_ran_flat(r[me], a, b);

	return res;

/*	return gsl_ran_flat(r, a, b);*/
}


/* random numbers chosen from Student's t distribution with nu degrees of freedom */
double trnd(double nu)
{
	double res;

	int me = torc_i_worker_id();
	res = gsl_ran_tdist(r[me], nu);

	return res;

/*	return gsl_ran_flat(r, a, b);*/
}



void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[])
{
	int me = torc_i_worker_id();
	gsl_ran_multinomial (r[me], K, N, q, nn);

	return;
}


void shuffle(int *perm, int N)
{
	int i;
	gsl_permutation * p = gsl_permutation_alloc (N);

	gsl_permutation_init (p);
#if VERBOSE
	gsl_permutation_fprintf (stdout, p, " %u");
	printf ("\n");
#endif
	int me = torc_i_worker_id();
	gsl_ran_shuffle (r[me], p->data, N, sizeof(size_t));
#if VERBOSE
	printf (" random permutation:");
	gsl_permutation_fprintf (stdout, p, " %u");
#endif
	for (i = 0; i <	N; i++)	perm[i] = p->data[i];
	gsl_permutation_free (p);
}


/**********************************************/
/* Multivariate Normal density function */
/**********************************************/

/*
 *  Multivariate Normal density function and random number generator
 *  Multivariate Student t density function and random number generator
 *  Wishart random number generator
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Copyright (C) 2006  Ralph dos Santos Silva
 */
int mvnrnd_silva(const gsl_rng *r, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result)
{
	size_t n = mean->size;
	/* multivariate normal distribution random number generator */
	/*
	*	n	dimension of the random vetor
	*	mean    vector of means of size n
	*	var     variance matrix of dimension n x n
	*	result  output variable with a sigle random vector normal distribution generation
	*/

	unsigned int k;
	gsl_matrix *work = gsl_matrix_alloc(n,n);

	gsl_matrix_memcpy(work,var);
/*	printf("I am here\n");*/
	gsl_linalg_cholesky_decomp(work);

	for(k=0; k<n; k++)
        	gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
	gsl_vector_add(result,mean);

	gsl_matrix_free(work);

	return 0;
}


/*
 *  @title multivariate normal random variables
 *  @author Carl Boettiger, <cboettig@gmail.com>
 *
 *  Based on the R function rmvnorm, from the mvtnorm package
 *  by Friedrich Leisch and Fabian Scheipl, implemented
 *  using the GSL libraries
 */

int mvnrnd_cboet(gsl_rng * rng, const gsl_vector * mean, gsl_matrix * covar, gsl_vector * ANS)
{
	unsigned int i;
	size_t n = mean->size;

	/* Calculate eigenvalues and eigenvectors of covar matrix */
	gsl_vector *eval = gsl_vector_alloc (n);
	gsl_matrix *evec = gsl_matrix_alloc (n, n);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
	gsl_eigen_symmv (covar, eval, evec, w);
	gsl_eigen_symmv_free (w);
/*	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);*/



	/* Setup for: evec * matrix(diag(eval)) * transpose(evec)  */
	gsl_matrix *eval_mx = gsl_matrix_calloc (n, n);
	gsl_matrix * x_M = gsl_matrix_alloc (n,n);
	gsl_matrix * x_M_x = gsl_matrix_alloc (n,n);


	gsl_vector_view diagonal = gsl_matrix_diagonal(eval_mx);
	gsl_vector_memcpy(&diagonal.vector, eval);
	for(i=0;i<n;i++)
	{
		gsl_vector_set( &diagonal.vector, 
						i,  
						sqrt( gsl_vector_get(&diagonal.vector, i) )
					  );
	}



	/* evec * matrix(diag(eval)) * transpose(evec)  */
/*	gsl_blas_dsymm (CblasLeft, CblasUpper, 1.0, evec, eval_mx, 0.0, x_M);*/

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 
					1.0, evec, eval_mx, 0.0, x_M);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 
					1.0, x_M, evec, 0.0, x_M_x);


	gsl_matrix_free(x_M);
	gsl_matrix_free(eval_mx);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);

	gsl_vector * rnorms = gsl_vector_alloc(n);
	for(i=0;i<n;i++)
	{ 
		gsl_vector_set 
			( rnorms, i, 
			  gsl_ran_gaussian_ziggurat(rng, 1)
			);
	}

	gsl_blas_dgemv( CblasTrans, 1.0, x_M_x, rnorms, 0, ANS);
	gsl_vector_add(ANS, mean);
	gsl_matrix_free(x_M_x);
	gsl_vector_free(rnorms);

	return 0;
	/* answer provided through pass by reference */
}


int mvnrnd_gsl(gsl_rng * rng, const gsl_vector * mean, gsl_matrix * covar, gsl_vector * ANS)
{
#if 1
	return mvnrnd_silva(rng, mean, covar, ANS);
#else
	return mvnrnd_cboet(rng, mean, covar, ANS);
#endif
}

int mvnrnd(double *mean, double *sigma, double *out, int N)
{
	int res;

	gsl_vector_view mean_view = gsl_vector_view_array(mean, N);
	gsl_matrix_view sigma_view = gsl_matrix_view_array(sigma, N,N);
	gsl_vector_view out_view = gsl_vector_view_array(out, N);

	int me = torc_i_worker_id();
	res = mvnrnd_gsl(r[me], &mean_view.vector, &sigma_view.matrix, &out_view.vector);

	return res;
}


#if 0
void aux_init()
{
	torc_register_task(reset_nfc_task);
	torc_register_task(get_nfc_task);
}
#endif


/* Copyright (C) 2006  Ralph dos Santos Silva */
/* http://lists.gnu.org/archive/html/help-gsl/2006-04/txtdb8Hdlx9uA.txt */

static double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var, int lognorm)
{
	/* multivariate normal density function    */
	/*
	*	n	dimension of the random vector
	*	x	random vector
	*	mean	vector of means of size n
	*	var	variance matrix of dimension n x n
	*/

	int s;
	double ax,ay;
	gsl_vector *ym, *xm;

	gsl_matrix *work = gsl_matrix_alloc(n,n), *winv = gsl_matrix_alloc(n,n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	gsl_matrix_memcpy( work, var );
	gsl_linalg_LU_decomp( work, p, &s );
	gsl_linalg_LU_invert( work, p, winv );
	ax = gsl_linalg_LU_det( work, s );		// ax = det(var);
	gsl_matrix_free( work );
	gsl_permutation_free( p );

	xm = gsl_vector_alloc(n);
	gsl_vector_memcpy( xm, x);
	gsl_vector_sub( xm, mean );
	ym = gsl_vector_alloc(n);
	gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
	gsl_matrix_free( winv );
	gsl_blas_ddot( xm, ym, &ay);
	gsl_vector_free(xm);
	gsl_vector_free(ym);

//	printf("xxx: ax = %e, log(ax) = %e, ay  = %e\n", ax, log(ax), ay);
#if 0
	if (ax == 0.0) ax = DBL_MIN;
#endif

	if (!lognorm)
		ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );
	else
		ay = 0.5*(-ay - n*log(2*M_PI) - log(ax));	/* PA */

	return ay;
}

static void dmvnorm_v(const int m, double *r, const int n, const gsl_vector *x, gsl_vector **mean, const gsl_matrix *var, int lognorm)
{
	/* multivariate normal density function    */
	/*
	*	m	number of evaluation points
	*	r	vector of results
	*	n	dimension of the random vector
	*	x	random vector
	*	mean	m vectors of mean of size n
	*	var	variance matrix of dimension n x n
	*/

	int s;
	double ax,ay;
	gsl_vector *ym, *xm;

	gsl_matrix *work = gsl_matrix_alloc(n,n), *winv = gsl_matrix_alloc(n,n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	gsl_matrix_memcpy( work, var );
	gsl_linalg_LU_decomp( work, p, &s );
	gsl_linalg_LU_invert( work, p, winv );
	ax = gsl_linalg_LU_det( work, s );		// ax = det(var);
	gsl_matrix_free( work );
	gsl_permutation_free( p );

	xm = gsl_vector_alloc(n);
	ym = gsl_vector_alloc(n);

	for (s = 0; s < m; s++)
	{
		gsl_vector_memcpy( xm, x);
		gsl_vector_sub( xm, mean[s] );
		gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
		gsl_blas_ddot( xm, ym, &ay);

		if (!lognorm)
			ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );
		else
			ay = 0.5*(-ay - n*log(2*M_PI) - log(ax));	/* PA */

		r[s] = ay;
	}

	gsl_matrix_free( winv );
	gsl_vector_free(xm);
	gsl_vector_free(ym);


	return;
}

static double dmvnorm2(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var, const gsl_matrix *rmat, double ds, int lognorm)
{
	/* multivariate normal density function    */
	/*
	*	n	dimension of the random vector
	*	x	random vector
	*	mean	vector of means of size n
	*	var	variance matrix of dimension n x n
	*/

	int s;
	double ax,ay;
	gsl_vector *ym, *xm;
	gsl_matrix *work = gsl_matrix_alloc(n,n), *winv = gsl_matrix_alloc(n,n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	gsl_matrix_memcpy( work, rmat );
	gsl_linalg_LU_decomp( work, p, &s );
	ax = gsl_linalg_LU_det( work, s );		// ax = det(var);
	//gsl_matrix_free( work );
	//gsl_permutation_free( p );
	double log_ax = log(ax);

	gsl_matrix_memcpy( work, var );
	gsl_linalg_LU_decomp( work, p, &s );
	gsl_linalg_LU_invert( work, p, winv );
	//double ax2 = gsl_linalg_LU_det( work, s );		// ax = det(var);

	gsl_matrix_free( work );
	gsl_permutation_free( p );

	xm = gsl_vector_alloc(n);
	gsl_vector_memcpy( xm, x);
	gsl_vector_sub( xm, mean );
	ym = gsl_vector_alloc(n);
	gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
	gsl_matrix_free( winv );
	gsl_blas_ddot( xm, ym, &ay);
	gsl_vector_free(xm);
	gsl_vector_free(ym);

#if 0
	if (ax == 0.0) ax = DBL_MIN;
#endif

	if (!lognorm)
		ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );
	else
	//	ay = 0.5*(-ay - n*log(2*M_PI) - log(ax));	/* PA */
		ay = 0.5*(-ay - n*log(2*M_PI) - log_ax - ds);

	return ay;
}

static void gsl_dmvnorm_v(int m, double *r, int n, double *xv, double *mv, double *vm, int lognorm)
{
	int i, j, k;
	
	gsl_vector *x = gsl_vector_calloc(n), **mean; // = gsl_vector_calloc(n);
	gsl_matrix *var = gsl_matrix_calloc(n,n);

	mean = (gsl_vector **) malloc(m*sizeof(gsl_vector *));
	for (k = 0; k < m; k++) mean[k] = gsl_vector_calloc(n);
 

	for (i = 0; i < n; i++)
		gsl_vector_set(x,i,xv[i]);
		

	if (mv != NULL) {
		for (k = 0; k < m; k++)
			for (i = 0; i < n; i++)
				gsl_vector_set(mean[k],i,mv[k*n+i]);
	}

		
	if (vm != NULL) {
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				gsl_matrix_set(var,i,j,vm[i*n+j]);
	}
	else {
		for (i = 0; i < n; i++)
			gsl_matrix_set(var,i,i, 1.0);
	}

	dmvnorm_v(m,r,n,x,mean,var,lognorm);
	
	gsl_vector_free(x);
	for (k = 0; k < m; k++) gsl_vector_free(mean[k]);
	free(mean);
	gsl_matrix_free(var);

	return; 
}

static double gsl_dmvnorm(int n, double *xv, double *mv, double *vm, int lognorm)
{
	double result;
	int i, j;
	
	gsl_vector *x = gsl_vector_calloc(n), *mean = gsl_vector_calloc(n);
	gsl_matrix *var = gsl_matrix_calloc(n,n);

	for (i = 0; i < n; i++)
		gsl_vector_set(x,i,xv[i]);
		
	if (mv != NULL) {
		for (i = 0; i < n; i++)
			gsl_vector_set(mean,i,mv[i]);
	}
		
	if (vm != NULL) {
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				gsl_matrix_set(var,i,j,vm[i*n+j]);
	}
	else {
		for (i = 0; i < n; i++)
			gsl_matrix_set(var,i,i, 1.0);
	}

	result = dmvnorm(n,x,mean,var,lognorm);
	
	gsl_vector_free(x);
	gsl_vector_free(mean);
	gsl_matrix_free(var);

	return result;
}

static double gsl_dmvnorm2(int n, double *xv, double *mv, double *vm, double *r, double ds, int lognorm)
{
	double result;
	int i, j;
	
	gsl_vector *x = gsl_vector_calloc(n), *mean = gsl_vector_calloc(n);
	gsl_matrix *var = gsl_matrix_calloc(n,n);
	gsl_matrix *rmat = gsl_matrix_calloc(n,n);

	for (i = 0; i < n; i++)
		gsl_vector_set(x,i,xv[i]);
		
	if (mv != NULL) {
		for (i = 0; i < n; i++)
			gsl_vector_set(mean,i,mv[i]);
	}
		
	if (vm != NULL) {
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				gsl_matrix_set(var,i,j,vm[i*n+j]);
	}
	else {
		for (i = 0; i < n; i++)
			gsl_matrix_set(var,i,i, 1.0);
	}

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			gsl_matrix_set(rmat,i,j,r[i*n+j]);

	result = dmvnorm2(n,x,mean,var,rmat,ds,lognorm);
	
	gsl_vector_free(x);
	gsl_vector_free(mean);
	gsl_matrix_free(var);
	gsl_matrix_free(rmat);

	return result;
}

double mvnpdf(int n, double *xv, double *mv, double *vm)
{
	double result = gsl_dmvnorm(n, xv, mv, vm, 0);
	return result;
}


void mvnpdf_v(int m, double *r, int n, double *xv, double *mv, double *vm)
{
	gsl_dmvnorm_v(m, r, n, xv, mv, vm, 0);
	return;
}


double logmvnpdf(int n, double *xv, double *mv, double *vm)
{
	double result = gsl_dmvnorm(n, xv, mv, vm, 1);
	return result;
}

double logmvnpdf2(int n, double *xv, double *mv, double *vm, double *rmat, double ds)
{
	double result = gsl_dmvnorm2(n, xv, mv, vm, rmat, ds, 1);
	return result;
}



long randi(int N)
{
	int me = torc_i_worker_id();
	long ri = 1 + gsl_rng_uniform_int(r[me], N-1);

	return ri;
}


void randsample_sorted(int x0, int x1, double *rset, long Neta)
{
	int me = torc_i_worker_id();
 
	//rset = sort( randsample(2:Np,Neta) );

	int Np = x1-x0+1;

	double b[Np];
	for (int i = 0; i < Np; i++)
	{
		b[i] = x0 + (double) i;
	}
	gsl_ran_choose (r[me], rset, Neta, b, Np, sizeof(double));
}

#include <gsl/gsl_linalg.h>


int inverse(double *mat, double *invmat, int N)
{
	int s;

	gsl_matrix_view m   = gsl_matrix_view_array(mat, N, N);
	gsl_matrix_view inv = gsl_matrix_view_array(invmat,N,N);
	gsl_permutation * p = gsl_permutation_alloc (N);

#if 0
	printf("The matrix is\n");
	for (int i = 0; i < N; ++i)
	for (int j = 0; j < N; ++j)
	printf(j==N-1?"%6.3f\n":"%6.3f ", gsl_matrix_get(&m.matrix,i,j));
#endif

	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);

#if 0
	printf("The inverse is\n");
	for (int i = 0; i < N; ++i)
	for (int j = 0; j < N; ++j)
	printf(j==N-1?"%6.3f\n":"%6.3f ",gsl_matrix_get(&inv.matrix,i,j));
#endif

	gsl_permutation_free (p);

	return 0;
}


#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_randist.h>

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

