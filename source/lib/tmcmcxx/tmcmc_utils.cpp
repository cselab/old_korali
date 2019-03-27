#include <cassert>
#include <math.h>
#include <limits>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>

#include "tmcmc_utils.hpp"

namespace tmcmc
{

int in_rect(double *v1, double *v2, double *diam, double sc, int D)
{
    int d;
    for (d = 0; d < D; ++d) {
        if (fabs(v1[d]-v2[d]) > sc*diam[d]) return 0;
    }
    return 1;
}


int compar_desc(const void* p1, const void* p2)
{
    int dir = +1;   /* -1: ascending order, +1: descending order */
    sort_t *s1 = (sort_t *) p1;
    sort_t *s2 = (sort_t *) p2;

    if (s1->nsel < s2->nsel) return dir;
    if (s1->nsel > s2->nsel) return -dir;

    return 0;
}

double compute_sum(double *v, int n)
{
    double s = 0;
    for (int i = 0; i < n; i++) s += v[i];
    return s;
}


double compute_dot_product(double row_vector[], double vector[], int dim)
{
	double sum = 0.0;
	for(int row=0; row<dim; row++) sum += row_vector[row] * vector[row];
  
	return sum;
}


void compute_mat_product_vect(double *mat/*2D*/, double vect[], double res_vect[], double coef, int dim)
{
    int row, column;
    double current_dot_product;

	for(row=0; row<dim; ++row){
		current_dot_product = 0.0;
        for(column=0; column<dim; ++column) current_dot_product += mat[row*dim+column] * vect[column];
        res_vect[row] = coef * current_dot_product;
    }
    return;
}


void inv_matrix(double *current_hessian/*2D*/, double *inv_hessian/*2D*/, int dim)
{
    gsl_matrix_view m   = gsl_matrix_view_array(current_hessian, dim, dim);
    gsl_matrix_view inv = gsl_matrix_view_array(inv_hessian, dim, dim);
    gsl_permutation * p = gsl_permutation_alloc (dim);

    int s;
    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);

    gsl_permutation_free (p);
    return;
}


double scale_to_box(const double* point, double sc, const double* add_vec, const double *elbds, const double *eubds, int dims)
{
	double pp[dims];
	for(int i=0; i<dims; ++i) pp[i] = point[i]+sc*add_vec[i];

	sc = fabs(sc);
	double c;
	for (int l=0; l<dims; l++)
	{
		if (pp[l]<elbds[l])
		{
			c = fabs( (elbds[l]-point[l]) / add_vec[l] );
			sc = fmin(sc,c);
		}
		if (pp[l]>eubds[l])
		{
			c = fabs( (point[l]-eubds[l]) / add_vec[l] );
			sc = fmin(sc,c);
		}
	}
	return sc;
}


void get_nfc_task(int *x)
{
    *x = l_nfeval;
}


int get_nfc()
{
    int c[1024]; /* MAX_NODES*/
#ifdef _USE_TORC_
    for (int i = 0; i < torc_num_nodes(); ++i) {
        torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())get_nfc_task, 1,
                       1, MPI_INT, CALL_BY_RES, &c[i]);
    }
    torc_waitall();
#else
    get_nfc_task(&c[0]);
#endif

    unsigned int s = 0;
#ifdef VERBOSE
    printf("get_nfc:");
#endif
    for (int i = 0; i < torc_num_nodes(); ++i) {
        s += c[i];
#ifdef VERBOSE
        printf("+%d", c[i]);
#endif
    }
    g_nfeval = s;
#ifdef VERBOSE
    printf("=%d\n", s);
#endif
    t_nfeval += g_nfeval;
    return g_nfeval;
}

int get_tfc()
{
    return t_nfeval;
}


void inc_nfc()
{
    pthread_mutex_lock(&feval_m);
    l_nfeval++;
    pthread_mutex_unlock(&feval_m);
}


void reset_nfc_task()
{
    l_nfeval = 0;
}


void reset_nfc()
{
#ifdef _USE_TORC_
    for (int i = 0; i < torc_num_nodes(); ++i) {
        torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())reset_nfc_task, 0);
    }
    torc_waitall();
#else
    reset_nfc_task();
#endif
}


//==================================================//
//=================Force_Pos_Def====================//
//==================================================//
/**
 *
 * input(1): non positive definite matrix
 * input(2): resulted forced positive definite matix (output)
 *
 * remarks: Force matrix to be positive definite following
 *  publication "Efficient stohastic generation of multi-site synthetic precipitation data"
 *  by F.P.Brissette et al.
 * "Modification of a negative eigenvalues to create a positive
 * definite matrices and approximation for standard errors of correlation estimates by L.R. Schaeffer"
 */

void force_pos_def0(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int PROBDIM)
{
    int row, column;
    int exist_zero_eig_value;
    double current_eig_value, sum_neg_eig_values, min_pos_eig_value;
    double normalization_factor, small_pos_eig_value;

    double eps  = 1e-12;	//set a value to place instead of zero or negative eigenvalues
    double zero = 1e-15;	//set a value to look for zero

    gsl_matrix *working_mat, *temp_mat, *eig_vectors, *diag_mat, *intermediate_mat;//, *pos_def_mat;
    gsl_vector *eig_values;
    gsl_eigen_symmv_workspace *work_v;

    working_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
    temp_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
    diag_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
    intermediate_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);

    eig_vectors = gsl_matrix_alloc(PROBDIM, PROBDIM);
    eig_values = gsl_vector_alloc(PROBDIM);

    //copy input non positive matrix twice
    gsl_matrix_memcpy(working_mat,non_pos_def_mat);
    gsl_matrix_memcpy(temp_mat, non_pos_def_mat);

    //Get the eigenystem of the input matrix
    //allocate the space for the eigenvalues and copy temp_mat to gsl matrix
    work_v = gsl_eigen_symmv_alloc(PROBDIM);

    gsl_eigen_symmv(temp_mat, eig_values, eig_vectors, work_v); //get eigenvalues and eigenvectors

    min_pos_eig_value = std::numeric_limits<double>::max();
    //	exist_pos_eig_values = 0; //set the flag for existence of positive eigenvalues to FALSE
    sum_neg_eig_values = 0.0;
    exist_zero_eig_value = 0; //set the flag for zero eigenvalue existence to FALSE
    //first_neg_value = 1;

    /**
     * loop through the eigenvalues to find:
     * the minimum positive value
     * the sum of negative values
     */

    if(gsl_vector_isneg(eig_values)) { //if all eigenvalues are negative
        min_pos_eig_value = -gsl_vector_max(eig_values); //if all values are negative assumed that the min positive value is the maximum of the negative values
        for(column=0; column<PROBDIM; column++) {	//sum of negative values equals the sum of all vector elements
            sum_neg_eig_values += gsl_vector_get(eig_values, column);
        }
        //exist_pos_eig_values = 1;
    } else {
        for(column=0; column < PROBDIM; column++) {
            current_eig_value = gsl_vector_get(eig_values, column);
            if(current_eig_value > zero) { // if the current eigenvalue is positive, try to update the minimum positive eigenvalue
                if(current_eig_value <= min_pos_eig_value) {
                    min_pos_eig_value = current_eig_value;
                }
                //exist_pos_eig_values = 1; // positive eigenvalue was found, set the flag to TRUE
            } else if (current_eig_value < 0.0) {
                sum_neg_eig_values += current_eig_value;
            } else {
                //handling zero eigenvalue
                printf("Warning: Zero eigenvalue\n");
                exist_zero_eig_value = 1; //as you find zero eigenvalue, set flag to true
            }
        }
    }

    if (sum_neg_eig_values == 0.0 && exist_zero_eig_value) { //let's check if only zero eigenvalues hinder matrix to be positive
        sum_neg_eig_values = 0.01; //set this value for correcting zero eigenvalues

        if (min_pos_eig_value == 1000000.0) {
            printf("Warning: Correction of all zero eigenvalues \n");
            min_pos_eig_value = 101.0;
            sum_neg_eig_values = 1.0;
            //exist_pos_eig_values = 1; //modifying flag as we recover the situtation of all zero eigenvalues
        }
    }

    assert(sum_neg_eig_values != 0.0); //test that the sum of negative values is negative(!)
    //	assert(exist_pos_eig_values && min_pos_eig_value != 1000000.0); //assert that the starting value for minimum is enough big

    sum_neg_eig_values = 2.0 * sum_neg_eig_values;
    //follow the publication on remarks to force the matrix to be positive definite
    normalization_factor = ( pow((sum_neg_eig_values),2) * 100.0) + 1.0; //compute the squared sum of negative values 2.0 * sum_neg
    //update the negative eigenvalues to create new small positive ones
    for(column=0; column< PROBDIM; column++) {
        current_eig_value = gsl_vector_get(eig_values, column);

        if (current_eig_value <= zero) { //<= for zero or negative eigenvalue
            small_pos_eig_value = min_pos_eig_value * ( (pow(sum_neg_eig_values - current_eig_value,2)) / normalization_factor );

            while (small_pos_eig_value < eps) { //iterative find the closeste value to given eps
                small_pos_eig_value *= 10.0;
            }
            gsl_vector_set(eig_values, column, small_pos_eig_value); //update the eigenvalue
        }
    }

    //follows a Eigen_Vectors * diag(Eigen_Values) * (Eigen_Vectors)^T
    gsl_matrix_set_identity(diag_mat); //construct a diagonal matrix with eigenvalues in the first diagonal
    //hard copying the eigenvalues to the identity matrix to construct a diag[eig1 eig2 .. eign]
    for(row=0; row<PROBDIM; row++) {
        gsl_matrix_set(diag_mat, row, row, gsl_vector_get(eig_values, row));
    }

    //now compute intermediate matrix = [v1 v2 .. vn] * diag[eig1 eig2 .. eign]
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eig_vectors, diag_mat, 0.0, intermediate_mat);

    //finally compute pos definite matrix = [v1 v2 .. vn]  * diag[eig1 eig2 .. eign] * [v1 v2 .. vn]^T
    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, intermediate_mat, eig_vectors, 0.0, forced_pos_def_mat); //pos_def_mat

    //free the workspace
    gsl_eigen_symmv_free(work_v);
    //free the matrices
    gsl_matrix_free(working_mat);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(diag_mat);
    gsl_matrix_free(intermediate_mat);
    gsl_matrix_free(eig_vectors);
    //free the vector
    gsl_vector_free(eig_values);
}


void force_pos_def1(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int PROBDIM)
{
    //	[v,d] = eig(a);
    //	a_psd = v * diag(max(diag(d), eps))/v;

    int row;
    //	double eps = 1e-3;  //set a value to place instead of zero or negative eigenvalues
    double eps = 1e-1;  //set a value to place instead of zero or negative eigenvalues

    gsl_matrix *working_mat, *temp_mat, *eig_vectors, *diag_mat, *intermediate_mat;
    gsl_vector *eig_values;
    gsl_eigen_symmv_workspace *work_v;

    working_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
    temp_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
    diag_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);
    intermediate_mat = gsl_matrix_alloc(PROBDIM, PROBDIM);

    eig_vectors = gsl_matrix_alloc(PROBDIM, PROBDIM);
    eig_values = gsl_vector_alloc(PROBDIM);

    //copy input non positive matrix twice
    gsl_matrix_memcpy(working_mat,non_pos_def_mat);
    gsl_matrix_memcpy(temp_mat, non_pos_def_mat);

    //Get the eigenystem of the input matrix
    //allocate the space for the eigenvalues and copy temp_mat to gsl matrix
    work_v = gsl_eigen_symmv_alloc(PROBDIM);

    gsl_eigen_symmv(temp_mat, eig_values, eig_vectors, work_v); //get eigenvalues and eigenvectors

    //follows a Eigen_Vectors * diag(Eigen_Values) * (Eigen_Vectors)^T
    gsl_matrix_set_identity(diag_mat); //construct a diagonal matrix with eigenvalues in the first diagonal
    //hard copying the eigenvalues to the identity matrix to construct a diag[eig1 eig2 .. eign]

    for(row=0; row<PROBDIM; row++) {
        double eig_val = gsl_vector_get(eig_values, row);
        if (eig_val < eps) eig_val = eps;
        gsl_matrix_set(diag_mat, row, row, eig_val);
    }

    //now compute intermediate matrix = [v1 v2 .. vn] * diag[eig1 eig2 .. eign]
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eig_vectors, diag_mat, 0.0, intermediate_mat);

    //finally compute pos definite matrix = [v1 v2 .. vn]  * diag[eig1 eig2 .. eign] * [v1 v2 .. vn]^T
    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, intermediate_mat, eig_vectors, 0.0, forced_pos_def_mat); //pos_def_mat

    //free the workspace
    gsl_eigen_symmv_free(work_v);
    //free the matrices and vector
    gsl_matrix_free(working_mat);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(diag_mat);
    gsl_matrix_free(intermediate_mat);
    gsl_matrix_free(eig_vectors);
    gsl_vector_free(eig_values);
}

int check_mat_pos_def(gsl_matrix *mat_to_check, int PROBDIM);	// forward declaration

//  Copyright 2003, Pekka Paalanen <pekka.paalanen@lut.fi>
void force_pos_def2(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int PROBDIM)
{
    //	D = size(sigma, 1);
    //	fixrate = 0.01;
    //	covfixmat = ones(D) + fixrate*eye(D);
    //	loops = 0;
    //	min_limit = eps*10;

    int i; //, row;
    double fixrate = 0.01;
    int loops = 0;
    double min_limit = 1e-15;
    double convfix = 1 + fixrate;

    gsl_matrix *nsigma;

    nsigma = gsl_matrix_alloc(PROBDIM, PROBDIM);
    gsl_matrix_memcpy(nsigma, non_pos_def_mat);

    double d[PROBDIM];
    //	gsl_matrix_fprintf (stdout, nsigma, "%10.3f");

    while (check_mat_pos_def(nsigma, PROBDIM) == 0) {
        loops++;
        printf("loops = %d\n", loops);

        int below_limit = 0;
        for (i = 0; i < PROBDIM; i++) {
            d[i] = gsl_matrix_get(nsigma, i, i);
            if (d[i] <= min_limit) below_limit = 1;
        }

        double maxE_abs = fabs(d[0]);
        double minE = d[0];
        for (i = 1; i < PROBDIM; i++) {
            if (fabs(d[i]) > maxE_abs) maxE_abs = fabs(d[i]);
            if (d[i] < minE) minE = d[i];
        }

        double m = maxE_abs * fixrate;
        double neg = minE;
        double addit;

        //		printf("m = %f, neg = %f\n", m, neg);
        if (below_limit) {
            if (neg < 0) {
                addit = (m - neg);
            } else {
                if (m < min_limit) m = min_limit;
                addit = m;
            }
            for (i = 0; i < PROBDIM; i++) {
                double newv = gsl_matrix_get(nsigma, i, i);
                newv = newv + addit;
                gsl_matrix_set(nsigma, i, i, newv);
            }
        } else {
            for (i = 0; i < PROBDIM; i++) {
                double newv = gsl_matrix_get(nsigma, i, i);
                newv = newv * convfix;
                gsl_matrix_set(nsigma, i, i, newv);
            }
        }
    }

    gsl_matrix_memcpy(forced_pos_def_mat, nsigma);
    gsl_matrix_free(nsigma);
}


void force_pos_def3(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int PROBDIM)
{
    int row, column;
    for (row = 0; row < PROBDIM ; row++) {
        for (column = 0; column < PROBDIM; column++)
            gsl_matrix_set(forced_pos_def_mat, row, column, 0.0);
    }

    for (row = 0; row < PROBDIM ; row++)
        gsl_matrix_set(forced_pos_def_mat, row, row, 1.0);
}


void force_pos_def(gsl_matrix *non_pos_def_mat, gsl_matrix *forced_pos_def_mat, int method, int PROBDIM)
{
    switch (method) {
    case  0:
        force_pos_def0(non_pos_def_mat, forced_pos_def_mat, PROBDIM);
        break;
    case  1:
        force_pos_def1(non_pos_def_mat, forced_pos_def_mat, PROBDIM);
        break;
    case  2:
        force_pos_def2(non_pos_def_mat, forced_pos_def_mat, PROBDIM);
        break;
    case  3:
    default:
        force_pos_def3(non_pos_def_mat, forced_pos_def_mat, PROBDIM);
        break;
    }
}


int check_mat_pos_def(gsl_matrix *mat_to_check, int PROBDIM)
{
    gsl_set_error_handler_off();

    gsl_eigen_symm_workspace *work;
    gsl_vector *eig_values;
    int is_pos_def;

    //Get the eigenvalues of hessian and check if there are positive
    work = gsl_eigen_symm_alloc(PROBDIM);
    gsl_matrix *mat_copy = gsl_matrix_alloc(PROBDIM, PROBDIM);
    gsl_matrix_memcpy(mat_copy, mat_to_check);

    eig_values = gsl_vector_alloc(PROBDIM);
    gsl_eigen_symm(mat_copy, eig_values, work);
    is_pos_def = gsl_vector_ispos(eig_values);

    gsl_vector_free(eig_values);
    gsl_eigen_symm_free(work);
    gsl_matrix_free(mat_copy);

    return is_pos_def;
}


int make_posdef(double *mat, int dim, int method)
{
    int row, column;
    int res = 1;

    gsl_matrix *hessian_mat = gsl_matrix_alloc(dim, dim);
    gsl_matrix *pos_def_mat = gsl_matrix_alloc(dim, dim);
    gsl_matrix *hessian_mat_cpy = gsl_matrix_alloc(dim, dim);

    for(row=0; row<dim; row++) {
        for(column=0; column<dim; column++) {
            gsl_matrix_set(hessian_mat, row, column, mat[row*dim+column]);
        }
    }
    gsl_set_error_handler_off();

    gsl_matrix_memcpy(pos_def_mat, hessian_mat);    // create a copy
    gsl_matrix_memcpy(hessian_mat_cpy, hessian_mat);// create a copy

    int is_pos_def = check_mat_pos_def(hessian_mat_cpy, dim); //check if the matrix is positive definite

    if(!is_pos_def) {       // if not, apply the fix, store	the result in pos_def_mat
        force_pos_def(hessian_mat, pos_def_mat, method, dim);
        res = 1;
    } else {
        res = 0;
    }

    for(row=0; row<dim; row++) {
        for(column=0; column<dim; column++) {
            mat[row*dim+column] = gsl_matrix_get(pos_def_mat, row, column);
        }
    }

    gsl_matrix_free(hessian_mat);
    gsl_matrix_free(pos_def_mat);
    gsl_matrix_free(hessian_mat_cpy);

    return res;

}

void print_matrixi(const char *name, int *x, int n)
{
    printf("\n%s =\n\n", name);
    for (int i = 0; i < n; ++i) printf("   %20.15d\n", x[i]);
    printf("\n");
}


void print_matrix(const char *name, double *x, int n)
{
    printf("\n%s =\n\n", name);
    for (int i = 0; i < n; ++i) printf("   %20.15lf\n", x[i]);
    printf("\n");
}

void print_matrix(const char *name, double *x, int n1, int n2)
{
    printf("\n%s =\n\n", name);
    for (int i = 0; i < n1; ++i) {
        for(int j = 0; j < n2; ++j) printf("   %20.15lf", x[i*n1+j]);
        printf("\n");
    }
    printf("\n");
}


void print_matrix_i(char *title, int *v, int n)
{
    printf("\n%s =\n\n", title);
    for (int i = 0; i < n; i++) printf("  %8d\n", v[i]);
    printf("\n");
}


void print_matrix_2d(const char *name, double **x, int n1, int n2)
{
    printf("\n%s =\n\n", name);
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            printf("   %20.15lf", x[i][j]);
        }
        printf("\n");
    }
    printf("\n");

}


} //namespace
