#ifndef FITFUN_HPP
#define FITFUN_HPP

#include <functional>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "Ifitfun.hpp"

namespace fitfun{

typedef struct {
    double loglike;
    int error_flg; //0: all good, 1: only llk, 2: llk & grad ok
    bool posdef;
    gsl_vector* grad;
    gsl_matrix* cov;
    gsl_matrix* evec;
    gsl_vector* eval;
} return_type;

using fmodel_ptr = std::function<double(const double *x, int n)>;
using grad_model_ptr = std::function<gsl_vector*(const double* x, int n)>;
using hess_model_ptr = std::function<gsl_matrix*(const double* x, int n)>;

class Fitfun : public IFitfun
{

public:
    
	Fitfun(fmodel_ptr model, grad_model_ptr grad = nullptr, hess_model_ptr hess = nullptr) : _model(model), _grad(grad), _hess(hess) {};
    
    // returns loglikelihood, output untouched
    double evaluate (const double* x, size_t n, void* output, int* info);
    
    // returns loglikelihood, output populated with grad & cov of loglikelihood
    double evaluateM (const double* x, size_t n, void* output, int* info);
    
    void initialize(int argc, const  char **argv) {};

    void finalize() {};

private:
 	fmodel_ptr _model;
    
    grad_model_ptr _grad;
    
    hess_model_ptr _hess;

};

inline double Fitfun::evaluate (const double* x, size_t n, void* output, int* info) { return _model(x,n); };

inline double Fitfun::evaluateM (const double* x, size_t n, void* output, int* info) { 
    
    double llk = _model(x,n);

    if (_grad == nullptr || _hess == nullptr) {
        printf("WARNING: Fitfun::evaluateM() _grad and or _hess not defined!! returning llk..\n");
        return llk;
    }
    
    return_type* result = static_cast<return_type*>(output);
    result->loglike   = llk;

    result->grad      = _grad(x,n);

    /*
    gsl_vector* invgrad = gsl_vector_calloc(n);
    gsl_vector_memcpy(invgrad, grad);

    gsl_matrix* igrad2 = gsl_matrix_calloc(n,n);
    for(size_t i = 0; i<n; ++i)
        for(size_t j = 0; j<=i; ++j) {
            double c = -gsl_vector_get(invgrad,i)*gsl_vector_get(grad,j);
            gsl_matrix_set(igrad2,i,j,c);
            gsl_matrix_set(igrad2,j,i,c);
        }

    gsl_vector_free(invgrad);

    gsl_matrix* hess  = _hess(x,n);
    gsl_matrix_scale(hess, 1.0/llk);

    gsl_matrix_add(hess, igrad2);
    gsl_matrix_free(igrad2);
    */


    gsl_permutation * permutation = gsl_permutation_alloc(n);
    gsl_matrix * inv_neg_hess = gsl_matrix_alloc(n,n);

    gsl_matrix* hess  = _hess(x,n);

    int LU_dec_err = 0, signum;
    LU_dec_err = gsl_linalg_LU_decomp(hess, permutation, &signum);
    
    int LU_inv_err = 0;
    if(LU_dec_err == 0)
    	LU_inv_err = gsl_linalg_LU_invert(hess, permutation, inv_neg_hess);
   
    gsl_permutation_free(permutation);
    
    if (LU_dec_err != 0 || LU_inv_err != 0) {
        if(LU_dec_err != 0) printf("Fitfun::evaluate : Error in LU decomp. \n");
	    else printf("Fitfun::evaluate : Error in LU invert. \n");
        result->error_flg = 2;
        result->posdef    = false;
        gsl_matrix_free(hess);
        gsl_matrix_free(inv_neg_hess);
        return llk;
    }

    gsl_matrix_scale(inv_neg_hess, -1.0);

    gsl_matrix * inv_hess_work = gsl_matrix_alloc(n,n);
    gsl_matrix_memcpy (inv_hess_work, inv_neg_hess);

    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
    gsl_eigen_symmv(inv_hess_work, eval, evec, w);

    gsl_eigen_symmv_free (w);
    gsl_matrix_free(inv_hess_work);
    gsl_matrix_free(hess);

    bool posdef = (gsl_vector_min(eval) > 0.0);
    int sigma_err = 0;
    for(size_t i = 0; i<n; ++i) {
        if (!isfinite(gsl_vector_get(eval,i))) {
            sigma_err = 1;
            printf("Fitfun::evaluate : Error in inv_neg_hess (not posdef). \n");
	        break;
        }
    }

    for(size_t i = 0; i<n; ++i)
        for(size_t j = 0; j<n; ++j) {
            if (!isfinite(gsl_matrix_get(evec,i,j))) {
                printf("Fitfun::evaluate : Error in evec (evec not finite).\n"); 
            	sigma_err = 1;
		        break;
            }
            if (!isfinite(gsl_matrix_get(inv_neg_hess,i,j))) {
                printf("Fitfun::evaluate : Error in inv_neg_hess (inv_hess not finite).\n"); 
                sigma_err = 1;
	    	    break;
	        }
        }

    if (sigma_err != 0) {
        result->error_flg = 2;
        result->posdef    = false;
        gsl_matrix_free(inv_neg_hess);
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        return llk;
    }
   
    result->error_flg = false;
    result->posdef    = posdef;
    result->cov       = inv_neg_hess;
    result->eval      = eval;
    result->evec      = evec;
    
    return llk;

};

}//namespace fitfun

#endif//FITFUN_HPP

