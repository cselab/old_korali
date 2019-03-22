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
    int posdef;
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

    double evaluate (const double* x, int n, void* output, int* info);
    
    void initialize(int argc, const  char **argv) {};

    void finalize() {};

private:
 	fmodel_ptr _model;
    
    grad_model_ptr _grad;
    
    hess_model_ptr _hess;

};

inline double Fitfun::evaluate (const double* x, int n, void* output, int* info) { 
    if (_grad == nullptr || _hess == nullptr) return _model(x, n); // TODO: (DW) find a way if only one of it is availale
    double llk = _model(x,n);
    
    return_type* result = static_cast<return_type*>(output);
    result->loglike   = llk;
    result->grad      = _grad(x,n);
   
    gsl_matrix* hess  = _hess(x,n);

    gsl_permutation * permutation = gsl_permutation_alloc(n);
    gsl_matrix * inv_hess = gsl_matrix_alloc(n,n);
   
    int LU_dec_err = 0, signum;
    LU_dec_err = gsl_linalg_LU_decomp(hess, permutation, &signum);
    
    int LU_inv_err = 0;
    if(LU_dec_err == 0)
    	LU_inv_err = gsl_linalg_LU_invert(hess, permutation, inv_hess);
   
    if (LU_dec_err != 0 || LU_inv_err != 0) {
        if(LU_dec_err != 0) printf("Fitfun::evaluate : Error in LU decomp. \n");
	    else printf("Fitfun::evaluate : Error in LU invert. \n");
        
        result->error_flg = 2;
        result->posdef    = 0;
        gsl_matrix_free(hess);
        gsl_matrix_free(inv_hess);
        gsl_permutation_free(permutation);
        return llk;
    }

    gsl_matrix * inv_hess_work = gsl_matrix_alloc(n,n);
    gsl_matrix_memcpy (inv_hess_work, inv_hess);

    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
    gsl_eigen_symmv(inv_hess_work, eval, evec, w);

    gsl_eigen_symmv_free (w);
    gsl_matrix_free(inv_hess_work);
    gsl_matrix_free(hess);
    gsl_permutation_free(permutation);

    bool posdef = (gsl_vector_min(eval) > 0.0);
    int sigma_err = 0;
    for(int i = 0; i<n; ++i) {
        if (!isfinite(gsl_vector_get(eval,i))) {
            sigma_err = 1;
            printf("Fitfun::evaluate : Error in inv_hess (not posdef). \n");
	        break;
        }
    }

    for(int i = 0; i<n; ++i)
        for(int j = 0; j<n; ++j) {
            if (!isfinite(gsl_matrix_get(evec,i,j))) {
                printf("Fitfun::evaluate : Error in inv_hess (evec not finite).\n"); 
            	sigma_err = 1;
		        break;
            }
            if (!isfinite(gsl_matrix_get(inv_hess,i,j))) {
                printf("Fitfun::evaluate : Error in inv_hess (inv_hess not finite).\n"); 
                sigma_err = 1;
	    	    break;
	        }
        }

    if (sigma_err != 0) {
        result->error_flg = 2;
        result->posdef    = 0;
        gsl_matrix_free(inv_hess);
        gsl_permutation_free(permutation);
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        return llk;
    }
   
    result->error_flg = false;
    result->eval      = eval;
    result->cov       = inv_hess;
    result->evec      = evec;
    result->posdef    = posdef;
    
    return llk;

};

}//namespace fitfun

#endif//FITFUN_HPP

