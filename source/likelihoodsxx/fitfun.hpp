#ifndef FITFUN_HPP
#define FITFUN_HPP

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
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

typedef double (*model_ptr) (const double* x, int n);
typedef gsl_vector* (*grad_model_ptr) (const double* x, int n);

class Fitfun : public IFitfun
{

public:
	Fitfun(model_ptr model, grad_model_ptr grad = nullptr) : _model(model), _grad(grad) {};

    double evaluate (const double* x, int n, void* output, int* info);
    
    void initialize(int argc, const  char **argv) {};

    void finalize() {};

private:
 	model_ptr _model;
    
    grad_model_ptr _grad;

};

inline double Fitfun::evaluate (const double* x, int n, void* output, int* info) { 
    if (_grad == nullptr) return _model(x, n);
    double llk = _model(x,n);
    
    return_type* result = static_cast<return_type*>(output);
    result->loglike   = llk;
    result->error_flg = 2; //TODO: check (DW)
    result->posdef    = 1;
    result->grad      = _grad(x,n);
    //TODO: finish output assignments (DW)
    return llk;

};

}//namespace fitfun

#endif//FITFUN_HPP

