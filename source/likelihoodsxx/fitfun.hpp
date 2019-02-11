#ifndef FITFUN_HPP
#define FITFUN_HPP

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "Ifitfun.hpp"

namespace fitfun{

typedef struct {
    double loglike;
    int error_flg;
    int posdef;
    gsl_vector* grad;
    gsl_matrix* cov;
    gsl_matrix* evec;
    gsl_vector* eval;
} return_type;

typedef double (*model_ptr) (const double* x, int n);

class Fitfun : public IFitfun
{

public:
	Fitfun(model_ptr model) : _model(model) {};

    double evaluate (const double* x, int n, void* output, int* info) { return _model(x, n); };

    void initialize(int argc, const  char **argv) {};

    void finalize() {};

private:
 	model_ptr _model;

};

}//namespace fitfun

#endif//FITFUN_HPP

