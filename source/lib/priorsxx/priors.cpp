#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "density.hpp"
#include "priors.hpp"

namespace priors {

	Prior::Prior(const char * fname) {

		density_factory_func(fname, densities_, &dim_);
	}

	Prior::~Prior() {
		delete[] densities_;
	}


	double Prior::eval_pdf(double *x){
		double res=1.;
		for(int i=0; i<dim_; ++i) res *= densities_[i].eval(x[i]);
		return res;
	}


	double Prior::eval_logpdf(double *x){
		double res=1.;
		for(int i=0; i<dim_; ++i) res *= densities_[i].log_eval(x[i]);
		return res;
	}

	void Prior::print() {
		printf("==============================\n");
		printf("===  Prior Distribution    ===\n");
		printf("==============================\n");

		for(int i=0; i<dim_; ++i) densities_[i].print();

		printf("==============================\n");
	}

}//namspace priors
