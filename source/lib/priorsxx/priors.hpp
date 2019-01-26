#ifndef PRIORS_HPP
#define PRIORS_HPP

#include "density.hpp"

namespace priors {

	class Prior {

	public:
		Prior(const char * fname);
		~Prior();

		double eval_pdf(double *x);
		double eval_logpdf(double *x);
        double rand(int dim);

		void  print();

	private:
		int dim_;
		Density * densities_;

		// (TODO: do we need this below? (DW))
        //void reassign_prior( Density *p, int Np, double *psi );
		//void new_prior_from_prior( Density **new_prior, Density *from_prior, int Npr );

	};

}//namespace priors

#endif//PRIORS_HPP
