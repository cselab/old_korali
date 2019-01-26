#ifndef DENSITY_HPP
#define DENSITY_HPP

namespace priors {

	class Density {

		public:
            
			void init(const char * name, double (*density) (double, double*),  
									   double (*logdensity) (double, double*), 
									   double (*ran) (double*), double * param, int n);

			~Density();
			double eval(double x);
			double log_eval(double x);
			double rand();
			void  print();

		private:
			char name_[32];			 		 // Name of the distribution	
			int npar_;						 // Number of parameters
			double *par_;					 // Parameters of density

			double (*f_)(double, double *);	 // Density function
			double (*lf_)(double, double *);	 // Log of density
			double (*r_)(double *);			 // Random number

	};

	void density_factory_func(const char *file, Density **out_densities, int &out_dim );
	void check_n(int N);

}//namespace priors

#endif//DENSITY_HPP
