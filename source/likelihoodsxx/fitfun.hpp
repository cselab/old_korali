#ifndef FITFUN_HPP
#define FITFUN_HPP

namespace fitfun {

	typedef double (*fitfun_t)(const double*, int);

    void fitfun_initialize(int argc, const  char **argv);
	
	double fitfun(const double *x, int N, void *output, int *info);
    
    void fitfun_finalize();

}

#endif //FITFUN_HPP
