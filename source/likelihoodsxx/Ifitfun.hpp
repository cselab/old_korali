#ifndef IFITFUN_HPP
#define IFITFUN_HPP

namespace fitfun{


class IFitfun
{

public:
 
    virtual double evaluate (const double * x, size_t n, void* output,int * info) = 0;
    
    virtual double evaluateM (const double * x, size_t n, void* output,int * info) = 0;
    
    virtual void initialize(int argc, const  char **argv) = 0;

    virtual void finalize() = 0;

};

}//namespace fitfun

#endif//IFITFUN_HPP
