#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iterator>

#include "myrand.hpp"
#include "density.hpp"

#define BLANKS " \t\n"

namespace priors
{

void Density::init(const char * name, double (*density) (double, double*),
                   double (*logdensity) (double, double*),
                   double (*ran) (double*),
                   double * param, int n)
{

    strcpy(name_, name);

    f_  = density;
    lf_ = logdensity;
    r_  = ran;

    npar_ = n;
    par_  = new double[n];
    for(int i = 0; i < n; ++i) par_[i] = param[i];

}

Density::~Density()
{
    delete[] par_;
}

double Density::eval(double x)
{
    return f_(x,par_);
}

double Density::log_eval(double x)
{
    return lf_(x,par_);
}

double Density::rand()
{
    return r_(par_);
}

void Density::print()
{
    if( strcmp(name_,"uniform")==0 )
        printf("%s: %lf  -  %lf \n", name_, par_[0], par_[1]);
    else if( strcmp(name_,"normal")==0 )
        printf("%s: %lf  -  %lf \n", name_, par_[0], par_[1]);
    else if( strcmp(name_,"exponential")==0 )
        printf("%s: %lf \n", name_, par_[0]);
    else if( strcmp(name_,"gamma")==0 )
        printf("%s: %lf  -  %lf \n", name_, par_[0], par_[1]);
    else
        printf("Did not recognize density: %s\n", name_);
}

void density_factory_func(const char *file, Density **out_densities, int &out_dim )
{

    FILE *fp = fopen( file,"r");
    if(fp==nullptr) {
        printf("\n%s could not be opened. Exit...\n", file );
        exit(EXIT_FAILURE);
    }

    int N=-1;

    size_t bufsize = 1024;
    char * buffer  = new char[bufsize];

    Density * densities = nullptr;
    int cnt = 0;

    while(  -1 != getline( &buffer, &bufsize, fp) ) {

        char * pch = strtok (buffer, BLANKS );
        while( pch != nullptr ) {

            if( pch[0]=='#' || pch[0]=='\n' )
                break;

            if( strcmp(pch,"N")==0 ) {
                pch = strtok (nullptr, BLANKS );
                N   = atoi(pch);
                check_n(N);
                densities = new Density[N];
                break;
            }

            if( strcmp(pch,"uniform")==0 || strcmp(pch,"uni")==0 ) {
                double * par_ = new double[2];
                par_[0] = atof ( strtok (nullptr, BLANKS) );
                par_[1] = atof ( strtok (nullptr, BLANKS) );
                densities[cnt].init("uniform", &uniform_pdf, &uniform_log_pdf, &uniform_rnd, par_, 2);
                cnt++;
                break;
            }


            if( strcmp(pch,"normal")==0 || strcmp(pch,"gaussian")==0 ) {
                double * par_ = new double[2];
                par_[0] = atof ( strtok (nullptr, BLANKS) );
                par_[1] = atof ( strtok (nullptr, BLANKS) );
                densities[cnt].init("normal", &normal_pdf, &normal_log_pdf, &normal_rnd, par_, 2);
                cnt++;
                break;
            }

            if( strcmp(pch, "exp")==0 || strcmp(pch,"exponential")==0  ) {
                double * par_ = new double[1];
                par_[0] = atof ( strtok (nullptr, BLANKS) );
                densities[cnt].init("exponential", &exp_pdf, &exp_log_pdf, &exp_rnd, par_, 1);
                cnt++;
                break;
            }

            if( strcmp(pch, "gam")==0 || strcmp(pch,"gamma")==0  ) {
                double * par_ = new double[2];
                par_[0] = atof ( strtok (nullptr, BLANKS) );
                par_[1] = atof ( strtok (nullptr, BLANKS) );
                densities[cnt].init("gamma", &gamma_pdf, &gamma_log_pdf, &gamma_rnd, par_, 2);
                cnt++;
                break;
            }

            puts(buffer);
            printf("\nSomething went wrong while reading the prior parameter file. Exit...\n");
            exit(EXIT_FAILURE);

        }
    }

    if ( N != cnt ) {
        printf("\nNumber of priors does not match prior definition count. Exit...\n");
        exit(EXIT_FAILURE);
    }

    fclose(fp);
    delete[] buffer;
    *out_densities = densities;
    out_dim = N;
}


void check_n( int N )
{
    if(N<=0) {
        puts("Check the priors parameter file. Error with reading N. Exit...\n");
        exit(EXIT_FAILURE);
    }
}


}//namespace priors
