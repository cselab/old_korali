#include <fstream>
#include <stdio.h>
#include <functional>

#include <libgp/libgp/include/core/gp.h>
#include <libgp/libgp/include/core/gp_utils.h>
#include <libgp/libgp/include/core/rprop.h>

#include "blood.hpp"
#include "fitfun.hpp"
#include "engine_cmaes.hpp"

using namespace libgp;
using namespace fitfun;

int P = 6;

int M = 1000000;

double bds_q1 [] = {0, 1};
double bds_q2 [] = {0, 1};
double bds_q3 [] = {0, 1};
double bds_q4 [] = {0, 1};
double bds_mu [] = {1, 5};
double bds_sig [] = {0.0, 0.02};

double *bds;

void write_best(int idx, CmaesEngine& engine);
void testGp(GaussianProcess& gp, double low, double up, size_t N);

const char* fname = "amresults6.txt";

int main(int argc, char** argv)
{
    switch(P) {
        case 1 : bds = bds_q1; break;
        case 2 : bds = bds_q2; break;
        case 3 : bds = bds_q3; break;
        case 4 : bds = bds_q4; break;
        case 5 : bds = bds_mu; break;
        case 6 : bds = bds_sig; break;
    }
    
    FILE * pFile = fopen(fname, "a");
    if (pFile == NULL) { printf("ERROR: could not open file!!!"); return 1; }
    fprintf(pFile, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    fprintf(pFile, " Profiling Parameter %d with bounds (%f, %f) and %d Steps\n", P, bds[0], bds[1], M);
    fprintf(pFile, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    fclose (pFile); 

    for(int i = 0; i <= M; ++i) {
        
        double dx = (double)i*(bds[1]-bds[0])/(double)M + bds[0];
        printf("\n\nITERATION %d of %d (dx = %.9f)\n", i, M, dx);
        auto f = [dx] (double* theta, int N, void*, int*) { double x[N+1]; for(int j = 0; j<N+1; ++j) { if ( j < (P-1) ) x[j] = theta[j]; else if (j == (P-1)) x[j] = dx; else x[j] = theta[j-1]; }  return -gpllk(x, N+1); };
        auto engine = CmaesEngine(f, "./arun/");
        engine.run();
        write_best(i , engine);
    }

}


void write_best(int idx, CmaesEngine& engine)
{
    FILE * pFile = fopen(fname, "a");
    if (pFile == NULL) { printf("ERROR: could not open file!!!"); return; }
    
    auto evo = engine.getEvo();
    double bestFunVal = cmaes_Get(evo,"fbestever");

    double dx = (double)idx*(bds[1]-bds[0])/(double)M + bds[0];;

    fprintf(pFile, "%d\t", idx);
    double x;
    for(int i = 0; i < evo->sp.N+1; ++i) {
        if ( i < (P-1) )     x = evo->rgxbestever[i]; 
        else if (i == (P-1)) x = dx; 
        else                 x = evo->rgxbestever[i-1];  
        
        fprintf(pFile, "%.15f\t", x);
    }
    fprintf(pFile, "%.15f\n", bestFunVal);

    fclose (pFile); 
    return;
}
