#ifndef ERROR_HELPERS_HPP
#define ERROR_HELPERS_HPP

void print_err1(double* mean, double **SS, double LB, double UB, int N)
{

    double err = 0;
    for(int i = 0; i < N; ++i) {
        err += fabs(mean[i])/N;
        for(int j = 0; j < N; ++j) {
            if (i == j) err += fabs(1 - SS[i][j])/(N*N);
            else err += fabs(SS[i][j])/(N*N);
        }
    }

    delete mean;
    for(int i = 0; i < 5; ++i) delete SS[i];
    delete SS;

    double explogevidence = -0.5*N*log(2.0*M_PI) - N*log(UB-LB);
    printf("err: %f\n",err);
    printf("expected log evidence: %f\n", explogevidence);
 
}

void print_err_le(double *le, int N, double LB, double UB, int NRUNS)
{
    double explogevidence = -0.5*N*log(2.0*M_PI) - N*log(UB-LB);
    printf("expected log evidence: %f\n", explogevidence);
    
    double l2err = 0;
    for(int i = 0; i < NRUNS; ++i) {
        printf("%f\n",le[i]);
        l2err += pow((le[i]-explogevidence),2);
    }
    l2err /= NRUNS;

    printf("L2 ERR: %f\n", l2err);

}

#endif//ERROR_HELPERS_HPP
