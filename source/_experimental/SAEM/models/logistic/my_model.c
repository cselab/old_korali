//function y = my_model(x,theta)


void my_model(double *theta, double *x, int FLAG, double *y)
{
        double res = 0.0;

    double f = exp(theta[2]*x[0]);
    res = ( theta[0]*theta[1]*f )/( theta[0] + theta[1]*(f-1) );

	if (isinf(res)) res = 1e300;

#if VERBOSE
        printf("res = %lf\n", res);
#endif
        *y = res;
}

