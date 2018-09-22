//function y = my_model(x,theta)


void my_model(double *theta, double *x, int FLAG, double *y)
{
        double res = 0.0;

#if VERBOSE
        printf("theta = %lf %lf %lf\n", theta[0], theta[1], theta[2]);
        printf("x[0] = %lf\n", x[0]);
#endif
	
	//y = theta^2;
	//y(~isfinite(y)) = 1e300;

	res = theta[0]*theta[0];
	if (isinf(res)) res = 1e300;

#if VERBOSE
        printf("res = %lf\n", res);
#endif
        *y = res;
}

