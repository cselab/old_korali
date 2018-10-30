//function omega_chol = chol_eval( omega )

void chol_eval(double *omega, double *omega_chol, int N)
{
//      [omega_chol,p] = chol( omega, 'lower' );
//      if( p )
//              disp(p)
//              disp(omega)
//              error(' The random effects covariance matrix is not positive definite');
//      end

        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
                omega_chol[i*N+j] = 0.0;

        gsl_matrix *work = gsl_matrix_calloc(N,N);

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                gsl_matrix_set(work, i, j, omega[i*N+j]);
            }

        gsl_set_error_handler_off();

        int p = gsl_linalg_cholesky_decomp(work);
        if (p != GSL_SUCCESS) {
		printf(" The random effects covariance matrix is not positive definite\n");
//              disp(p)
//              disp(omega)
		exit(1);
        }

        for (int i = 0; i < N; i++)
        for (int j = 0; j <= i; j++)
                omega_chol[i*N+j] = work->data[i*N+j];

        gsl_matrix_free(work);
}


