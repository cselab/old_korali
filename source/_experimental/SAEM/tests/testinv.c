#include <stdio.h>
#include <gsl/gsl_linalg.h>


int inverse(double *mat, double *invmat, int N)
{
     int s;

     gsl_matrix_view m   = gsl_matrix_view_array(mat, N, N);
     gsl_matrix_view inv = gsl_matrix_view_array(invmat,N,N);
     gsl_permutation * p = gsl_permutation_alloc (N);

#if 0
     printf("The matrix is\n");
     for (int i = 0; i < N; ++i)
         for (int j = 0; j < N; ++j)
             printf(j==N-1?"%6.3f\n":"%6.3f ", gsl_matrix_get(&m.matrix,i,j));
#endif

     gsl_linalg_LU_decomp (&m.matrix, p, &s);
     gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);

#if 0
     printf("The inverse is\n");
     for (int i = 0; i < N; ++i)
         for (int j = 0; j < N; ++j)
             printf(j==N-1?"%6.3f\n":"%6.3f ",gsl_matrix_get(&inv.matrix,i,j));
#endif

     gsl_permutation_free (p);

     return 0;
}

int main(int argc, char *argv[])
{

     double F[] = { 1.0, 0.6, 0.0,
                         0.0, 1.5, 1.0,
                         0.0, 1.0, 1.0 };

     double invF[] = { 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0 };

     int N = 3;     
     /*
      * Inverse is
      *    1  -1.2   1.2
      *    0   2.0  -2.0
      *    0  -2.0   3.0
      */

     printf("The matrix is\n");
     for (int i = 0; i < N; ++i)
         for (int j = 0; j < N; ++j)
             printf(j==N-1?"%6.3f\n":"%6.3f ", F[i*N+j]);

     inverse(F, invF, N);

     printf("The inverse is\n");
     for (int i = 0; i < N; ++i)
         for (int j = 0; j < N; ++j)
             printf(j==N-1?"%6.3f\n":"%6.3f ", invF[i*N+j]);
	
     return 0;

}
