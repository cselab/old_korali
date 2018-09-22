#include <stdio.h>
#include "../auxil.c" 

int main(int argc, char *argv[])
{
	#define M_ll	4
	#define Ni	4

	double LL[M_ll][Ni] = { {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};

	for (int i = 0; i < M_ll; i++)
	for (int j = 0; j < Ni; j++)
		LL[i][j] = exp(LL[i][j]);

	for (int i = 1; i < M_ll; i++)
	for (int j = 0; j < Ni; j++)
		LL[i][j] += LL[i-1][j];

	for (int i = 0; i < M_ll; i++)
	for (int j = 0; j < Ni; j++)
		LL[i][j] /= (i+1);

	double LLout[Ni];

	for (int j = 0; j < Ni; j++)
	{
		LLout[j] = 0;
		for (int i = 0; i < M_ll; i++)
			LLout[j] += log(LL[i][j]);
	}


	//LL = bsxfun( @rdivide, cumsum(exp(LL),1), (1:size(LL,1))' ); //'
	//LL = sum(log(LL),2);

	print_matrix_2d_linear("LL", (double *)LL, M_ll, Ni);

	print_matrix("LLout", LLout, Ni);

	return 0;
}


