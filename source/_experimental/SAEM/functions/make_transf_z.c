//%function transf_z = make_transf_z( paramtransform, deriv )

//	model->par_transf     = make_transf_z( par_transf );
//	model->err_par_transf = make_transf_z( error_par_transf );

void transf_zselector(double *phi, int size_phi, int *tr, double *psi)
{
	for (int i = 0; i < size_phi; i++)
		psi[i] = phi[i];

	for (int i = 0; i < size_phi; i++)
	{
		if (tr[i]==1)
			psi[i]=exp(phi[i]);		// lognormal
		if (tr[i]==2)
			psi[i]=normcdf(phi[i]);		// probit
		if (tr[i]==3)
			psi[i]=1/(1+exp(-phi[i]));	// logit
	}
}

//% Utility for computing the deriviative of the transform when there is a mix of types
void dtransf_zselector(double *phi, int size_phi, int *tr, double *d_psi)
{
	for (int i = 0; i < size_phi; i++) d_psi[i] = 1;
	for (int i = 0; i < size_phi; i++)
	{
		if (tr[i]==1)
			d_psi[i]=exp(phi[i]);				// lognormal
		if (tr[i]==2)
			d_psi[i]=normpdf(phi[i]);			// probit
		if (tr[i]==3)
			d_psi[i]=1/(2+exp(-phi[i])+exp(phi[i]));	// logit
	}
}


void transf_z(double *phi, int size_phi, int deriv, double *out)
{
	int *paramtransform = par.transf;

#if VERBOSE
	for (int i = 0; i < size_phi; i++)
		printf("paramtransform[%d] = %d\n", i, paramtransform[i]);
#endif

	if (deriv==0)
		transf_zselector(phi, size_phi, paramtransform, out);
	else
		dtransf_zselector(phi, size_phi, paramtransform, out);

//	if nargin<2 || ~deriv
//		transf_z = @(x,ind) transf_zselector(x,paramtransform);
//	else
//		transf_z = @(x,ind) dtransf_zselector(x,paramtransform);
}


void transf_z_ex(double *phi, int *paramtransform, int size_phi, int deriv, double *out)
{
//	int *paramtransform = par.transf;

#if VERBOSE
	for (int i = 0; i < size_phi; i++)
		printf("paramtransform[%d] = %d\n", i, paramtransform[i]);
#endif

	if (deriv==0)
		transf_zselector(phi, size_phi, paramtransform, out);
	else
		dtransf_zselector(phi, size_phi, paramtransform, out);

//	if nargin<2 || ~deriv
//		transf_z = @(x,ind) transf_zselector(x,paramtransform);
//	else
//		transf_z = @(x,ind) dtransf_zselector(x,paramtransform);
}


void error_transf_z(double *phi, int size_phi, int deriv, double *out)
{
	int *paramtransform = &par.error_transf;
	if (deriv==0)
		transf_zselector(phi, size_phi, paramtransform, out);
	else
		dtransf_zselector(phi, size_phi, paramtransform, out);

//	if nargin<2 || ~deriv
//		transf_z = @(x,ind) transf_zselector(x,paramtransform);
//	else
//		transf_z = @(x,ind) dtransf_zselector(x,paramtransform);
}
