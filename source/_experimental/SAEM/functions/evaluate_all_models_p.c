//function [ ll, sserr ] = evaluate_all_models_p( i, z, par, data, model )


#ifndef max
#define max(a,b) (a)>(b)?(a):(b)
#endif

void evaluate_all_models_p(int idx, double *z, par_t *par, data_t *data, model_t *model, double *ll, double *sserr)
{
//	%
//	% Evaluate all models and return the log-likelihood for each one.
//	%
//	%   i : the i-th individual
//	%

	//%theta = z(1:par.Nmp);
	double theta[par->Nmp];

	for (int i = 0; i < par->Nmp; i++)
		theta[i] = z[i];

#if VERBOSE
	print_matrix("theta", theta, par->Nmp);
#endif
	double error;
	if(strcmp(model->error,"ind") == 0)
	{
		error = z[par->N];
	}
	else
	{
		error = par->alpha;
	}

//	theta = model->par_transf( theta );
//	error = model->err_par_transf( error );
	int deriv = 0;
	model->par_tranf( theta, par->Nmp, deriv, theta);
	model->err_par_tranf( &error, 1, deriv, &error);

//	print_matrix("theta2", theta, par->Nmp);

//	ind = ( data->id == data->uid[i] );
	int ind_size = 0;
	for (int i = 0; i < data->Nd; i++)
	{
		if (data->id[i] == data->uid[idx])
			ind_size++;
	}

#if VERBOSE
	printf("ind_size = %d\n", ind_size);
#endif

	int ind[ind_size];
	for (int i = 0, k = 0; i < data->Nd; i++)
	{
		if (data->id[i] == data->uid[idx])
		{
			ind[k] = i;
			k++;
		}
	}

//	print_matrix_i("ind", ind, ind_size);

	double y[ind_size];

#if VERBOSE
	double x[ind_size];
	for (int i = 0; i < ind_size; i++)
		x[i] = data->x[ind[i]];

	print_matrix("x", x, ind_size);
#endif

//	y = model->fun( theta, data->x(ind) );
	for (int i = 0; i < ind_size; i++)
	{
		model->fun(theta, &data->x[ind[i]], 1, &y[i]);	// peh: check second and last arguments
	}

#if VERBOSE
	print_matrix("y", y, ind_size);
#endif

//	if( any(~isfinite(y)) )
//		error('Not finite model evaluation')
//	end
	for (int i = 0; i < ind_size; i++)
	{
		if (isinf(y[i]))
		{
			printf("Error: Not finite model evaluation\n");
			abort();
		}
	}

	double sigma2 = error*error;

	double tmp[ind_size];
	if (strcmp(model->error_model, "constant") == 0)
	{
		//*sserr = sum( (y-data.y(ind)).^2 );
		*sserr = 0;
		for (int i = 0; i < ind_size; i++)
			*sserr += pow(y[i]-data->y[ind[i]],2);

#if DEBUG
		printf("sserr = %lf\n", *sserr);
#endif
		// tmp = ones(1,data.Ndi(i));
		for (int i = 0; i < ind_size; i++)
			tmp[i] = sigma2*1.0; 	// peh: ind_size == data->Ndi[idx]

	}
	else if (strcmp(model->error_model, "proportional") == 0)
	{
		double eps = 1e-6;

		//y2 = max(eps,y*y);
		double y2[ind_size];
		for (int i = 0; i < ind_size; i++)
		{
			y2[i] = max(eps,y[i]*y[i]);
		}

		//*sserr = sum(  ((y-data.y(ind)).^2) ./ y2   );
		*sserr = 0;
		for (int i = 0; i < ind_size; i++)
			*sserr += (pow(y[i]-data->y[ind[i]],2) / y2[i]);

		//tmp = sigma2*y2;
		for (int i = 0; i < ind_size; i++)
			tmp[i] = sigma2*y2[i];

	}
	else
	{
		printf("Unknown error model\n");
	}

	double sum_log_tmp = 0.0;

	for (int i = 0; i < ind_size; i++)
	{
		sum_log_tmp += log(tmp[i]);
	}

#if VERBOSE
	printf("sum_log_tmp = %e\n", sum_log_tmp);
	printf("sigma2 = %e\n", sigma2);
#endif

	double log2pi = 0.5*log(2*M_PI);

#if VERBOSE
	printf("log2pi = %e\n", log2pi);
	printf("data_Ndi = %d\n", data->Ndi[idx]);
#endif

	if( !isinf(*sserr))
	{
		*ll = -data->Ndi[idx]*log2pi - 0.5*sum_log_tmp - 0.5*(*sserr)/sigma2;

#if VERBOSE
		printf("*ll = %lf\n", *ll);
#endif
	}
	else
	{
		*ll = -1e200;
	}

	return;
}
