# Optimizing a posterior distribution
In this example we will show how to compute the maximum a posteriori (MAP) estimate of the posterior distribution in an uncertainty quantification problem. First, we will create synthetic data using the model,

$$
f(x;\varphi) = \varphi_1 \sin(\varphi_2  x + \varphi_3 )
$$

We fix $\varphi^{\star} = (2,3,1)$ and create $100$ data point using the equation,

$$
d_i = f(x_i,\varphi^{\star}) + \sigma \epsilon, \quad \epsilon \sim \mathcal{N}(0,1) \, ,
$$

where $x_i = 0.02 i,\; i=1,\ldots,100$ and $\sigma=0.3$. We want to optimize the the posterior distribution of $\vartheta=(\varphi,\sigma)$ conditioned on the data $d$. The prior distribution is uniform for each parameter,

\begin{align}
    p(\vartheta_1) &= \mathcal{U}( \vartheta_1 | 0,5) \\
    p(\vartheta_2) &= \mathcal{U}( \vartheta_2 | 0,10) \\
	p(\vartheta_3) &= \mathcal{U}( \vartheta_3 | -3.14,3.14) \\
	p(\vartheta_4) &= \mathcal{U}( \vartheta_4 | 0,5)  \, ,
\end{align}

and the likelihood function in given by,

\begin{align}
    p(d | \vartheta) & = \prod_{i=1}^{100} p(d_i | \vartheta) \\
					 &=  \prod_{i=1}^{100} \mathcal{N}( d_i | f(x;\varphi),\sigma ) \, .
\end{align}


Finally, we want to compute,
$$
\vartheta^\star = \mathop{\arg\max}\limits_{\vartheta} \,\, p( \vartheta | d) \,.
$$



## Optimizing with CMA-ES
Here we use the [Covariance Matrix Adaptation - Evolution Strategy](https://arxiv.org/pdf/1604.00772.pdf) (CMA-ES) algorithm in order to optimize the posterior distribution $p(d | \vartheta)$.

### Compile and run

From the base folder run
```sh
cd build
make cmaes_theta_internal
```

Make sure that we set use_torc=0 and use_omp=0 in the Makefile since we don't want to run the code in parallel. In the build folder execute following commands in order to setup the example:

```sh
cd ../examples/optimization/internal/
./setup_fast_optimize.sh
cd runs/run_001/
```

Run the CMA-ES optimization algorithm:
```sh
./cmaes_theta_internal
```

The output in the terminal will look like this:

![](images/cmaes-terminal.png)

The algorithm has converged to
```
2.0392035199474865  2.9457365262310957  1.0625543252360856  0.3022538646514062
```

Note that the values $\vartheta^{\star} = (2,3,1,0.3)$ have been used to create the synthetic data.


Finally, visualize the samples:
```sh
cp ../../../../../source/tools/postprocessing_tools/cmaes/cmaplt.py .
python plotmatrix_hist.py final.txt
```

![](images/cmaes.png)






### Behind the scripts

The script `./setup_tmcmc.sh` creates new folders and copies the executable and the required files to target directory:

1. the data file `data.txt` containing $d_i$,
2. the file [priors.par](../developing/par_files.md#priors.par) that contains the prior information $p(\vartheta_i)$ (i.e number of prior distributions, distribution functions, function parameters),
3. the parameter file for cmaes, [cmaes_bounds.par](../developing/par_files.md#cmaes_bounds.par), [cmaes_initials.par](../developing/par_files.md#cmaes_initials.par), [cmaes_signals.par](../developing/par_files.md#cmaes_signals.par).

In this example, Korali has been linked to the file `loglike_theta_fast.c`. Inside this file you can find the model $f$ as well as the likelihood function $p(d | \vartheta)$. More information on the implementation of a likelihood can be found [here](../developing/likelihoods.md).
