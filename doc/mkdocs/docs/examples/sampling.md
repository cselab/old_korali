# Sampling a posterior distribution

In this example we will show how to sample the posterior distribution in an uncertainty quantification problem. First, we will create synthetic data using the model,

$$
f(t;\varphi) = \varphi_1 \sin(\varphi_2  x + \varphi_3 )
$$

We fix $\varphi^* = (2,3,1)$ and create $100$ data point using the equation,

$$
d_i = f(x_i,\varphi^{\star}) + \sigma \epsilon, \quad \epsilon \sim \mathcal{N}(0,1) \, ,
$$

where $x_i = 0.02 i,\; i=1,\ldots,100$ and $\sigma=0.3$. We will sample the posterior distribution of $\vartheta=(\varphi,\sigma)$ conditioned on the data $d$. The prior distribution is uniform for each parameter,

\begin{align}
    p(\vartheta_1) &= \mathcal{U}( \vartheta_1 | 0,5) \\
    p(\vartheta_2) &= \mathcal{U}( \vartheta_2 | 0,10) \\
	p(\vartheta_3) &= \mathcal{U}( \vartheta_3 | -3.14,3.14) \\
	p(\vartheta_4) &= \mathcal{U}( \vartheta_4 | 0,5)  \, ,
\end{align}

and the likelihood function in given by,

\begin{align}
    p(d | \vartheta) & = \prod_{i=1}^4 p(d_i | \vartheta) \\
					 &=  \prod_{i=1}^4 \mathcal{N}( d_i | f(x;\varphi),\sigma ) \, .
\end{align}



## Sampling with TMCMC

### Compile and run

From the base folder run
```sh
cd build
make tmcmc_theta_internal
```

Make sure that `use_torc=0` and `use_omp=0` in the Makefile since we don't want to run parallel in this example. Go back to the the base folder and run

```sh
cd ../examples/sampling/internal/tmcmc/
./setup_tmcmc.sh
cd runs/run_001/
```

Run the TMCMC sampling algorithm:
```sh
./tmcmc_theta_internal
```

Finally, visualize the samples:
```sh
cp ../../../../../../source/tools/display/plotmatrix_hist.py .
./plotmatrix_hist.py final.txt
```

![](fig_tmcmc.png)



### Behind the scripts
The script `./setup_tmcmc.sh` makes a new running folder and copies inside the executable and the needed files:

1. the data file data.txt,
1. the file that contains the prior information, [priors.par](../developing/par_files.md#priors.par),
1. the parameter file for tmcmc, [tmcmc.par](../developing/par_files.md#tmcmc.par).

For this example, TMCMC has been linked to the likelihood function `loglike_theta_fast.c`. Inside this file the model $f$ as well as the likelihood function $p(d | \vartheta)$ has beem implemented. More information on the the likelihood implementation and how to write your own likelihood you can find [here](../developing/likelihoods.md).






## Sampling with DRAM

### Compile and run

From the base folder run
```sh
cd build
make dram_theta_internal
```

Make sure that `use_torc=0` and `use_omp=0` in the Makefile since we don't want to run parallel in this example. Go back to the the base folder and run

```sh
cd ../examples/sampling/internal/dram/
./setup_dram.sh
cd runs/run_001/
```

Run the TMCMC sampling algorithm:
```sh
./tmcmc_dram_internal
```

Finally, visualize the samples:
```sh
cp ../../../../../../source/tools/display/plotmatrix_hist.py .
tail -n +1000 chain.txt > tmp
./plotmatrix_hist.py tmp
```
With this command `tail -n +1000 chain.txt` we discard the first $1000$ samples which we consider as the burn-in period.

![](fig_dram.png)



### Behind the scripts
Same as in the TMCMC example. The only difference is that instead of the tmcmc.par is substituted by [dram.par](../developing/par_files.md#dram.par).
