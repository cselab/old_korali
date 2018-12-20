# Hierarchical Bayesian Model


## Synthetic Data

In this example we show how to sample the posterior distribution in an uncertainty quantification problem. First, we create 5 sets of synthetic data using the model below,

$$
f(x;\varphi) = \varphi_1 \sin(\varphi_2  x + \varphi_3 ).
$$

For each set we fix the parameter $\varphi_i^{\star}$ and create $20$ data points using the equation,

$$
d_{i,j} = f(x_j,\varphi_i^{\star}) + \sigma \epsilon, \quad \epsilon \sim \mathcal{N}(0,1) \, ,
$$

where $x_j = 0.1 j,\; j=1,\ldots,20$ and $\sigma=0.3$. The index $i$ indicates the 5 different sets. The parameters for each set were chosen randomly according to

$$
\varphi_i^{\star} = (2,3,1) + 0.4\zeta, \quad \zeta \sim \mathcal{N}(0,1).
$$

These are the 5 data sets:

<img src="../images/data.png" width="50%" height="50%"/>

 We want to sample the posterior distribution of each $\vartheta_i=(\varphi_i,\sigma)$ conditioned on all the data $d=\{d_1,\ldots,d_5\}$.

<img src="../images/Fig_DAG_HB_plate.png" width="40%" height="40%"/>

The prior distribution is uniform for each parameter,

\begin{align}
   p(\vartheta_1) &= \mathcal{U}( \vartheta_1 | 0,5) \\
   p(\vartheta_2) &= \mathcal{U}( \vartheta_2 | 0,10) \\
   p(\vartheta_3) &= \mathcal{U}( \vartheta_3 | -3.14,3.14) \\
   p(\vartheta_4) &= \mathcal{U}( \vartheta_4 | 0,5)  \, ,
\end{align}

<a href="../predictions.png">Notes on HB</a>


## Build executables

In the build folder `/korali/build` run following commands
```sh
make tmcmc_theta_internal
make tmcmc_psi
make tmcmc_posterior_theta_internal
```
in order to execute the steps within Phase 1 - 3.

## Phase 1: sample from the instrumental distribution

In the build folder execute following commands in order to change directory and run phase 1: 

```sh
cd ../examples/hierarchical/internal/
./setup_phase_1.sh
./run_phase_1.sh
```

And then visualize the results analogous to the other examples.

<img src="../images/theta.png" width="100%" height="100%"/>


## Phase 2: sample hyper-parameters

In directory `korali/examples/hierarchical/internal/` run the commands 
```sh
./setup_phase_2.sh
./run_phase_2.sh
```
and then visualize the results. 

<img src="../images/psi.png" width="80%" height="80%"/>


## Phase 3: sample posterior parameters

Run the commands 
```sh
./setup_phase_3.sh
./run_phase_3.sh
```
and then visualize the results. 

<img src="../images/posterior_theta_001.png" width="80%" height="80%"/>
