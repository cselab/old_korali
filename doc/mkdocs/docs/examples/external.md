#Coupling with external code

In this example we sample the posterior distribution of the parameters of an externally defined model and likelihood function, i.e. the TMCMC library calls a shell script that executes two python scripts for the evaluation of the model and the loglikelihood.
Here we build the code with the tasking library TORC for clusters.

**To build the example:**

```sh
  cd build
  make tmcmc_theta_external use_torc=1
```

**To setup the example:**
```sh
	cd ../examples/sampling/external/tmcmc
	./setup_tmcmc.sh
	cd runs/run_001
```

**Contents of `runs/run_001`:**  

- The korali executable: `tmcmc_theta_external`  
- The files `priors.par` and `tmcmc.par`. These files contain information about the choice of the prior distributions and the TMCMC parameters respectively.  

- A directory named `model`. Inside this directory, korali expects to find:  
    - a file with the experimental data named `data.txt`
    - a user-provided script named `doall.sh`, that is called from within the TMCMC library and which (i) runs the external simulation, (ii) compares the output with the experimental data, and (iii) saves the log-likelihood inside a file called `loglike.txt`. 
	The value inside `loglike.txt` is then read from korali.


**To run the example:**

This example creates intermediate folders of the form `tmpdir.*.*.*.*` whilst running on a set of parameters saved in `params.txt`. After the execution the directory is deleted. If you want to keep the directories set `REMOVEDIRS  0` in `source/likelihoods/loglike_theta.c`.

```sh
	mpirun -np 4 ./tmcmc_theta_external
```

or
```sh
	export TORC_WORKERS=4
	mpirun -np 1 ./tmcmc_theta_external
```

Note that the evaluation of externally defined models slow down the execution of the code. You might want to reduce the population size in `tmcmc.par` in order to quickly test this example.
We recommend to use externally defined models for model evaluations that have a longer execution time such that the runtime is not dominated by switching inbetween processes.