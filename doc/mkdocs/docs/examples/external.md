#Coupling with external code

In this example we run the sampling with two external parallel libraries, namely MPI and TORC. 

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
    - a user-provided script named `doall.sh`, which (i) runs the external simulation, (ii) compares the output with the experimental data, and (iii) saves the log-likelihood inside a file called `loglike.txt`. The value inside `loglike.txt` is then read from korali.


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
