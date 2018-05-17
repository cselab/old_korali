#Coupling with external code

```sh
  cd build
  make tmcmc_theta_external use_torc=1
```

This version creates intermediate folders of the form `tmpdir.*.*.*.*` where the external code is running on a set of parameters saved in `params.txt`. After the execution the directory is deleted. If you want to keep the directories set `REMOVEDIRS  0` in `source/likelihoods/loglike_theta.c`.


**To setup the example:**
```sh
	cd examples/sampling/external/tmcmc
	./setup_tmcmc.sh
	cd runs/run_001
```

**Contents of `runs/run_001`:**  

- The Π4U executable: `tmcmc_theta_external`  
- The files `priors.par` and `tmcmc.par`. These files contain information about the choice of the prior distributions and the TMCMC parameters respectively.  

- A directory named `model`. Inside this directory, Π4U expects to find:  
    - a file with the experimental data named `data.txt`
    - a user-provided script named `doall.sh`, which (i) runs the external simulation, (ii) compares the output with the experimental data, and (iii) saves the log-likelihood inside a file called `loglike.txt`. The value inside `loglike.txt` is then read from Π4U.


**To run the example:** 
```sh
	mpirun -np 4 ./tmcmc_theta_external
```

or
```sh
	export TORC_WORKERS=4
	mpirun -np 1 ./tmcmc_theta_external
```
