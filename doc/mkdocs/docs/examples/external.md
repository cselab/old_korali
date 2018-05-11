#Coupling with external code

```sh
  cd build
  make tmcmc_theta_external use_torc=1
```

This version creates intermediate folders of the form `tmpdir.*.*.*.*` where the external code is running on a set of parameters saved in `params.txt`. After the execution the directory is deleted. If you want to keep the directories set `REMOVEDIRS  0` in `source/likelihoods/loglike_theta.c`.


To run the an example:
```sh
	cd examples/sampling/external/tmcmc
	./setup_tmcmc.sh
	cd runs/run_001
```

and then
```sh
	mpirun -np 4 ./tmcmc_theta_external
```

or
```sh
	export TORC_WORKERS=4
	mpirun -np 1 ./tmcmc_theta_external
```
