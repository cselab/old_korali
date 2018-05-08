# Π4U
Π4U is a high performance framework for Bayesian uncertainty quantification of large scale computational models.


# INSTALLATION

## PREREQUISITES
- An MPI implementation must be installed on your system (preferably with full thread safety)
- GSL-2.4 (http://www.gnu.org/software/gsl/) is required by some engines and must be installed on your system


## INSTALLATION STEPS
### 1. Torc library
A tasking library that allows to write platform-independent code.  
We assume that the MPI compiler is named mpicc:

	cd lib/torc_lite  
	autoreconf  
	./configure CC=mpicc --prefix=$HOME/usr/torc  
	make; make install  
	export PATH=$HOME/usr/torc/bin:$PATH  

After installing torc, the following flags are available:  
`torc_cflags`  
`torc_libs`  


### 2. UQ and Optimization Algorithms
Enter the `build` directory:  

	cd build  

Before compiling, the following need to be checked: 
- Path to GSL-2.4 inside the `Makefile`.  
- Name for the MPI compiler in the `Makefile`, since this can be named differently on different platforms (e.g. CC=cc on Piz Daint).  

**Compilation options:**  

build the default option (uses use_torc=1):
	
	make
	
build the OpenMP version:
	
	make use_omp=1
	
build the TORC-based version:

	make use_torc=1
	
build the serial version:
	
	make use_omp=0 use_torc=0



# EXAMPLES
Exmple setups can be found inside the `examples` directory.  
More info will follow.


# TESTING
More info will follow.


# NOTES

Do not hesitate to ask for help and report any problems at:
- chatzidp AT ethz.ch
- garampat AT ethz.ch
