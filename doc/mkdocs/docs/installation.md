# Installation

## Prerequisities
- An MPI implementation must be installed on your system, preferably with full thread safety, e.g., [MPICH](http://www.mpich.org).
- [GSL-2.4](http://www.gnu.org/software/gsl/) or later must be installed on your system.


## Installation steps


### 1. Torc library
A tasking library that allows to write platform-independent code. We assume that the MPI compiler is named mpicc:
```sh
	cd lib/torc_lite  
	autoreconf  
	./configure CC=mpicc --prefix=$HOME/usr/torc  
	make; make install  
	export PATH=$HOME/usr/torc/bin:$PATH  
```
After installing torc, the following flags are available:  
`torc_cflags`  
`torc_libs`  






### 2. GSL library

The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers.

- Download the latest version of GSL. In a terminal write
```sh
	wget http://mirror.switch.ch/ftp/mirror/gnu/gsl/gsl-latest.tar.gz
	tar -xvzf gsl-latest.tar.gz
```

- Configure. Set the install folder to be `/$HOME/usr`. If you want to install in the default diretcory, `/usr/local`, delete the  `--prefix=/$HOME/usr`.
```sh
	./configure   --prefix=/$HOME/usr
```

- Compile and install
```sh
	make -j2
	make install
```
This step will take some time. If you have more available cores, change the 2 in the -j2 to a bigger number.




### 3. Sampling and Optimization Algorithms

Enter the `build` directory:  
```sh
	cd build  
```




Before compiling, the following need to be checked:

- Path to GSL-2.4 inside the `Makefile`.  
- Name for the MPI compiler in the `Makefile`, since this can be named differently on different platforms (e.g. CC=cc on Piz Daint).  

**Compilation options:**  

build the default option (uses use_torc=0):
```sh
	make
```

build the OpenMP version:
```sh
	make use_omp=1
```

build the TORC-based version:
```sh
	make use_torc=1
```

build the serial version:
```sh
	make  use_omp=0  use_torc=0
```

## Test

coming soon


## Notes

Please send your questions to:

- chatzidp AT ethz.ch
- garampat AT ethz.ch
