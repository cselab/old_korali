# Installation

## Prerequisities
- An MPI implementation must be installed on your system, preferably with full thread safety, e.g., [MPICH](http://www.mpich.org).
- [GSL-2.4](http://www.gnu.org/software/gsl/) or later must be installed on your system.


## Installation steps


### 0. Download
Download the Korali project from [GitHub](https://github.com/cselab/korali).



### 1. Torc Library
The Torc library is shipped with the Korali code. Torc is tasking library that enables platform-independent code. In this example we assume that the MPI compiler is named mpicc:
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



### 2. GSL Library

The [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) is a numerical library for C and C++ programmers.

- Download and unpack the latest version of GSL (outside of korali). Execute following commands in your terminal:
```sh
	wget http://mirror.easyname.at/gnu/gsl/gsl-latest.tar.gz  
	tar -xvzf gsl-latest.tar.gz
```

- Configure GSL. Set the install folder to be `/$HOME/usr`. If you want to install in the default diretcory, `/usr/local`, delete the  `--prefix=/$HOME/usr`.
```sh
	./configure   --prefix=/$HOME/usr
```

- Compile and install GSL.
```sh
	make -j2
	make install
```
This step will take some time. If your machine has more than two cores, set the 2 in the -j2 to that number number.

#### GSL Installation Hints

- Read the INSTALL notes of GSL.
- Update the path variable of your terminal with the location of the binaries (gsl-config, gsl-histogram and gsl-randist).
- Set the LD_LIBRARY_PATH variable of your terminal to the newly installed lib folder.
 


### 3. Korali Installation

Enter the `build` directory:  
```sh
	cd build  
```


Before compiling, the following needs to be checked:

- Name for the MPI compiler in the `source/make/common.mk`, since this can be named differently on different platforms and different compilation options (e.g. mpicc vs gcc) .  

**Compilation options:**  

Build the default option (uses use_torc=0):
```sh
	make
```

Build the OpenMP version:
```sh
	make use_omp=1
```

Build the TORC-based version:
```sh
	make use_torc=1
```

Build the serial version:
```sh
	make  use_omp=0  use_torc=0
```

## Test

coming soon


## Notes

Please send your questions to:

- garampat AT ethz.ch
- wadaniel AT ethz.ch
