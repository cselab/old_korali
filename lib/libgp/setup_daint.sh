#!/bin/bash


# load the latest gcc
module unload gcc
module load gcc/7.1.0



# complie the libgp library
(

cd libgp



if [ -d build ]; then
	echo
	echo "Run cleanup.sh first"
	echo
	exit 1
fi




mkdir build 

cd build 

# eigen 3 has been downloaded 
cmake ../  -DEIGEN3_INCLUDE_DIR=$HOME/test_surrogates/surrogates/eigen3/include/eigen3


# make the library
if make gp; then
	echo "OK"
else
	exit 1
fi


cd ../

# move the library in the lib file
mkdir -p lib
cp build/libgp.a lib/.


# run a test
cd gpreg
if make ; then
	echo "OK"
else
	exit 1
fi
./mytest 100 10 ./data/sincos.txt

)

