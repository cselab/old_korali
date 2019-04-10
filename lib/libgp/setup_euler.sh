#!/bin/bash

tmp_dir=$(pwd)

# load the latest gcc
module unload gcc
module load new
module load gcc/6.3.0

echo "WARNING: this script does not work properly, compile code manually with
cmake + make and the move libgp.a to build  folder (see code here)"

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
echo "HallO"
echo $(pwd)

cd build
echo "BA"
echo $(pwd)

# eigen 3 has been downloaded 
cmake ../  -DEIGEN3_INCLUDE_DIR=$tmp_dir/eigen3

# make the library
if make ..; then
	echo "OK"
else
	exit 1
fi

wait

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

