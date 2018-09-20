#! /bin/bash

BASE_DIR="runs"

CURRENT_DIR=$(pwd)

export TORC_WORKERS=4


RUN_DIR="${BASE_DIR}/psi"



if [ ! -d "$RUN_DIR" ]; then
	echo "${RUN_DIR} does not exist"
	exit
fi


cd $RUN_DIR


mpirun -np 1   ./tmcmc_psi
#./sample_psi

mkdir -p ../../data/psi
cp final.txt ../../data/psi/psi.txt


cd $CURRENT_DIR

