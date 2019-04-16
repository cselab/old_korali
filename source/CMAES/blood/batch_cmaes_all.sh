#!/bin/bash

for d in {3..7}
do

rm cmaes_gp_blood.o 
rm cmaes_gp_blood

make DATA=${d} euler=1

for i in {1..6}
do
    echo "HI"
    #bsub -n 1 -W 04:00 ./cmaes_gp_blood ${i}
done

done
