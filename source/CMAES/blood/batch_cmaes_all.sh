for i in {1..6}
do
    rm cmaes_gp_blood
    make
    bsub -n 1 -W 04:00 ./cmaes_gp_blood ${i}
done

