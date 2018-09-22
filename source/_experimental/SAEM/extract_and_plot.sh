#!/bin/bash

rm -f par_transf_vec_0.txt par_transf_vec_1.txt par_transf_vec_2.txt \
	sqrt_omega_par_transf_vec_0.txt sqrt_omega_par_transf_vec_1.txt sqrt_omega_par_transf_vec_2.txt \
	err_par_transf_vec.txt sqrt_omega_err_par_transf_vec.txt loglikehood.txt

grep -A4 " par_transf_vec_i =" out.txt | grep "\." | awk 'NR % 3 == 1' > par_transf_vec_0.txt 
grep -A4 " par_transf_vec_i =" out.txt | grep "\." | awk 'NR % 3 == 2' > par_transf_vec_1.txt
grep -A4 " par_transf_vec_i =" out.txt | grep "\." | awk 'NR % 3 == 0' > par_transf_vec_2.txt

grep -A4 "sqrt(omega_par_transf_vec_i) =" out.txt | grep "\." | awk 'NR % 3 == 1' > sqrt_omega_par_transf_vec_0.txt 
grep -A4 "sqrt(omega_par_transf_vec_i) =" out.txt | grep "\." | awk 'NR % 3 == 2' > sqrt_omega_par_transf_vec_1.txt
grep -A4 "sqrt(omega_par_transf_vec_i) =" out.txt | grep "\." | awk 'NR % 3 == 0' > sqrt_omega_par_transf_vec_2.txt

grep -A2 "err_par_transf_vec_i =" out.txt | grep "\." > err_par_transf_vec.txt
grep -A2 "sqrt(omega_err_par_transf_vec_i) =" out.txt | grep "\." > sqrt_omega_err_par_transf_vec.txt

grep loglikelihood out.txt | awk '{print 1*$3}' > loglikehood.txt

gnuplot < plot.gp
