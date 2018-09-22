CC=gcc
CFLAGS= -std=c99 -Wall 
CFLAGS += -O3 
#CFLAGS += -O0 -g -ggdb
LDFLAGS = -lgsl -lgslcblas -lm

all:saem_try compute_likelihood_FIM

saem_try: saem_try.c auxil.c global.h load_logistic.c functions/init_saem.c functions/evaluate_all_models_p.c functions/draw_z_samples_p.c functions/make_transf_z.c functions/logp_eta_eval.c functions/update_parameters.c functions/estimate_loglike_p.c functions/stats.c
	$(CC) $(CFLAGS)  -o saem_try saem_try.c $(LDFLAGS)

compute_likelihood_FIM: compute_likelihood_FIM.c auxil.c global.h load_logistic.c functions/init_saem.c functions/evaluate_all_models_p.c functions/draw_z_samples_p.c functions/make_transf_z.c functions/logp_eta_eval.c functions/update_parameters.c functions/estimate_loglike_p.c functions/stats.c functions/log_tpdf.c functions/estimate_FIM.c functions/d_loglike.c functions/stats.c
	$(CC) $(CFLAGS) -o compute_likelihood_FIM compute_likelihood_FIM.c $(LDFLAGS)

clean:
	rm -f saem_try compute_likelihood_FIM
	rm -rf saem_try.dSYM compute_likelihood_FIM.dSYM
