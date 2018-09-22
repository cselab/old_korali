#if 0
model = 

              path: '../models/logistic/'
               fun: @(theta,t)my_model(t,theta,1)
             error: 'common'
       error_model: 'constant'
        par_transf: @(x,ind)transf_zselector(x,paramtransform)
    err_par_transf: @(x,ind)transf_zselector(x,paramtransform)
#endif

typedef struct model_s
{
	char path[256];
	void (*fun)(double *theta, double *t, int x, double *y);
	char error[256];
	char error_model[256];
//	void (*par_transf)(double *x, double ind);
//	void (*err_par_transf)(double *x, double ind);
	void (*par_tranf)(double *phi, int size_phi, int deriv, double *out);
	void (*err_par_tranf)(double *phi, int size_phi, int deriv, double *out);
} model_t;


#if 0
data = 

     id: [210x1 double]
      x: [210x1 double]
      y: [210x1 double]
    uid: [10x1 double]
     Ni: 10
     Nd: 210
    Ndi: [21 21 21 21 21 21 21 21 21 21]
#endif

typedef struct data_s
{
	int *id;
	double *x;
	double *y;
	int *uid;
	int Ni;
	int Nd;
	int *Ndi;
} data_t;


#if 0
par = 

          beta: [3x1 double]
         omega: [3x3 double]
         alpha: 1
           Nmp: 3
    vectorized: 0
        transf: [0 0 0]
             N: 3
            Nd: 210
    omega_chol: [3x3 double]
            Ni: 10
             L: 5
          M_ll: 5000
         sigma: [4-D double]
#endif

typedef struct par_s
{
	double *beta;
	double *omega;
	int length_omega;	// peh: addition
	double alpha;
	int Nmp;
	int vectorized;

	int *transf;
	int error_transf;	// peh: addition

	int N;
	int Nd;
	double *omega_chol;
	int Ni;
	int L;
	int M_ll;
	double **sigma;
} par_t;


#if 0
stats = 

  struct with fields:
             num: [5x10x3 double]
           denom: [5x10x3 double]
              ar: [5x10x3 double]
          m_cond: [3x10 double]
          v_cond: [3x10 double]
         samples: [14x3x10x5 double]
      compute_mv: 1
    save_samples: 1
              ll: [5x10 double]
           sserr: [5x10 double]
       sserr_tot: [14x10x5 double]
      N_data_tot: 210
#endif

typedef struct stats_s
{
	double *num;
	double *denom;
	double *ar;
	double *m_cond;
	double *v_cond;
	double *samples;
	int compute_mv;
	int save_samples;
	double *ll;
	double *sserr;
	double *sserr_tot;
	int N_data_tot;
	int Nsteps; // new
} stats_t;


#if 0
MCMC = 

  struct with fields:

	steps: [2 2 2]
	outer_steps: 1
#endif

typedef struct MCMC_s
{
	int steps[3];
	int outer_steps;
} MCMC_t;


typedef struct accum_s
{
	double *s1;
	double *s2; 
	double s3;
} accum_t;


extern model_t model;
extern data_t data;
extern par_t par;
extern MCMC_t MCMC;
extern accum_t accum; 

