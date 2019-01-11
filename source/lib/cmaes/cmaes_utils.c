#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cmaes_utils.h"

#if defined(_USE_TORC_)

#include <mpi.h>
#include <torc.h>

double get_time() {
    return torc_gettime();
}

#else

#include <sys/time.h>
double get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}

#endif

#include <sys/stat.h>

int cmaes_utils_file_exists(const char *fname) {
	struct stat buffer;
	return (stat (fname, &buffer) == 0);
}

void cmaes_utils_read_bounds(int verbose, const char *fname, double **p_lower_bound, double **p_upper_bound, int dim) {
    double *lower_bound = malloc(dim*sizeof(double));
    double *upper_bound = malloc(dim*sizeof(double));

    FILE *f = fopen(fname, "r");
    
    if (f != NULL){
      
        printf("Reading the bounds from '%s'\n", fname);

      	char line[256];
      	int found;
      	int line_no = 0;
      	for (int i = 0; i < dim; i++) {
        	
            found = 0;
            while (fgets(line, 256, f)!= NULL) {
          	line_no++;

          	if ((line[0] == '#')||(strlen(line)==0)) continue;

                char bound[32];
                sprintf(bound, "B%d", i);
                if (strstr(line, bound) != NULL) {
                    sscanf(line, "%*s %lf %lf", &lower_bound[i], &upper_bound[i]);
                    found = 1;
                    break;
                }
            }
            if (!found) {
                printf("Bounds for parameters %d not found in '%s'. Exit...'\n", i, fname);
                exit(1);
            }
            rewind(f);
            line_no = 0;
        }
      	fclose(f);
    }
    else {
        printf("Parameters file '%s' could not be opened. Exit...\n", fname);
        exit(1);
    }


    if( verbose ){
    	printf("Parameter Bounds:\n");
    	for (int i = 0; i < dim; i++) {
            printf("B%d: %15.6f %15.6f\n", i, lower_bound[i], upper_bound[i]);
    	}
    }

    (*p_lower_bound) = lower_bound;
    (*p_upper_bound) = upper_bound;
}

double cmaes_utils_load_pop_from_file(int verbose, int step, double * const* pop, double *arFunvals, int dim, int lambda, int * checkp) {	
    char filename[256];
    sprintf(filename, "curgen_db_%03d.txt", step);
    FILE *fp = fopen(filename, "r");
    double tt0, tt1 ;
     
    tt0 = get_time();

    if( fp != NULL ) {
   
        for( int i = 0; i < lambda; i++ ){
            for( int j = 0; j < dim; j++ ){
                int r = fscanf(fp, "%le", &pop[i][j]);
   
                if (verbose) printf("[%d] pop[%d][%d]=%f\n", r, i, j, pop[i][j]);
                 
                if (r < 0){
                    printf("Error occured while reading the (%d,%d) element from %s. Exit...\n",i,j,filename);
                    exit(1);
                }
            }
         	
            int r = fscanf(fp, "%le", &arFunvals[i]);
         	
            if (verbose) printf("[%d] arFunvals[%d] = %f\n", r, i, arFunvals[i]);
			
            if (r < 0){
                printf("Error occured while reading the (%d) function value from %s. Exit...\n",i, filename);
                exit(1);
            }
        }
    	fclose(fp);
    }
    else
        *checkp = 0;

    tt1 = get_time();

    return tt1-tt0;
}

static int is_feasible(double *pop, double *lower_bound, double *upper_bound, int dim) {
    int i, good;
    for (i = 0; i < dim; i++) {
        good = (lower_bound[i] <= pop[i]) && (pop[i] <= upper_bound[i]);
        if (!good) {
            return 0;
        }
    }
    return 1;
}

void cmaes_utils_make_all_points_feasible( cmaes_t *evo, double* const *pop, double * lower_bound, double * upper_bound ){

	int lambda = cmaes_Get( evo, "lambda");
    int dim    = cmaes_Get( evo, "dim");

	for( int i=0; i<lambda; ++i)
    	while( !is_feasible( pop[i],lower_bound,upper_bound,dim ) )
            cmaes_ReSampleSingle( evo, i );

}

void cmaes_utils_print_the_best( cmaes_t evo, int step ) {
    int dim    = cmaes_Get( &evo, "dim");
    	
    const double *xbever = cmaes_GetPtr(&evo, "xbestever");
    double        fbever = cmaes_Get(   &evo, "fbestever");

    printf("BEST @ %5d: ", step);
    for( int i = 0; i < dim; i++ )
        printf("%25.16lf ", xbever[i]);
    printf("%25.16lf\n", fbever);
}

void cmaes_utils_write_pop_to_file( cmaes_t evo, double *arFunvals, double * const* pop, int step ){
    int dim    = cmaes_Get( &evo, "dim");
    int lambda = cmaes_Get( &evo, "lambda");
    	
    char filename[256];
    sprintf(filename, "curgen_db_%03d.txt", step);
    FILE *fp = fopen(filename, "w");
    for (int i = 0; i < lambda; i++) {
        for ( int j = 0; j < dim; j++) 
            fprintf(fp, "%.6le ", pop[i][j]);
        fprintf(fp, "%.6le\n", arFunvals[i]);
    }
    fclose(fp);
}

int cmaes_utils_is_there_enough_time(long job_max_time, double gt0, double dt ) {
	
    if (job_max_time <= 0) return 1;
	
    double lastgen_time = dt;
    long runt, remt;
    	
    long maxt = job_max_time;    //get_job_maxTime(); // from lsf or provided by the user
    printf("job maxtime = %ld\n", maxt);

    runt = get_time()-gt0;    //runt = get_job_runTime(); // from lsf or provided by the application: runt = omp_get_wtime()-gt0;
    remt = maxt - runt;
    printf("job runtime = %ld remtime = %ld\n", runt, remt);

    if ((lastgen_time*1.1) > remt){
        printf("No more available time, exiting...\n");
        return 0;
    }

    return 1;
}
