#define _GNU_SOURCE

#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <fitfun.h>
#include <priors.h>
#include "loglike_posterior_thetanew.h"



#define BLANKS " \t\n="
#define BUFLEN 1024


typedef struct database_s
{	
	char psi_folder[128];
	char psi_file[128];
	int  Npsi;
	int  Dpsi;
}database;

static database db;

static Density **priors;
static int Npr;



//=============================================================================
//
//
// count the number of lines in file fp  and check if less than N
static int enough_lines(FILE *fp, int N)
{
    char ch;
    int lines = 0;

	rewind(fp);

    while (!feof(fp)){
            ch = fgetc(fp);
            if (ch == '\n')
                lines++;
    }
    rewind(fp);

    return lines >= N;
}



//=============================================================================
//
//
void read_db_file( const char *file )
{
	FILE *fp = fopen( file, "r" );
	if(fp==NULL){
		printf("\n%s could not be opened. Exit...\n", file );
		exit(EXIT_FAILURE);
	}


	size_t bufsize = 1024;
	char * buffer = (char *)malloc(bufsize * sizeof(char));	

	while(  -1 != getline( &buffer, &bufsize, fp) ){

		char * pch = strtok (buffer, BLANKS );
  		while( pch != NULL ){

			if( pch[0]=='#' || pch[0]=='\n' ) //go to the next line. 
				break;

			if( strcmp(pch,"psi_folder")==0 ){
    			pch = strtok (NULL, BLANKS );
				strcpy( db.psi_folder, pch );
				break;
			}

			if( strcmp(pch,"psi_file")==0 ){
    			pch = strtok (NULL, BLANKS );
				strcpy( db.psi_file, pch );
				break;
			}

			if( strcmp(pch,"Npsi")==0 ){
    			pch = strtok (NULL, BLANKS );
				db.Npsi = atoi( pch );
				break;
			}

			if( strcmp(pch,"Dpsi")==0 ){
    			pch = strtok (NULL, BLANKS );
				db.Dpsi = atoi( pch );
				break;
			}

			puts(buffer);
			printf("\nSomething went wrong while reading the theta database file %s. Exit...\n",file);
			exit(EXIT_FAILURE);
		}
	}

	fclose(fp);
	free(buffer);
}



//=============================================================================
//
//
void print_db_file( const char *file )
{	
	printf("================================\n");
	printf("contents of %s\n",file);

	printf("psi folder            : %s \n", db.psi_folder );
	printf("psi file              : %s \n", db.psi_file );
	printf("Number of psi samples : %d \n", db.Npsi );
	printf("Number of psi dim.    : %d \n", db.Dpsi );
	printf("\n");
	printf("================================\n");
}



//=============================================================================
//
//
void loglike_posterior_thetanew_initialize( )
{
	read_db_file("db_psi.par");
	print_db_file("db_psi.par");
	
	// allocate global and local arrays
    double **psi    = (double **)malloc( sizeof(double *) * db.Npsi );
    psi[0] = (double * )malloc( sizeof(double) * db.Npsi * db.Dpsi);
 	for( int i=0; i<db.Npsi; i++)
        psi[i] = ( *psi + db.Dpsi * i);



    char filename[BUFLEN];
    FILE *fp;


    // read psi
    printf("\n3) Reading %d data from psi database. \n", db.Npsi);

    snprintf(filename, BUFLEN, "%s/%s", db.psi_folder, db.psi_file );
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s does not exist! exiting...\n", filename);
        exit( EXIT_FAILURE );
    }

    if( !enough_lines(fp, db.Npsi) ){
		printf("\n\n Error: Number of samples less than %d in file %s. Exit... \n\n", db.Npsi, filename);
 		exit(EXIT_FAILURE);
    }

    for( int i=0; i < db.Npsi; i++) {
        for( int j=0; j<db.Dpsi; j++)
            fscanf(fp, "%lf", &psi[i][j]);
        double ignore;
        fscanf(fp, "%lf", &ignore);
		fscanf(fp, "%lf", &ignore);
    }
    fclose(fp);
    printf("\nSuccesfull reading data from psi database.\n\n");


	// populate the priors
	priors = (Density **) malloc( db.Npsi*sizeof(Density*) );

	read_priors( "priors_theta.par", &priors[0], &Npr );
	reassign_prior( priors[0], Npr, psi[0] );

	for( int i=1; i<db.Npsi; i++ ){
		new_prior_from_prior( &priors[i], priors[0], Npr );
		reassign_prior( priors[i], Npr, psi[i] );
	}


	free(psi[0]);
    free(psi);
}



//=============================================================================
//
//
void loglike_posterior_thetanew_finalize()
{
	for(int i=0; i<db.Npsi; i++)
		delete_prior( priors[i], Npr );
	free(priors);
}



//=============================================================================
//
//
double loglike_posterior_thetanew(double *theta, int n, void *output, int *info)
{
	double sum = 0;
	for (int i = 0; i < db.Npsi; i++){

		double pr_hb  = exp( prior_log_pdf( priors[i], Npr, theta) );

		sum += pr_hb;
	}

	if(sum==0)	return -1e12;

	double out = -log((double)db.Npsi) + log(sum);

	if (isinf(out) || isnan(out)) out = -1e12;

	return out;
}
