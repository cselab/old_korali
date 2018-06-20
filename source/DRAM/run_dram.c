/*
 *  run_dram.c
 *  Pi4U
 *
 *  Copyright 2018 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "dram.h"
#include "fitfun.h"


void fitfun_init(Params *par){
	char str[12];
	sprintf(str, "%d", par->Npar);
	fitfun_initialize(str);
}




int main(int argc, char *argv[]){

    Params par;


	dram_init(&par);


	fitfun_init(&par);


	dram(&par);


	dram_finalize();
	fitfun_finalize();

	return 0;
}
