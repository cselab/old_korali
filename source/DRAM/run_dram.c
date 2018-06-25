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

int main(int argc, char *argv[]) {

	dram_init();
        fitfun_initialize(argc, argv);

	dram();

	dram_finalize();
	fitfun_finalize();

	return 0;
}
