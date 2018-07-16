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
#include <fitfun.h>

typedef struct {
    // empty for now
} Args;

/* shift arguments */
static int shift(int *c, char ***v) {
    (*c)--; (*v)++;
    return (*c) > 0;
}

/* parse optional arguments and "eats" them */
static void parse(int *c, char ***v, Args *a) {
    shift(c, v); // skip executable
}


int main(int argc, char *argv[]) {
    Args a;
    Params par;
    dram_init(&par);

    parse(&argc, &argv, &a);
    fitfun_initialize(argc, argv);

    dram(&par);

    dram_finalize(&par);
    fitfun_finalize();

    return 0;
}
