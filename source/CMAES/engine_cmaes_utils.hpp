#ifndef ENGINE_CMAES_UTILS_HPP
#define ENGINE_CMAES_UTILS_HPP

/* shift arguments */
static int shift(int *c, char ***v) {
    (*c)--; (*v)++;
    return (*c) > 0;
}

/* parse optional arguments and "eats" them */
static void parse(int *c, char ***v, Args *a) {
    shift(c, v); // skip executable
    if (*c && (0 == strcmp(**v, "-cr"))) {
        a->restart = 1;
        shift(c, v);
    } else {
        a->restart = 0;
    }
}

#endif /* ENGINE_CMAES_UTILS_HPP */
