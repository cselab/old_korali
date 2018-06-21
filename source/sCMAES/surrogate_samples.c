#include <stdlib.h>
#include <string.h>

#include "surrogate_samples.h"

typedef struct {
    int nmax, dim;
    double ** pop;
} CoordsArray ;

struct Archive {
    int offset, size;
    CoordsArray *ca;
};

struct Surrogate_pop {
    CoordsArray *ca;
};

static int min2i(int a, int b) {return a < b ? a : b;}

static void coords_array_ini(int dim, int nmax, CoordsArray **cap) {
    CoordsArray *ca;
    int i;
    ca = malloc(sizeof(*ca));
    *cap = ca;
    ca->pop = malloc(nmax * sizeof(ca->pop[0]));
    for (i = 0; i < nmax; ++i) {
        ca->pop[i] = malloc(dim * sizeof(double));
    }
    ca->nmax = nmax;
    ca->dim = dim;    
}

static void coords_array_fin(CoordsArray *ca) {
    int i;
    for (i = 0; i < ca->nmax; ++i)
        free(ca->pop[i]);
    free(ca->pop);
    free(ca);
}

void archive_ini(int dim, int nmax, Archive **ap) {
    Archive *a;
    a = malloc(sizeof(*a));
    coords_array_ini(dim, nmax, &a->ca);
    a->size = 0;
    a->offset = 0;
    *ap = a;
}

void archive_fin(Archive *a) {
    coords_array_fin(a->ca);
    free(a);
}

void archive_add(Archive *a, int n, double *const* pop) {
    int i, j;
    CoordsArray *ca = a->ca;
    size_t sz = ca->dim * sizeof(double);
    
    for (i = 0; i < n; ++i) {
        j = (a->offset + i) % ca->nmax;
        memcpy(ca->pop[j], pop[i], sz);
    }
    a->offset = (a->offset + n) % ca->nmax;
    a->size   = min2i(a->size + n, ca->nmax);
}

void surrogate_pop_ini(int dim, int nmax, Surrogate_pop **sp) {
    Surrogate_pop *s;
    s = malloc(sizeof(*s));
    coords_array_ini(dim, nmax, &s->ca);
    *sp = s;
}

void surrogate_pop_fin(Surrogate_pop *s) {
    coords_array_fin(s->ca);
    free(s);
}
