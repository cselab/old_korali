#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "surrogate_samples.h"

typedef struct {
    int nmax, dim;
    double ** pop;
} CoordsArray ;

struct Archive {
    int offset, size;
    CoordsArray *ca;
    double *funvals;
    int *candidates;
};

struct Surrogate_pop {
    CoordsArray *ca;
    double *funvals;
    int n;
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
    a->candidates = malloc(nmax*sizeof(int));
    a->funvals    = malloc(nmax*sizeof(double));
    a->size = 0;
    a->offset = 0;
    *ap = a;
}

void archive_fin(Archive *a) {
    coords_array_fin(a->ca);
    free(a->candidates);
    free(a->funvals);
    free(a);
}

void archive_add(Archive *a, int n, double *const* pop, double *funvals) {
    int i, j;
    CoordsArray *ca = a->ca;
    size_t sz = ca->dim * sizeof(double);
    
    for (i = 0; i < n; ++i) {
        j = (a->offset + i) % ca->nmax;
        memcpy(ca->pop[j], pop[i], sz);
        a->funvals[j] = funvals[i];
    }
    a->offset = (a->offset + n) % ca->nmax;
    a->size   = min2i(a->size + n, ca->nmax);
}

void archive_mark_candidates(Archive *a, double r, cmaes_t *t) {
    double sigma, d;
    int i;
    CoordsArray *ca = a->ca;
    
    sigma = cmaes_Get(t, "sigma");
    
    for (i = 0; i < a->size; ++i) {
        d = cmaes_transform_distance(t, ca->pop[i]);
        d *= sigma;
        a->candidates[i] = d < r;
    }
}



void surrogate_pop_ini(int dim, int nmax, Surrogate_pop **sp) {
    Surrogate_pop *s;
    s = malloc(sizeof(*s));
    coords_array_ini(dim, nmax, &s->ca);
    s->funvals = malloc(nmax * sizeof(double));
    s->n = 0;        
    *sp = s;
}

void surrogate_pop_fin(Surrogate_pop *s) {
    coords_array_fin(s->ca);
    free(s->funvals);
    free(s);
}

void surrogate_pop_select_from_archive(Surrogate_pop *s, Archive *a) {
    int i, n;
    size_t sz = s->ca->dim * sizeof(double);
    
    for (i = n = 0; i < a->size; ++i) {
        if (n >= s->ca->nmax) break;
        if (a->candidates[i]) {
            s->funvals[n] = a->funvals[i];
            memcpy(s->ca->pop[n], a->ca->pop[i], sz);
            ++n;
        }
    }
    /* printf("keep %d points\n", n); */
    s->n = n;
}

int surrogate_pop_get_n (Surrogate_pop *s) {return s->n;}
double *const* surrogate_pop_get_pop (Surrogate_pop *s) {return s->ca->pop;}
double * surrogate_pop_get_funvals(Surrogate_pop *s) {return s->funvals;}
