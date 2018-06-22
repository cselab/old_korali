#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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
    double *w; // work, dimension: dim
    int n;
};

static int min2i(int a, int b) {return a < b ? a : b;}
static double min2d(double a, double b) {return a < b ? a : b;}

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
    double d;
    int i;
    CoordsArray *ca = a->ca;
    
    for (i = 0; i < a->size; ++i) {
        d = cmaes_transform_distance(t, ca->pop[i]);
        /* printf("%g\n", d); */
        a->candidates[i] = d < r;
    }
}



void surrogate_pop_ini(int dim, int nmax, Surrogate_pop **sp) {
    Surrogate_pop *s;
    s = malloc(sizeof(*s));
    coords_array_ini(dim, nmax, &s->ca);
    s->funvals = malloc(nmax * sizeof(double));
    s->w = malloc(dim * sizeof(double));
    s->n = 0;
    *sp = s;
}

void surrogate_pop_fin(Surrogate_pop *s) {
    coords_array_fin(s->ca);
    free(s->funvals);
    free(s->w);
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
    s->n = n;
}

void surrogate_pop_set(Surrogate_pop *s, int n, double *const* pop) {
    int i;
    size_t sz = s->ca->dim * sizeof(double);
    
    for (i = 0; i < n; ++i)
        memcpy(s->ca->pop[i], pop[i], sz);
    s->n = n;
}

static double dot(int n, const double *a, const double *b) {
    double d = 0;
    int i;
    for (i = 0; i < n; ++i) d += a[i] * b[i];
    return d;
}

static void project_eigen_space(cmaes_distr_t *d, double *xr, double *xe) {
    int i, j, n;
    double **Q, *D, *mu, *w;

    n = d->dim;
    Q = d->Q;
    D = d->D;
    mu = d->mu;
    w = d->w;
    
    if (d->flgdiag) {
        for (i = 0; i < n; ++i)
            xe[i] = (xr[i] - mu[i]) / sqrt(D[i]);
    } else {
        for (i = 0; i < n; ++i) {
            w[i] = 0;
            for (j = 0; j < n; ++j) 
                w[i] += (xr[j] - mu[j]) * Q[j][i];
            w[i] /= sqrt(D[i]);
        }
        for (i = 0; i < n; ++i)
            xe[i] = dot(n, w, Q[i]);
    }
}


void surrogate_pop_transform_coords(Surrogate_pop *s, cmaes_distr_t *d) {
    double *w = s->w;
    int i, dim, n;
    dim = s->ca->dim;
    n = s->n;
    size_t sz = dim * sizeof(double);
    
    for (i = 0; i < n; ++i) {
        project_eigen_space(d, s->ca->pop[i], w);
        memcpy(s->ca->pop[i], w, sz);
    }
}

int surrogate_pop_get_n (Surrogate_pop *s) {return s->n;}
double *const* surrogate_pop_get_pop (Surrogate_pop *s) {return s->ca->pop;}
double * surrogate_pop_get_funvals(Surrogate_pop *s) {return s->funvals;}


static double minnd(int n, const double *a) {
    int i;
    if (n <= 0) return 1e9;
    double m = a[0];
    for (i = 1; i < n; ++i) m = min2d(m, a[i]);
    return m;
}

static void shift_fvals(double mina, double mins, int n, double *fvals) {
    double shift = mina - mins;
    int i;
    if (mina <= mins) return;
    for (i = 0; i < n; ++i) fvals[i] += shift;
}

void archive_shift_fvals(Archive *a, int n, double *fvals) {
    double mina, mins;
    mina = minnd(a->size, a->funvals);
    mins = minnd(n, fvals);

    shift_fvals(mina, mins, n, fvals);
}
