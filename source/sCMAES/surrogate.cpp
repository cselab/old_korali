#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gp.h>
#include <rprop.h>

extern "C" {
#include "surrogate.h"
}

typedef libgp::GaussianProcess GP;

struct Surrogate {
    GP *gp;
    int dim;
};

void surrogate_ini(int dim, Surrogate **p_s) {
    Surrogate *s;
    s = (Surrogate*) malloc(sizeof(Surrogate));
    *p_s = s;

    s->gp = new GP(dim, "CovSum ( CovSEiso, CovNoise)");
    s->dim = dim;
}

void surrogate_fin(Surrogate *s) {
    delete s->gp;
    free(s);
}

void surrogate_reset(Surrogate *s) {
    s->gp->clear_sampleset();
}

void surrogate_add_point(const double *x, double y, Surrogate *s) {
    s->gp->add_pattern(x, y);
}

void surrogate_optimize(Surrogate *s) {
    libgp::RProp rprop;
    int verbose = 0;
    rprop.init();
    rprop.maximize(s->gp, 100, verbose);
}

double surrogate_eval(const double *x, Surrogate *s) {
    return s->gp->f(x);
}

double surrogate_error(int n, double* const* xx, const double *yyexact, Surrogate *s) {
    int i;
    double err, dy, *x;

    err = 0;
    for (i = 0; i < n; ++i) {
        x = xx[i];
        dy = surrogate_eval(x, s);
        dy -= yyexact[i];
        err += dy * dy;

        // err += s->gp->var(x);
    }
    return sqrt(n ? err / n : 0.0);
}

