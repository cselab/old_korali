#include <stdlib.h>
#include <gp.h>

#include "surrogate.h"

typedef libgp::GaussianProcess GP;

struct Surrogate {
    GP *gp;
};

void surrogate_ini(int dim, Surrogate **p_s) {
    Surrogate *s;
    s = (Surrogate*) malloc(sizeof(Surrogate));
    *p_s = s;

    s->gp = new GP(dim, "CovSum ( CovSEiso, CovNoise)");
}

void surrogate_fin(Surrogate *s) {
    delete s->gp;
    free(s);
}

void surrogate_add_point(const double *x, double y, Surrogate *s) {
    s->gp->add_pattern(x, y);
}

double surrogate_eval(const double *x, Surrogate *s) {
    return s->gp->f(x);
}

double surrogate_error(const double *x, Surrogate *s) {
    return s->gp->var(x);
}

