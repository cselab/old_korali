#ifndef SUROOGATE_STORAGE_H
#define SUROOGATE_STORAGE_H

#include <cmaes.h>

typedef struct Surrogate_pop Surrogate_pop;
typedef struct Archive Archive;

void archive_ini(int dim, int nmax, Archive **ap);
void archive_fin(Archive *a);

void archive_add(Archive *a, int n, double *const* pop, double *funvals);
void archive_mark_candidates(Archive *a, double r, cmaes_t *t);


void surrogate_pop_ini(int dim, int nmax, Surrogate_pop **sp);
void surrogate_pop_fin(Surrogate_pop *sp);

void surrogate_pop_select_from_archive(Surrogate_pop *s, Archive *a);

int            surrogate_pop_get_n      (Surrogate_pop *s);
double *const* surrogate_pop_get_pop    (Surrogate_pop *s);
double *       surrogate_pop_get_funvals(Surrogate_pop *s);

void archive_shift_fvals(Archive *a, int n, double *fvals);

void surrogate_pop_set(Surrogate_pop *s, int n, double *const* pop);
void surrogate_pop_transform_coords(Surrogate_pop *s, cmaes_distr_t *d);

#endif // SUROOGATE_STORAGE_H
