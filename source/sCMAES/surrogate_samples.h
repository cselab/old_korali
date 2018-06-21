#ifndef SUROOGATE_STORAGE_H
#define SUROOGATE_STORAGE_H

typedef struct Surrogate_pop Surrogate_pop;
typedef struct Archive Archive;

void archive_ini(int dim, int nmax, Archive **ap);
void archive_fin(Archive *a);

void archive_add(Archive *a, int n, double *const* pop);

void surrogate_pop_ini(int dim, int nmax, Surrogate_pop **sp);
void surrogate_pop_fin(Surrogate_pop *sp);


#endif // SUROOGATE_STORAGE_H
