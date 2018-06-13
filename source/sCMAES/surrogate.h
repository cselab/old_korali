struct Surrogate;

void surrogate_ini(Surrogate **s);
void surrogate_fin(Surrogate *s);

void surrogate_add_point(const double *x, double y, Surrogate *s);
double surrogate_eval(const double *x, Surrogate *s);
double surrogate_error(const double *x, Surrogate *s);

