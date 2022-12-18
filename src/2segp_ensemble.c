#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <string.h>

#include "alloc.h"
#include "util.h"
#include "parsetree.h"


struct ensemble {
    int size;
    int capacity;

    struct parsetree **model;

    double *a;
    double *b;

    double *y; /* buffer space for ensemble predictions */
};









struct ensemble *ensemble_create()
{
    int i;
    struct ensemble *E;

    E = MALLOC(1, sizeof(struct ensemble));
    E->size = 0;
    E->capacity = 500;

    E->model = MALLOC(E->capacity, sizeof(struct parsetree *));
    for (i = 0; i < E->capacity; ++i) E->model[i] = NULL;

    E->a = MALLOC(E->capacity, sizeof(double));
    E->b = MALLOC(E->capacity, sizeof(double));

    E->y = MALLOC(E->capacity, sizeof(double));

    return E;
}





void ensemble_free(struct ensemble *E)
{
    if (E == NULL) return;

    free(E->y);

    free(E->b);
    free(E->a);

    while (E->size > 0) parsetree_release(E->model[--E->size]);
    free(E->model);

    free(E);
}




void ensemble_add_model(struct ensemble *E, struct parsetree *model, double a, double b)
{
    if (E->size == E->capacity) {
        E->capacity += 500;
        E->model = REALLOC(E->model, E->capacity, sizeof(struct parsetree *));
        E->a = REALLOC(E->a, E->capacity, sizeof(double));
        E->b = REALLOC(E->b, E->capacity, sizeof(double));
        E->y = REALLOC(E->y, E->capacity, sizeof(double));
    }

    E->model[E->size] = parsetree_copy(model);
    E->a[E->size] = a;
    E->b[E->size] = b;
    E->size++;
}





double ensemble_predict(const struct ensemble *E, const double *X)
{
    int i, n;

    n = 0;
    for (i = 0; i < E->size; ++i) {
        E->y[n] = E->a[i] + E->b[i] * parsetree_predict(E->model[i], X);
        if (sane(E->y[n])) n++; /* remove +/-Inf and NaN elements before calculation */
    }

    return mean(E->y, n);
}





double ensemble_mse(const struct ensemble *E, const double * const *X, const double *t, int n)
{
    int j;
    double r, mse = 0;

    for (j = 0; j < n; ++j) {
        r = t[j] - ensemble_predict(E, X[j]);
        mse += (r * r - mse) / (j + 1);
    }

    return mse;
}





int ensemble_total_size(const struct ensemble *E)
{
    int b, n = 0;

    for (b = 0; b < E->size; ++b) n += parsetree_size(E->model[b]);

    return n;
}





int ensemble_unique_size(const struct ensemble *E)
{
    int i, j, n = E->size;
    bool *dup = MALLOC(E->size, sizeof(bool));

    for (i = 0; i < E->size; ++i) dup[i] = false;

    for (i = 0; i < E->size; ++i) {
        if (dup[i]) continue;
        for (j = i + 1; j < E->size; ++j) {
            if (parsetree_equal(E->model[i], E->model[j])) dup[j] = true;
        }
    }

    for (i = 0; i < E->size; ++i) if (dup[i]) n--;

    free(dup);

    return n;
}