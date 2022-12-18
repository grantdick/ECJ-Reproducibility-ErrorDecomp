#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#ifdef  __cplusplus
extern "C" {
#endif

    #include "parsetree.h"

    struct ensemble; /* opaque structure - definition not needed by external users */

    struct ensemble *ensemble_create();

    void ensemble_free(struct ensemble *E);

    void ensemble_add_model(struct ensemble *E, struct parsetree *model, double a, double b);

    double ensemble_predict(const struct ensemble *E, const double *X);

    double ensemble_mse(const struct ensemble *E, const double * const *X, const double *t, int n);

    int ensemble_total_size(const struct ensemble *E);

    int ensemble_unique_size(const struct ensemble *E);

#ifdef  __cplusplus
}
#endif

#endif
