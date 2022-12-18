#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#ifdef  __cplusplus
extern "C" {
#endif

    #include "baselearner.h"

    struct ensemble; /* opaque structure - definition not needed by external users */





    struct ensemble *ensemble_create(const int size,
                                     const double *trainy, int ntrain,
                                     const double **testX, const double *testy, int ntest);

    void ensemble_free(struct ensemble *E);

    void ensemble_makeboot(struct ensemble *E, double (*rnd)(void));

    bool ensemble_update(struct ensemble *E, int id, const struct baselearner *model, const double *yhat);

    double ensemble_predict(const struct ensemble *E, const double *X);

    double ensemble_oob_error(const struct ensemble *E);

    double ensemble_ib_error(const struct ensemble *E);

    double ensemble_test_error(const struct ensemble *E);

    double ensemble_mse(const struct ensemble *E, const double * const *X, const double *t, int n);
    
    int ensemble_size(const struct ensemble *E);

    int ensemble_total_size(const struct ensemble *E);

    int ensemble_unique_size(const struct ensemble *E);

    const int *ensemble_bootstrap(const struct ensemble *E, int id);

    void ensemble_print(FILE *out, const struct ensemble *E);

#ifdef  __cplusplus
}
#endif

#endif
