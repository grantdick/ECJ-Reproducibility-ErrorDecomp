#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <string.h>

#include "alloc.h"
#include "util.h"
#include "baselearner.h"
#include "parsetree.h"


struct ensemble {
    int size;

    int **boot;
    bool **oob;

    struct baselearner **model;

    int ntrain;
    double **model_yhat_train;
    const double *trainy;

    int ntest;
    const double **testX;
    const double *testy;
    double **ensemble_yhat_test;

    double *model_oob;
    double test_error;

    double *yhat; /* buffer space for ensemble predictions */
};










struct ensemble *ensemble_create(const int size,
                                 const double *trainy, int ntrain,
                                 const double **testX, const double *testy, int ntest)
{
    int i, j;
    struct ensemble *E;

    E = MALLOC(1, sizeof(struct ensemble));

    E->size = size;

    E->boot = MALLOC(size, sizeof(int *));
    E->boot[0] = CALLOC(size * ntrain, sizeof(int));
    for (i = 0; i < E->size; ++i) E->boot[i] = E->boot[0] + i * ntrain;

    E->oob = MALLOC(size, sizeof(bool *));
    E->oob[0] = CALLOC(size * ntrain, sizeof(bool));
    for (i = 0; i < E->size; ++i) {
        E->oob[i] = E->oob[0] + i * ntrain;
        for (j = 0; j < ntrain; ++j) E->oob[i][j] = true;
    }

    E->model = MALLOC(size, sizeof(struct baselearner *));
    for (i = 0; i < E->size; ++i) E->model[i] = baselearner_create();


    E->ntrain = ntrain;
    E->trainy = trainy;

    E->model_yhat_train = MALLOC(size, sizeof(double *));
    E->model_yhat_train[0] = MALLOC(size * ntrain, sizeof(double));
    for (i = 0; i < E->size; ++i) E->model_yhat_train[i] = E->model_yhat_train[0] + i * ntrain;

    E->ntest = ntest;
    E->testX = testX;
    E->testy = testy;

    if (ntest > 0) {
        E->ensemble_yhat_test = MALLOC(ntest, sizeof(double *));
        E->ensemble_yhat_test[0] = MALLOC(ntest * size, sizeof(double));
        for (i = 0; i < ntest; ++i) E->ensemble_yhat_test[i] = E->ensemble_yhat_test[0] + i * size;
    } else {
        E->ensemble_yhat_test = NULL;
    }

    E->model_oob = MALLOC(size, sizeof(double));
    for (i = 0; i < E->size; ++i) E->model_oob[i] = NAN;

    E->test_error = NAN;

    E->yhat = MALLOC(size, sizeof(double));

    return E;
}





void ensemble_free(struct ensemble *E)
{
    if (E == NULL) return;

    free(E->yhat);

    free(E->model_oob);

    if (E->ntest > 0) {
        free(E->ensemble_yhat_test[0]);
        free(E->ensemble_yhat_test);
    }

    free(E->model_yhat_train[0]);
    free(E->model_yhat_train);

    free(E->oob[0]);
    free(E->oob);

    free(E->boot[0]);
    free(E->boot);

    while (E->size > 0) baselearner_release(E->model[--E->size]);
    free(E->model);

    free(E);
}





void ensemble_makeboot(struct ensemble *E, double (*rnd)(void))
{
    int i, j;

    if (rnd) {
        for (i = 0; i < E->size; ++i) {
            memset(E->boot[i], 0, E->ntrain * sizeof(int));
            for (j = 0; j < E->ntrain; ++j) {
                E->boot[i][j] = (int)(rnd() * E->ntrain);
                E->oob[i][E->boot[i][j]] = false;
            }
        }
    } else {
        /* no bootstrapping - all ensemble members work on all training elements */
        for (i = 0; i < E->size; ++i) {
            for (j = 0; j < E->ntrain; ++j) {
                E->boot[i][j] = j;
                E->oob[i][j] = false;
            }
        }
    }
}





bool ensemble_update(struct ensemble *E, int id, const struct baselearner *model, const double *z)
{
    int j, cmp, n = 0;
    double R, loss = 0;
    bool replace;
    const struct parsetree *pmodel, *cmodel;

    for (j = 0; j < E->ntrain; ++j) {
        if (E->oob[id][j]) {
            R = (E->trainy[j] - baselearner_transform(model, z[j]));
            loss += (R*R - loss) / ++n;
        }
    }

    if (n == 0) { /* nothing left out of bag, so use in-bag error for comparison */
        for (j = 0; j < E->ntrain; ++j) {
            R = (E->trainy[E->boot[id][j]] - baselearner_transform(model, z[E->boot[id][j]]));
            loss += (R*R - loss) / (j + 1);
        }
    }

    cmp = cmp_double(loss, E->model_oob[id]);
    replace = !sane(E->model_oob[id]);
    if (!replace) replace = cmp < 0;
    if (!replace && cmp == 0) {
        cmodel = baselearner_get_model(E->model[id]);
        pmodel = baselearner_get_model(model);
        replace = parsetree_size(pmodel) < parsetree_size(cmodel);
    }

    if (replace) {
        baselearner_release(E->model[id]);

        E->model[id] = baselearner_copy(model);

        for (j = 0; j < E->ntrain; ++j) E->model_yhat_train[id][j] = baselearner_transform(model, z[j]);
        E->model_oob[id] = loss;

        for (j = 0; j < E->ntest; ++j) E->ensemble_yhat_test[j][id] = baselearner_predict(E->model[id], E->testX[j]);

        return true;
    } else {
        return false;
    }
}





double ensemble_predict(const struct ensemble *E, const double *X)
{
    int i, n = 0;

    for (i = 0; i < E->size; ++i) {
        E->yhat[n] = baselearner_predict(E->model[i], X);
        if (sane(E->yhat[n])) n++; /* remove +/-Inf and NaN elements before calculation */
    }

    return median(E->yhat, n);
}





double ensemble_oob_error(const struct ensemble *E)
{
    int b, j, n, noob = 0;
    double R, oob = 0;

    for (j = 0; j < E->ntrain; ++j) {
        n = 0;
        for (b = 0; b < E->size; ++b) {
            if (E->oob[b][j]) E->yhat[n++] = E->model_yhat_train[b][j];
        }

        if (n > 0) {
            R = median(E->yhat, n) - E->trainy[j];
            oob += ((R * R) - oob) / (++noob);
        }
    }

    return (noob == 0) ? NAN : oob; /* MSE */
}





double ensemble_ib_error(const struct ensemble *E)
{
    int b, j, n, nib = 0;
    double R, ib = 0;

    for (j = 0; j < E->ntrain; ++j) {
        n = 0;
        for (b = 0; b < E->size; ++b) {
            if (E->boot[b][j] > 0) E->yhat[n++] = E->model_yhat_train[b][j];
        }

        if (n > 0) {
            R = median(E->yhat, n) - E->trainy[j];
            ib += ((R * R) - ib) / (++nib);
        }
    }

    return (nib == 0) ? NAN : ib; /* MSE */
}





double ensemble_test_error(const struct ensemble *E)
{
    double R, L = 0;
    int b, j, n;

    if (E->ntest <= 0) return NAN;

    for (j = 0; j < E->ntest; ++j) {
        n = 0;
        for (b = 0; b < E->size; ++b) {
            E->yhat[n] = E->ensemble_yhat_test[j][b];
            if (sane(E->yhat[n])) n++; /* remove +/-Inf and NaN elements before calculation */
        }

        R = median(E->yhat, n) - E->testy[j];
        L += ((R * R) - L) / (j + 1);
    }

    return L; /* MSE */
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





int ensemble_size(const struct ensemble *E)
{
    return E->size;
}





int ensemble_total_size(const struct ensemble *E)
{
    int b, n = 0;

    for (b = 0; b < E->size; ++b) n += parsetree_size(baselearner_get_model(E->model[b]));

    return n;
}





int ensemble_unique_size(const struct ensemble *E)
{
    int b, j, n = E->size;
    bool *dup = MALLOC(E->size, sizeof(bool));

    for (b = 0; b < E->size; ++b) dup[b] = false;

    for (b = 0; b < E->size; ++b) {
        if (dup[b]) continue;
        for (j = b + 1; j < E->size; ++j) {
            if (parsetree_equal(baselearner_get_model(E->model[b]), baselearner_get_model(E->model[j]))) dup[j] = true;
        }
    }

    for (b = 0; b < E->size; ++b) if (dup[b]) n--;

    free(dup);

    return n;
}




const int *ensemble_bootstrap(const struct ensemble *E, int id)
{
    if (id < E->size) return E->boot[id];

    fprintf(stderr, "%s:%d - WARNING: requested invalid bootstrap sample %d of %d\n", __FILE__, __LINE__, id, E->size);

    return NULL;
}










void ensemble_print(FILE *out, const struct ensemble *E)
{
    int i, j, n;

    fprintf(out, "MDL:\n");
    for (i = 0; i < E->size; ++i) {
        fprintf(out, "%3d:", i);
        for (j = 0; j < E->ntest; ++j) fprintf(out, " %f", E->ensemble_yhat_test[j][i]);
        fprintf(out, "\n");
    }

    fprintf(out, "MED:");
    for (j = 0; j < E->ntest; ++j) {
        n = 0;
        for (i = 0; i < E->size; ++i) {
            E->yhat[n] = E->ensemble_yhat_test[j][i];
            if (sane(E->yhat[n])) n++; /* remove +/-Inf and NaN elements before calculation */
        }
        fprintf(out, " %f", median(E->yhat, n));
    }
    fprintf(out, "\n");

    for (i = 0; i < E->size; ++i) {
        baselearner_print(out, E->model[i]);
        fprintf(out, "\n");
    }
}
