#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <math.h>

#include "alloc.h"
#include "readline.h"

void read_data(const char *src,
               double ***X_ptr, double **t_ptr, double ***intervals_ptr,
               int *n_ptr, int *p_ptr)
{
    int n, p, i, j;

    FILE *data;
    char *tok, *line, *buffer;
    size_t bufsz;

    double **X, *t, **intervals;

    data = fopen(src, "r");

    buffer = NULL;
    bufsz = 0;
    line = next_line(&buffer, &bufsz, data);
    p = 0;
    tok = strtok(line, " \t");
    while (tok) {
        p++;
        tok = strtok(NULL, " \t");
    }
    rewind(data);

    n = 0; p--; /* to account for response */
    X = NULL;
    t = NULL;
    line = next_line(&buffer, &bufsz, data);
    while (!feof(data)) {
        if (strlen(line) > 0) {
            X = REALLOC(X, n + 1, sizeof(double *));
            X[n] = MALLOC(p, sizeof(double));
            t = REALLOC(t, n + 1, sizeof(double));

            for (i = 0, tok = strtok(line, " \t"); i < p; ++i, tok = strtok(NULL, " \t")) X[n][i] = atof(tok);
            t[n] = atof(tok);
            n++;
        }
        line = next_line(&buffer, &bufsz, data);
    }
    X = REALLOC(X, n + 1, sizeof(double *));
    X[n] = NULL;

    fclose(data);

    free(buffer);

    if (intervals_ptr) {
        intervals = MALLOC(p + 1, sizeof(double *));

        for (j = 0; j < p; ++j) {
            intervals[j] = MALLOC(2, sizeof(double));

            intervals[j][0] = INFINITY;
            intervals[j][1] = -INFINITY;
            for (i = 0; i < n; ++i) {
                if (X[i][j] < intervals[j][0]) intervals[j][0] = X[i][j];
                if (X[i][j] > intervals[j][1]) intervals[j][1] = X[i][j];
            }
        }

        intervals[p] = NULL;

        *intervals_ptr = intervals;
    }

    *X_ptr = X;
    *t_ptr = t;

    *n_ptr = n;
    *p_ptr = p;
}





void release_data(double **X, double *t, double **intervals)
{
    int i;

    if (intervals) {
        i = 0;
        while(intervals[i]) free(intervals[i++]);
        free(intervals);
    }

    i = 0;
    while (X[i]) free(X[i++]);
    free(X);

    free(t);
}




void scale_data(double **trainX, int ntrain, int p,
                double **testX, int ntest,
                bool bias_correct,
                double **mu_ptr, double **s_ptr)
{
    double delta, *mu = MALLOC(p, sizeof(double)), *s = MALLOC(p, sizeof(double));
    int i, j;

    for (j = 0; j < p; ++j) {
        mu[j] = s[j] = 0;
        for (i = 0; i < ntrain; ++i) {
            delta = trainX[i][j] - mu[j];
            mu[j] += delta / (i + 1);
            s[j] += delta * (trainX[i][j] - mu[j]);
        }

        s[j] = bias_correct ? sqrt(s[j] / (ntrain - 1)) : sqrt(s[j] / ntrain);

        for (i = 0; i < ntrain; ++i) trainX[i][j] = (trainX[i][j] - mu[j]) / s[j];
        for (i = 0; i < ntest; ++i) testX[i][j] = (testX[i][j] - mu[j]) / s[j];
    }

    if (mu_ptr) *mu_ptr= mu; else free(mu);
    if (s_ptr) *s_ptr= s; else free(s);
}





void configure_data(const char *split_file, int split_idx,
                    double **X, double *t, int n, int p,
                    double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                    double ***testX_ptr, double **testt_ptr, int *ntest_ptr)
{
    FILE *splits = fopen(split_file, "r");
    char *line = NULL, *buffer = NULL;
    size_t bufsz = 0;

    int i;
    double **trainX = NULL, *traint = NULL;
    double **testX = NULL, *testt = NULL;
    int ntrain = 0, ntest = 0;

    for (i = 0; i < split_idx; ++i) line = next_line(&buffer, &bufsz, splits);
    fclose(splits);

    for (i = 0; i < n; ++i) {
        if (line[i] == '0') {
            trainX = REALLOC(trainX, ntrain + 1, sizeof(double *));
            trainX[ntrain] = MALLOC(p, sizeof(double));
            memcpy(trainX[ntrain], X[i], p * sizeof(double));

            traint = REALLOC(traint, ntrain + 1, sizeof(double));
            traint[ntrain] = t[i];
            ntrain++;
        } else if (line[i] == '1') {
            testX = REALLOC(testX, ntest + 1, sizeof(double *));
            testX[ntest] = MALLOC(p, sizeof(double));
            memcpy(testX[ntest], X[i], p * sizeof(double));

            testt = REALLOC(testt, ntest + 1, sizeof(double));
            testt[ntest] = t[i];
            ntest++;
        } /* otherwise ignore */
    }

    free(buffer);

    trainX = REALLOC(trainX, ntrain + 1, sizeof(double *));
    trainX[ntrain] = NULL;

    testX = REALLOC(testX, ntest + 1, sizeof(double *));
    testX[ntest] = NULL;


    *trainX_ptr = trainX;
    *traint_ptr = traint;
    *ntrain_ptr = ntrain;

    *testX_ptr = testX;
    *testt_ptr = testt;
    *ntest_ptr = ntest;
}





void bootstrap_data(double **X, double *t, int n, int p,
                    double ***bootX_ptr, double **boott_ptr,
                    double (*rnd)(void))
{
    int i, b;

    double **bootX = MALLOC(n + 1, sizeof(double *)), *boott = MALLOC(n, sizeof(double));

    for (i = 0; i < n; ++i) {
        b = (int)(rnd() * n);

        bootX[i] = MALLOC(p, sizeof(double));
        memcpy(bootX[i], X[b], p * sizeof(double));
        boott[i] = t[b];
    }
    bootX[n] = NULL;

    *bootX_ptr = bootX;
    *boott_ptr = boott;
}





void split_data(double **X, double *t, bool *leave_in, int n, int p,
                double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                double ***testX_ptr, double **testt_ptr, int *ntest_ptr)
{
    int i;
    double **trainX, *traint;
    double **testX, *testt;
    int ntrain, ntest;

    trainX = NULL;
    traint = NULL;
    testX = NULL;
    testt = NULL;
    ntrain = ntest = 0;
    for (i = 0; i < n; ++i) {
        if (leave_in[i]) {
            trainX = REALLOC(trainX, ntrain + 1, sizeof(double *));
            trainX[ntrain] = MALLOC(p, sizeof(double));
            memcpy(trainX[ntrain], X[i], p * sizeof(double));

            traint = REALLOC(traint, ntrain + 1, sizeof(double));
            traint[ntrain] = t[i];
            ntrain++;
        } else {
            testX = REALLOC(testX, ntest + 1, sizeof(double *));
            testX[ntest] = MALLOC(p, sizeof(double));
            memcpy(testX[ntest], X[i], p * sizeof(double));

            testt = REALLOC(testt, ntest + 1, sizeof(double));
            testt[ntest] = t[i];
            ntest++;
        }
    }

    trainX = REALLOC(trainX, ntrain + 1, sizeof(double *));
    trainX[ntrain] = NULL;

    testX = REALLOC(testX, ntest + 1, sizeof(double *));
    testX[ntest] = NULL;

    *trainX_ptr = trainX;
    *traint_ptr = traint;
    *ntrain_ptr = ntrain;

    *testX_ptr = testX;
    *testt_ptr = testt;
    *ntest_ptr = ntest;
}
