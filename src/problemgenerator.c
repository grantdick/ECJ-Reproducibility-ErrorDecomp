#include <stdlib.h>

#include <float.h>
#include <math.h>

#include "alloc.h"
#include "util.h"





static double **create_data_matrix(int n, int p)
{
    int i;
    double **X;

    X = MALLOC(n + 1, sizeof(double *));
    for (i = 0; i < n; ++i) X[i] = MALLOC(p, sizeof(double));
    X[n] = NULL;

    return X;
}





static void friedman1(double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                      double ***testX_ptr, double **testt_ptr, int *ntest_ptr,
                      double ***intervals_ptr,
                      int *p_ptr, double (*rnd)(void), double (*rnd_gauss)(double, double))
{
    double x0, x1, x2, x3, x4;
    double **trainX, *traint;
    double **testX, *testt;
    double **intervals;
    int ntrain, ntest, p;
    int i;

    ntrain = 200;
    ntest = 2000;
    p = 10;

    trainX = create_data_matrix(ntrain, p);
    traint = MALLOC(ntrain, sizeof(double));
    testX  = create_data_matrix(ntest, p);
    testt = MALLOC(ntest, sizeof(double));

    for (i = 0; i < ntrain; ++i) {
        x0 = trainX[i][0] = rnd();
        x1 = trainX[i][1] = rnd();
        x2 = trainX[i][2] = rnd();
        x3 = trainX[i][3] = rnd();
        x4 = trainX[i][4] = rnd();
        trainX[i][5] = rnd();
        trainX[i][6] = rnd();
        trainX[i][7] = rnd();
        trainX[i][8] = rnd();
        trainX[i][9] = rnd();
        traint[i] = 10 * sin(M_PI * x0 * x1) + 20 * (x2 - 0.5) * (x2 - 0.5) + 10 * x3 + 5 * x4;
    }

    for (i = 0; i < ntest; ++i) {
        x0 = testX[i][0] = rnd();
        x1 = testX[i][1] = rnd();
        x2 = testX[i][2] = rnd();
        x3 = testX[i][3] = rnd();
        x4 = testX[i][4] = rnd();
        testX[i][5] = rnd();
        testX[i][6] = rnd();
        testX[i][7] = rnd();
        testX[i][8] = rnd();
        testX[i][9] = rnd();
        testt[i] = 10 * sin(M_PI * x0 * x1) + 20 * (x2 - 0.5) * (x2 - 0.5) + 10 * x3 + 5 * x4;
    }

    for (i = 0; i < ntrain; ++i) traint[i] += rnd_gauss(0, 1);
    for (i = 0; i < ntest; ++i) testt[i] += rnd_gauss(0, 1);

    if (intervals_ptr) {
        intervals = create_data_matrix(p, 2);
        for (i = 0; i < p; ++i) {
            intervals[i][0] = 0;
            intervals[i][1] = 1;
        }
    }


    *trainX_ptr = trainX;
    *traint_ptr = traint;
    *ntrain_ptr = ntrain;

    *testX_ptr = testX;
    *testt_ptr = testt;
    *ntest_ptr = ntest;

    *p_ptr = p;

    if (intervals_ptr) *intervals_ptr = intervals;
}





static void friedman2(double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                      double ***testX_ptr, double **testt_ptr, int *ntest_ptr,
                      double ***intervals_ptr,
                      int *p_ptr, double (*rnd)(void), double (*rnd_gauss)(double, double))
{
    double R, f, L, C, Z;
    double delta, meanZ, m2Z, sdZ, n;
    double **trainX, *traint;
    double **testX, *testt;
    double **intervals;
    int ntrain, ntest, p;
    int i;

    ntrain = 200;
    ntest = 2000;
    p = 4;

    trainX = create_data_matrix(ntrain, p);
    traint = MALLOC(ntrain, sizeof(double));
    testX  = create_data_matrix(ntest, p);
    testt = MALLOC(ntest, sizeof(double));

    meanZ = m2Z = 0;
    n = 0;
    for (i = 0; i < ntrain; ++i) {
        R = trainX[i][0] = rnd() * 100;
        f = trainX[i][1] = 20 + rnd() * 260;
        L = trainX[i][2] = rnd();
        C = trainX[i][3] = 1 + rnd() * 10;

        Z = traint[i] = sqrt(R*R + (2*M_PI*f*L - 1 / (2*M_PI*f*C))*(2*M_PI*f*L - 1 / (2*M_PI*f*C)));

        delta = Z - meanZ;
        meanZ += delta / ++n;
        m2Z += delta * (Z - meanZ);
    }

    for (i = 0; i < ntest; ++i) {
        R = testX[i][0] = rnd() * 100;
        f = testX[i][1] = 20 + rnd() * 260;
        L = testX[i][2] = rnd();
        C = testX[i][3] = 1 + rnd() * 10;

        Z = testt[i] = sqrt(R*R + (2*M_PI*f*L - 1 / (2*M_PI*f*C))*(2*M_PI*f*L - 1 / (2*M_PI*f*C)));

        delta = Z - meanZ;
        meanZ += delta / ++n;
        m2Z += delta * (Z - meanZ);
    }

    sdZ = sqrt(m2Z / (n - 1));
    for (i = 0; i < ntrain; ++i) traint[i] += rnd_gauss(0, sdZ) / 3; /* SNR of 3 */
    for (i = 0; i < ntest; ++i) testt[i] += rnd_gauss(0, sdZ) / 3; /* SNR of 3 */

    if (intervals_ptr) {
        intervals = create_data_matrix(p, 2);
        intervals[0][0] =  0; intervals[0][1] = 100;
        intervals[1][0] = 20; intervals[1][1] = 280;
        intervals[2][0] =  0; intervals[2][1] =   1;
        intervals[3][0] =  1; intervals[3][1] =  11;
    }

    *trainX_ptr = trainX;
    *traint_ptr = traint;
    *ntrain_ptr = ntrain;

    *testX_ptr = testX;
    *testt_ptr = testt;
    *ntest_ptr = ntest;

    *p_ptr = p;

    if (intervals_ptr) *intervals_ptr = intervals;
}





static void friedman3(double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                      double ***testX_ptr, double **testt_ptr, int *ntest_ptr,
                      double ***intervals_ptr,
                      int *p_ptr, double (*rnd)(void), double (*rnd_gauss)(double, double))
{
    double R, f, L, C, Z;
    double delta, meanZ, m2Z, sdZ, n;
    double **trainX, *traint;
    double **testX, *testt;
    double **intervals;
    int ntrain, ntest, p;
    int i;

    ntrain = 200;
    ntest = 2000;
    p = 4;

    trainX = create_data_matrix(ntrain, p);
    traint = MALLOC(ntrain, sizeof(double));
    testX  = create_data_matrix(ntest, p);
    testt = MALLOC(ntest, sizeof(double));

    meanZ = m2Z = 0;
    n = 0;
    for (i = 0; i < ntrain; ++i) {
        R = trainX[i][0] = rnd() * 100;
        f = trainX[i][1] = 20 + rnd() * 260;
        L = trainX[i][2] = rnd();
        C = trainX[i][3] = 1 + rnd() * 10;

        Z = traint[i] = atan((2*M_PI*f*L - 1 / (2*M_PI*f*C)) / R);

        delta = Z - meanZ;
        meanZ += delta / ++n;
        m2Z += delta * (Z - meanZ);
    }

    for (i = 0; i < ntest; ++i) {
        R = testX[i][0] = rnd() * 100;
        f = testX[i][1] = 20 + rnd() * 260;
        L = testX[i][2] = rnd();
        C = testX[i][3] = 1 + rnd() * 10;

        Z = testt[i] = atan((2*M_PI*f*L - 1 / (2*M_PI*f*C)) / R);

        delta = Z - meanZ;
        meanZ += delta / ++n;
        m2Z += delta * (Z - meanZ);
    }

    sdZ = sqrt(m2Z / (n - 1));
    for (i = 0; i < ntrain; ++i) traint[i] += rnd_gauss(0, sdZ) / 3; /* SNR of 3 */
    for (i = 0; i < ntest; ++i) testt[i] += rnd_gauss(0, sdZ) / 3; /* SNR of 3 */

    if (intervals_ptr) {
        intervals = create_data_matrix(p, 2);
        intervals[0][0] =  0; intervals[0][1] = 100;
        intervals[1][0] = 20; intervals[1][1] = 280;
        intervals[2][0] =  0; intervals[2][1] =   1;
        intervals[3][0] =  1; intervals[3][1] =  11;
    }

    *trainX_ptr = trainX;
    *traint_ptr = traint;
    *ntrain_ptr = ntrain;

    *testX_ptr = testX;
    *testt_ptr = testt;
    *ntest_ptr = ntest;

    *p_ptr = p;

    if (intervals_ptr) *intervals_ptr = intervals;
}





static void pagie1(double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                   double ***testX_ptr, double **testt_ptr, int *ntest_ptr,
                   double ***intervals_ptr,
                   int *p_ptr, double (*rnd)(void))
{
    double x0, x1;
    double **X, *y;
    double **trainX, *traint;
    double **testX, *testt;
    double **intervals;
    int ntrain, ntest, p;
    int n, i, j, k, *idx;

    n = 676;
    ntrain = n / 10;
    ntest = n - ntrain;
    p = 2;

    X = create_data_matrix(n, p);
    y = MALLOC(n, sizeof(double));

    trainX = MALLOC(ntrain + 1, sizeof(double *));
    traint = MALLOC(ntrain, sizeof(double));

    testX = MALLOC(ntest + 1, sizeof(double *));
    testt = MALLOC(ntest, sizeof(double));

    for (k = 0, i = 0; i < 26; ++i) {
        for (j = 0; j < 26; ++j, ++k) {
            x0 = X[k][0] = -5 + i * 0.4;
            x1 = X[k][1] = -5 + j * 0.4;
            y[k] = 1 / (1 + 1 / (x0*x0*x0*x0)) + 1 / (1 + 1 / (x1*x1*x1*x1));
        }
    }


    idx = sample_m_of_n(ntrain, n, rnd);
    for (i = 0; i < ntrain; ++i) {
        trainX[i] = X[idx[i]];
        traint[i] = y[idx[i]];
        X[idx[i]] = NULL;
    }
    trainX[ntrain] = NULL;

    for (ntest = i = 0; i < n; ++i) {
        if (X[i] == NULL) continue;
        testX[ntest] = X[i];
        testt[ntest++] = y[i];
        X[i] = NULL;
    }
    testX[ntest] = NULL;

    free(idx);
    free(y);
    free(X);

    if (intervals_ptr) {
        intervals = create_data_matrix(p, 2);
        intervals[0][0] = -5; intervals[0][1] = 5;
        intervals[1][0] = -5; intervals[1][1] = 5;
    }

    *trainX_ptr = trainX;
    *traint_ptr = traint;
    *ntrain_ptr = ntrain;

    *testX_ptr = testX;
    *testt_ptr = testt;
    *ntest_ptr = ntest;

    *p_ptr = p;

    if (intervals_ptr) *intervals_ptr = intervals;
}










void generate_data(const char *problem,
                   double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                   double ***testX_ptr, double **testt_ptr, int *ntest_ptr,
                   double ***intervals_ptr,
                   int *p_ptr, double (*rnd)(void), double (*rnd_gauss)(double, double))
{
    if (strncmp(problem, "F1", 2) == 0) {
        friedman1(trainX_ptr, traint_ptr, ntrain_ptr,
                  testX_ptr, testt_ptr, ntest_ptr,
                  intervals_ptr, p_ptr, rnd, rnd_gauss);
    } else if (strncmp(problem, "F2", 2) == 0) {
        friedman2(trainX_ptr, traint_ptr, ntrain_ptr,
                  testX_ptr, testt_ptr, ntest_ptr,
                  intervals_ptr, p_ptr, rnd, rnd_gauss);
    } else if (strncmp(problem, "F3", 2) == 0) {
        friedman3(trainX_ptr, traint_ptr, ntrain_ptr,
                  testX_ptr, testt_ptr, ntest_ptr,
                  intervals_ptr, p_ptr, rnd, rnd_gauss);
    } else if (strncmp(problem, "PAGIE1", 6) == 0) {
        pagie1(trainX_ptr, traint_ptr, ntrain_ptr,
                  testX_ptr, testt_ptr, ntest_ptr,
                  intervals_ptr, p_ptr, rnd);
    } else {
        fprintf(stderr, "%s:%d - ERROR - no generator defined for %s problem, quitting.\n",
                __FILE__, __LINE__, problem);
        exit(EXIT_FAILURE);
    }
}
