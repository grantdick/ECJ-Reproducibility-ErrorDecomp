#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <unistd.h>
#include <string.h>
#include <math.h>

#include "alloc.h"
#include "data.h"
#include "baselearner.h"
#include "baggp.h"
#include "readline.h"
#include "util.h"

#include "rng.h"





int main(int argc __attribute__((unused)), char **argv)
{
    double **X, *y;
    double **trainX, *traint;
    double **bootX, *boott;
    double **testX, *testt, **testy;

    int n, p, ntrain, ntest;
    int i, b, r, D, R, arg_offset;

    struct baselearner *mdl;
    int nboot, pop_size, ngen;
    double pcross, pmutation;
    bool bootstrap = true;

    seed_rng(argv[1], atoi(argv[2]));
    arg_offset = 2;

    D = atoi(argv[arg_offset + 3]);
    R = atoi(argv[arg_offset + 4]);

    read_data(argv[arg_offset + 1], &X, &y, NULL, &n, &p);
    configure_data(argv[arg_offset + 2], D,
                   X, y, n, p,
                   &trainX, &traint, &ntrain,
                   &testX, &testt, &ntest);
    release_data(X, y, NULL);
    scale_data(trainX, ntrain, p, testX, ntest, true, NULL, NULL);

    nboot       = atoi(argv[arg_offset + 5]);
    pop_size    = atoi(argv[arg_offset + 6]);
    ngen        = atoi(argv[arg_offset + 7]);

    pcross      = 0.5;
    pmutation   = 0.5;

    if (nboot == 0) { /* don't use ensemble, model is a single solution */
        nboot = 1;
        bootstrap = false;
    } else if (nboot < 0) { /* don't use bootstrapping, just perform a simple aggregating ensemble */
        nboot *= -1;
        bootstrap = false;
    }
    testy = MALLOC(ntest, sizeof(double *));
    for (i = 0; i < ntest; ++i) testy[i] = MALLOC(nboot, sizeof(double));

    for (r = 0; r < R; ++r) {
        for (b = 0; b < nboot; ++b) {
            if (bootstrap) {
                bootstrap_data(trainX, traint, ntrain, p, &bootX, &boott, next_rnd);
            } else {
                bootX = trainX;
                boott = traint;
            }

            mdl = gp_evolve((const double * const *)bootX, boott, ntrain, p,
                            (const double * const *)testX, testt, ntest,
                            pop_size, ngen, 6,
                            pcross, pmutation,
                            2, 6, 4, -1, 500,
                            NULL, NULL, -1, NULL, -1,
                            false,
                            next_rnd);

            for (i = 0; i < ntest; ++i) testy[i][b] = baselearner_predict(mdl, testX[i]);
            baselearner_release(mdl);
            if (bootstrap) release_data(bootX, boott, NULL);
        }

        for (i = 0; i < ntest; ++i) {
            cummean(testy[i], nboot);
            fprintf(stdout, "%2d %2d %4d %g", D, r + 1, i, testt[i]);
            for (b = 0; b < nboot; ++b) fprintf(stdout, " %g", testy[i][b]);
            fprintf(stdout, "\n");
            fflush(stdout);
        }
    }


    for (i = 0; i < ntest; ++i) free(testy[i]);
    free(testy);

    release_data(testX, testt, NULL);
    release_data(trainX, traint, NULL);

    return EXIT_SUCCESS;
}
