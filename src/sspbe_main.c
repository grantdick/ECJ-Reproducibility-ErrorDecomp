#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <unistd.h>
#include <string.h>
#include <math.h>

#include "data.h"
#include "sspbe_ensemble.h"
#include "sspbe.h"
#include "readline.h"
#include "util.h"

#include "rng.h"





int main(int argc __attribute__((unused)), char **argv)
{
    char *problem;

    double **X, *y;
    double **trainX, *trainy, *mu = NULL, *s = NULL;
    double **testX, *testy;
    bool *leave_in;
    int n, p, ntrain, ntest;
    bool verbose = false;
    int i, r, D, R, arg_offset;

    struct ensemble *ens;
    int nboot, pop_size, ngen;
    double pcross, pmutation;

    unsigned long seed;
    seed_rng(argv[1], atoi(argv[2]));
    seed = atoi(argv[2]);
    arg_offset = 2;

    if (argv[arg_offset + 1][0] == 'D') {
        D = atoi(argv[arg_offset + 4]);
        R = atoi(argv[arg_offset + 5]);

        read_data(argv[arg_offset + 2], &X, &y, NULL, &n, &p);
        configure_data(argv[arg_offset + 3], D,
                       X, y, n, p,
                       &trainX, &trainy, &ntrain,
                       &testX, &testy, &ntest);
        release_data(X, y, NULL);
        scale_data(trainX, ntrain, p, testX, ntest, true, &mu, &s);

        nboot       = atoi(argv[arg_offset + 6]);
        pop_size    = atoi(argv[arg_offset + 7]);
        ngen        = atoi(argv[arg_offset + 8]);

        pcross      = 0.5;
        pmutation   = 0.5;

        for (r = 0; r < R; ++r) {
            ens = sspbe_evolve((const double **)trainX, trainy,
                               ntrain, p,
                               (const double **)testX, testy, ntest,
                               pop_size, ngen, nboot, 1, 2,
                               pcross, pmutation,
                               2, 6, 4, -1, 500,
                               NULL, NULL, -1, NULL, -1,
                               verbose,
                               next_rnd);

            for (i = 0; i < ntest; ++i) fprintf(stdout, "%2d %2d %3d %4d %g %g\n", D, r + 1, nboot, i, testy[i], ensemble_predict(ens, testX[i]));
            fflush(stdout);
            ensemble_free(ens);
        }

        free(s);
        free(mu);
        release_data(testX, testy, NULL);
        release_data(trainX, trainy, NULL);
    } else {
        if (access(argv[arg_offset + 1], F_OK ) != -1) {
            problem = argv[arg_offset + 1];
            read_data(argv[arg_offset + 1], &X, &y, NULL, &n, &p);
            leave_in = sample_indices(n, atof(argv[arg_offset + 2]), next_rnd);
            split_data(X, y, leave_in, n, p,
                    &trainX, &trainy, &ntrain,
                    &testX, &testy, &ntest);
            free(leave_in);
            release_data(X, y, NULL);
            arg_offset += 2;
        } else if (argv[arg_offset + 1][0] == 'L') {
            problem = argv[arg_offset + 2];
            read_data(argv[arg_offset + 2], &X, &y, NULL, &n, &p);
            configure_data(argv[arg_offset + 3], atoi(argv[arg_offset + 4]),
                           X, y, n, p,
                           &trainX, &trainy, &ntrain,
                           &testX, &testy, &ntest);
            release_data(X, y, NULL);
            arg_offset += 4;
        } else {
            fprintf(stderr, "%s:%d - ERROR: Failed to find data file.\n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        scale_data(trainX, ntrain, p, testX, ntest, true, NULL, NULL);

        nboot       = atoi(argv[arg_offset + 1]);
        pop_size    = atoi(argv[arg_offset + 2]);
        ngen        = atoi(argv[arg_offset + 3]);
        if ((argc - arg_offset) > 4) verbose = argv[arg_offset + 4][0] == 'Y' || argv[arg_offset + 4][0] == 'y';

        pcross      = 0.5;
        pmutation   = 0.5;

        ens = sspbe_evolve((const double **)trainX, trainy,
                           ntrain, p,
                           (const double **)testX, testy, ntest,
                           pop_size, ngen, nboot, 1, 2,
                           pcross, pmutation,
                           2, 6, 4, -1, 500,
                           NULL, NULL, -1, NULL, -1,
                           verbose,
                           next_rnd);

        if (!verbose) {
            fprintf(stdout, "%s %lu %3d %3d %3d %g %g %d %d\n",
                    problem, seed, nboot, pop_size, ngen,
                    sqrt(ensemble_mse(ens, (const double * const *)trainX, trainy, ntrain)),
                    sqrt(ensemble_mse(ens, (const double * const *)testX, testy, ntest)),
                    ensemble_total_size(ens),
                    ensemble_unique_size(ens));
        }

        ensemble_free(ens);

        release_data(testX, testy, NULL);
        release_data(trainX, trainy, NULL);
    }

    return EXIT_SUCCESS;
}
