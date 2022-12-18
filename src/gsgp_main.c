#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <unistd.h>
#include <string.h>
#include <math.h>

#include "data.h"
#include "gsgp.h"
#include "readline.h"
#include "util.h"

#include "rng.h"


struct decomp_args
{
    bool include_initial_gen;
    int gen_freq;
    int D;
    int R;
    int mutation_depth;
    double mutation_scale;
};

static void report(int g,
                   const double *traint, double *trainy, int ntrain, double var_train,
                   const double *testt, double *testy, int ntest, double var_test,
                   __attribute__((unused)) void *args_ptr)
{
    double train = mse(traint, trainy, ntrain);
    double test  = mse(testt, testy, ntest);

    fprintf(stdout, "%3d %g %g   |   %g %g\n", g, train, test, sqrt(train / var_train), sqrt(test / var_test));
    fflush(stdout);
}

static void decomp(int g,
                   __attribute__((unused)) const double *traint, __attribute__((unused)) double *trainy, __attribute__((unused)) int ntrain, __attribute__((unused)) double var_train,
                   const double *testt, double *testy, int ntest, __attribute__((unused)) double var_test,
                   void *args_ptr)
{
    int i;
    struct decomp_args *args = args_ptr;

    if (g == 0 && !args->include_initial_gen) return;
    if ((g % args->gen_freq) != 0) return;
    for (i = 0; i < ntest; ++i) {
        fprintf(stdout, "%2d %2d %g %3d %4d %4d %g %g\n",
                args->D, args->R, args->mutation_scale, args->mutation_depth, g, i, testt[i], testy[i]);
    }
    fflush(stdout);
}

int main(int argc __attribute__((unused)), char **argv)
{
    char *problem;

    double **X, *y;
    double **trainX, *traint, *mu = NULL, *s = NULL;
    double **testX, *testt, *testy;
    bool *leave_in;
    int n, p, ntrain, ntest;
    bool verbose = false;
    int r, D, R, arg_offset;

    int pop_size, ngen, mutation_depth;
    double pcross, pmutation, mutation_scale;

    struct decomp_args dargs;

    unsigned long seed;
    seed_rng(argv[1], atoi(argv[2]));
    seed = atoi(argv[2]);
    arg_offset = 2;

    dargs.gen_freq = 10;
    dargs.include_initial_gen = true;

    if (argv[arg_offset + 1][0] == 'D') {
        D = atoi(argv[arg_offset + 4]);
        R = atoi(argv[arg_offset + 5]);

        read_data(argv[arg_offset + 2], &X, &y, NULL, &n, &p);
        configure_data(argv[arg_offset + 3], D,
                       X, y, n, p,
                       &trainX, &traint, &ntrain,
                       &testX, &testt, &ntest);
        release_data(X, y, NULL);
        arg_offset += 5;

        mutation_scale = atof(argv[arg_offset + 1]);
        mutation_depth = atoi(argv[arg_offset + 2]);
        pop_size       = atoi(argv[arg_offset + 3]);
        ngen           = atoi(argv[arg_offset + 4]);
        if ((argc - arg_offset) > 5) dargs.gen_freq =  atoi(argv[arg_offset + 5]);
        if ((argc - arg_offset) > 6) dargs.include_initial_gen = argv[arg_offset + 6][0] == 'Y' || argv[arg_offset + 6][0] == 'y';

        pcross      = 0.3;
        pmutation   = 0.7;

        dargs.D = D;
        dargs.mutation_depth = mutation_depth;
        dargs.mutation_scale = mutation_scale;

        for (r = 0; r < R; ++r) {
            dargs.R = r + 1;
            testy = gsgp_evolve((const double **)trainX, traint,
                                ntrain, p,
                                (const double **)testX, testt, ntest,
                                pop_size, ngen, 6,
                                pcross, pmutation,
                                2, 6, mutation_depth, mutation_scale,
                                NULL, NULL, -1, NULL, -1,
                                decomp,
                                next_rnd, &dargs);

            free(testy);
        }

        free(s);
        free(mu);
        release_data(testX, testt, NULL);
        release_data(trainX, traint, NULL);
    } else {
        if (access(argv[arg_offset + 1], F_OK ) != -1) {
            problem = argv[arg_offset + 1];
            read_data(argv[arg_offset + 1], &X, &y, NULL, &n, &p);
            leave_in = sample_indices(n, atof(argv[arg_offset + 2]), next_rnd);
            split_data(X, y, leave_in, n, p,
                       &trainX, &traint, &ntrain,
                       &testX, &testt, &ntest);
            free(leave_in);
            release_data(X, y, NULL);
            arg_offset += 2;
        } else if (argv[arg_offset + 1][0] == 'L') {
            problem = argv[arg_offset + 2];
            read_data(argv[arg_offset + 2], &X, &y, NULL, &n, &p);
            configure_data(argv[arg_offset + 3], atoi(argv[arg_offset + 4]),
                           X, y, n, p,
                           &trainX, &traint, &ntrain,
                           &testX, &testt, &ntest);
            release_data(X, y, NULL);
            arg_offset += 4;
        } else {
            fprintf(stderr, "%s:%d - ERROR: Failed to find data file.\n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        mutation_scale = atof(argv[arg_offset + 1]);
        mutation_depth = atoi(argv[arg_offset + 2]);
        pop_size    = atoi(argv[arg_offset + 3]);
        ngen        = atoi(argv[arg_offset + 4]);
        if ((argc - arg_offset) > 5) verbose = argv[arg_offset + 5][0] == 'Y' || argv[arg_offset + 5][0] == 'y';

        pcross      = 0.3;
        pmutation   = 0.7;

        testy = gsgp_evolve((const double **)trainX, traint,
                            ntrain, p,
                            (const double **)testX, testt, ntest,
                            pop_size, ngen, 6,
                            pcross, pmutation,
                            2, 6, mutation_depth, mutation_scale,
                            NULL, NULL, -1, NULL, -1,
                            verbose ? report : NULL,
                            next_rnd, NULL);

        if (!verbose) {
            fprintf(stdout, "%s %lu %3d %f %3d %3d %g\n",
                    problem, seed, mutation_depth, mutation_scale, pop_size, ngen, mse(testt, testy, ntest));
        }

        free(testy);
        release_data(testX, testt , NULL);
        release_data(trainX, traint, NULL);
    }

    return EXIT_SUCCESS;
}
