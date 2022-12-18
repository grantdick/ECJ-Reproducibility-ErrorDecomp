#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <string.h>

#include "alloc.h"
#include "baselearner.h"
#include "sspbe_ensemble.h"
#include "parsetree.h"
#include "util.h"









struct location {
    int bootstrap_id;
    const int *B;

    struct baselearner *current;
    struct baselearner *proposed;

    double *cz;
    double *pz;

    double cf;
    double pf;
};










static int map_location(int x, int side)
{
    while (x < 0) x += side;
    while (x >= side) x -= side;
    return x;
}





static double mean_tree_size(struct location ** pop, int rows, int cols)
{
    int r, c, n;
    double m;

    m = 0;
    n = 0;
    for (r = 0; r < rows; ++r) {
        for (c = 0; c < cols; ++c) {
            m += (parsetree_size(baselearner_get_model(pop[r][c].current)) - m) / ++n;
        }
    }

    return m;
}





static const struct location *tournament(struct location **pop,
                                         int r, int c, int rows, int cols,
                                         int dim, int K,
                                         const double *t, int n,
                                         double (*rnd)(void),
                                         const struct location *last)
{
    int *samples;
    int i, j, nr, nc, d, pick;
    double pick_fitness, challenge_fitness;

    d = 2 * dim + 1;
    samples = sample_m_of_n(K + 1, d * d, rnd);

    j = 0;
    do {
        pick = samples[j++];
        nr = map_location(r + pick / d - dim, rows);
        nc = map_location(c + pick % d - dim, cols);
    } while (&(pop[nr][nc]) == last);
    pick_fitness = fit(t, pop[nr][nc].cz, n, pop[r][c].B, NULL, NULL);

    i = 1;
    while (i < K) {
        nr = map_location(r + samples[j] / d - dim, rows);
        nc = map_location(c + samples[j] % d - dim, cols);
        if (&(pop[nr][nc]) != last) {
            i++;
            challenge_fitness = fit(t, pop[nr][nc].cz, n, pop[r][c].B, NULL, NULL);
            if (cmp_double(challenge_fitness, pick_fitness) <= 0) {
                pick_fitness = challenge_fitness;
                pick = samples[j];
            }
        }
        j++;
    }

    free(samples);

    nr = map_location(r + pick / d - dim, rows);
    nc = map_location(c + pick % d - dim, cols);
    return &(pop[nr][nc]);
}





static struct location **create_population(int rows, int cols,
                                           const struct ensemble *E, int ntrain,
                                           double (*rnd)(void))
{
    int i, j, r, c, *bootid;
    struct location **pop;

    bootid = MALLOC(rows * cols, sizeof(int));
    for (i = 0; i < (rows * cols); ++i) {
        j = (int)(rnd() * (i + 1));
        if (i != j) bootid[i] = bootid[j];
        bootid[j] = i % ensemble_size(E);
    }

    pop = MALLOC(rows, sizeof(struct location *));
    pop[0] = MALLOC(rows * cols, sizeof(struct location));

    for (r = 0; r < rows; ++r) {
        pop[r] = pop[0] + r * cols;
        for (c = 0; c < cols; ++c) {
            pop[r][c].bootstrap_id = bootid[r * cols + c];
            pop[r][c].B = ensemble_bootstrap(E, pop[r][c].bootstrap_id);

            pop[r][c].current = baselearner_create();
            pop[r][c].proposed = baselearner_create();

            pop[r][c].cz = MALLOC(ntrain, sizeof(double));
            pop[r][c].pz = MALLOC(ntrain, sizeof(double));

            pop[r][c].cf = NAN;
            pop[r][c].pf = NAN;
        }
    }

    free(bootid);

    return pop;
}





static void release_population(struct location **pop, int rows, int cols)
{
    int r, c;

    for (r = 0; r < rows; ++r) {
        for (c = 0; c < cols; ++c) {
            baselearner_release(pop[r][c].proposed);
            baselearner_release(pop[r][c].current);
            free(pop[r][c].pz);
            free(pop[r][c].cz);
        }
    }

    free(pop[0]);
    free(pop);
}




static void initial_generation(struct location **pop, int rows, int cols,
                               struct ensemble *E,
                               const double **X, const double *t, int n,
                               int min_init_height, int max_init_height,
                               const parsetree_nodeop * const functions, const int * const arity, int nfunc,
                               const parsetree_nodeop * const terminals, int nterm,
                               double (*rnd)(void))
{
    int idx, r, c, i, j, k, N = rows * cols, half_pop_size = N / 2, reqd;
    int attempt, max_attempts = 10 * N;
    int curr_max_height, next_height_interval, init_height_interval = N / (max_init_height - min_init_height + 1) / 2;
    bool grow, successful;

    struct parsetree *model, **models = MALLOC(N, sizeof(struct parsetree *));
    int *shuffle = MALLOC(N, sizeof(int));
    double a, b;

    for (k = 0; k < 2; ++k) {
        grow = k != 0;
        curr_max_height = min_init_height;
        next_height_interval = init_height_interval;

        reqd = half_pop_size + k * (N % 2);
        for (i = 0; i < reqd; ++i) {
            if (i >= next_height_interval) {
                next_height_interval += init_height_interval;
                curr_max_height++;
            }

            idx = i + (k * half_pop_size);
            models[idx] = NULL;
            for (attempt = 0; attempt < max_attempts; ++attempt) {
                parsetree_release(models[idx]);
                models[idx] = parsetree_create(0, curr_max_height, grow, functions, arity, nfunc, terminals, nterm, rnd);
                successful = true;
                for (j = 0; (j < idx) && successful; ++j) successful = !parsetree_equal(models[idx], models[j]);
                if (successful) break;
            }
        }
    }

    for (i = 0; i < N; ++i) {
        j = (int)(rnd() * (i + 1));
        if (i != j) shuffle[i] = shuffle[j];
        shuffle[j] = i;
    }

    for (i = 0; i < (rows * cols); ++i) {
        r = shuffle[i] / cols;
        c = shuffle[i] % cols;

        model = models[shuffle[i]];
        for (j = 0; j < n; ++j) pop[r][c].cz[j] = parsetree_predict(model, X[j]);
        pop[r][c].cf = fit(t, pop[r][c].cz, n, pop[r][c].B, &a, &b);
        baselearner_set_model(pop[r][c].current, model, a, b);
        ensemble_update(E, pop[r][c].bootstrap_id, pop[r][c].current, pop[r][c].cz);
    }

    free(shuffle);
    free(models);
}





static void breed_generation(struct location **pop, int rows, int cols,
                             int deme_size, int tourn_size,
                             const double **X, const double *t, int n,
                             double pm, double pc,
                             int max_mut_depth, int max_tree_height, int max_tree_size,
                             const parsetree_nodeop * const functions, const int * const arity, int nfunc,
                             const parsetree_nodeop * const terminals, int nterm,
                             double (*rnd)(void))
{
    int j, r, c;
    double p, a, b;
    const struct location *m, *f;
    struct parsetree *offspring;

    for (r = 0; r < rows; ++r) {
        for (c = 0; c < cols; ++c) {
            m = tournament(pop, r, c, rows, cols, deme_size, tourn_size, t, n, rnd, NULL);

            p = rnd();
            if (p < pm) {
                offspring = parsetree_mutation(parsetree_copy(baselearner_get_model(m->current)),
                                               max_mut_depth, functions, arity, nfunc, terminals, nterm, rnd);
            } else if (p < (pm + pc)) {
                f = tournament(pop, r, c, rows, cols, deme_size, tourn_size, t, n, rnd, m);
                offspring = parsetree_crossover(parsetree_copy(baselearner_get_model(m->current)), baselearner_get_model(f->current), rnd);
            } else {
                offspring = parsetree_copy(baselearner_get_model(m->current));
            }

            if ((max_tree_size >= 0 && parsetree_size(offspring) > max_tree_size) || (max_tree_height >= 0 && parsetree_height(offspring) > max_tree_height)) {
                parsetree_release(offspring);
                offspring = parsetree_copy(baselearner_get_model(pop[r][c].current));
            }

            for (j = 0; j < n; ++j) pop[r][c].pz[j] = parsetree_predict(offspring, X[j]);
            pop[r][c].pf = fit(t, pop[r][c].pz, n, pop[r][c].B, &a, &b);
            baselearner_set_model(pop[r][c].proposed, offspring, a, b);
        }
    }
}





static void elitism(struct location **pop, int rows, int cols,
                    struct ensemble *E, int n)
{
    int r, c;
    int csize, psize, cmp;
    bool replace;

    for (r = 0; r < rows; ++r) {
        for (c = 0; c < cols; ++c) {
            psize = parsetree_size(baselearner_get_model(pop[r][c].proposed));
            csize = parsetree_size(baselearner_get_model(pop[r][c].current));

            cmp = cmp_double(pop[r][c].pf, pop[r][c].cf);
            replace = (cmp < 0) || (cmp == 0 && psize < csize);
            if (replace) {
                baselearner_release(pop[r][c].current);
                pop[r][c].current = baselearner_copy(pop[r][c].proposed);
                memcpy(pop[r][c].cz, pop[r][c].pz, n * sizeof(double));
                pop[r][c].cf = pop[r][c].pf;

                ensemble_update(E, pop[r][c].bootstrap_id, pop[r][c].current, pop[r][c].cz);
            }
        }
    }
}





static void report(int g,
                   const struct ensemble *E,
                   struct location **pop, int rows, int cols,
                   double msr_train, double msr_test)
{
    double ib, oob, test;

    ib   = sqrt(ensemble_ib_error(E));
    oob  = sqrt(ensemble_oob_error(E));
    test = sqrt(ensemble_test_error(E));

    fprintf(stdout, "%3d %f %f %f %f %f %f %d\n",
            g, ib, oob, test, oob / sqrt(msr_train), test / sqrt(msr_test),
            mean_tree_size(pop, rows, cols), ensemble_unique_size(E)
        );
    fflush(stdout);
}





struct ensemble *sspbe_evolve(const double **trainX, const double *traint, int ntrain, int nvar,
                              const double **testX, const double *testt, int ntest,
                              int N, int G, int nboot, int deme_size, int tourn_size,
                              double pc, double pm,
                              int min_init_depth, int max_init_depth,
                              int max_mut_depth,
                              int max_tree_height, int max_tree_size,
                              const parsetree_nodeop * const functions, const int * const arity, int nfunc,
                              const parsetree_nodeop * const terminals, int nterm,
                              bool verbose,
                              double (*rnd)(void))
{
    static const parsetree_nodeop def_functions[] = { exec_add, exec_sub, exec_mul };
    static const int def_arity[] = { 2, 2, 2 };
    static const int def_nfunc = 3;

    static const parsetree_nodeop def_terminals[] = { exec_erc };
    static const int def_nterm = 1;

    int g = 0, i;

    struct ensemble *E;
    struct location **pop;
    int rows = (int)sqrt(N), cols = N / rows;
    double msr_train = var(traint, ntrain, false), msr_test = var(testt, ntest, false);

    const parsetree_nodeop * int_functions, *int_terminals;
    parsetree_nodeop *var_terminals;
    const int *int_arity;

    if (functions == NULL) {
        int_functions = def_functions;
        int_arity = def_arity;
        nfunc = def_nfunc;
    } else {
        int_functions = functions;
        int_arity = arity;
    }

    if (terminals == NULL) {
        nterm = nvar + def_nterm;
        var_terminals = MALLOC(nterm, sizeof(parsetree_nodeop));
        for (i = 0; i < nvar; ++i) var_terminals[i] = exec_var;
        for (i = 0; i < def_nterm; ++i) var_terminals[i + nvar] = def_terminals[i];
    } else {
        var_terminals = MALLOC(nvar + nterm, sizeof(parsetree_nodeop));
        for (i = 0; i < nvar; ++i) var_terminals[i] = exec_var;
        for (i = 0; i < nterm; ++i) var_terminals[i + nvar] = terminals[i];
        nterm += nvar;
    }
    int_terminals = var_terminals;

    if (nboot == 0) {
        E = ensemble_create(1, traint, ntrain, testX, testt, ntest);
        ensemble_makeboot(E, NULL);
    } else if (nboot < 1) {
        E = ensemble_create(-nboot, traint, ntrain, testX, testt, ntest);
        ensemble_makeboot(E, NULL);
    } else {
        E = ensemble_create(nboot, traint, ntrain, testX, testt, ntest);
        ensemble_makeboot(E, rnd);
    }

    N = rows * cols;

    pop = create_population(rows, cols, E, ntrain, rnd);

    initial_generation(pop, rows, cols,
                       E,
                       trainX, traint, ntrain,
                       min_init_depth, max_init_depth,
                       int_functions, int_arity, nfunc, int_terminals, nterm,
                       rnd);

    if (verbose) report(g, E, pop, rows, cols, msr_train, msr_test);

    for (g = 1; g <= G; ++g) {
        breed_generation(pop, rows, cols, deme_size, tourn_size,
                         trainX, traint, ntrain,
                         pm, pc, max_mut_depth, max_tree_height, max_tree_size,
                         int_functions, int_arity, nfunc, int_terminals, nterm,
                         rnd);
        elitism(pop, rows, cols, E, ntrain);
        if (verbose) report(g, E, pop, rows, cols, msr_train, msr_test);
    }

    release_population(pop, rows, cols);
    free(var_terminals);

    return E;
}
