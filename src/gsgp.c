#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <string.h>

#include "alloc.h"
#include "parsetree.h"
#include "util.h"









struct individual {
    double *trainy;
    double *testy;

    double fitness;
    double test_mse;
};










static const struct individual *tournament(struct individual *pop, int N, int K,
                                           double (*rnd)(void),
                                           const struct individual *last)
{
    int *samples;
    int i, j, pick;

    samples = sample_m_of_n(K + 1, N, rnd);

    j = 0;
    do { pick = samples[j++]; } while (&(pop[pick]) == last);

    i = 1;
    while (i < K) {
        if (&(pop[samples[j]]) != last) {
            i++;
            if (cmp_double(pop[samples[j]].fitness, pop[pick].fitness) <= 0) {
                pick = samples[j];
            }
        }
        j++;
    }

    free(samples);

    return &(pop[pick]);
}





static struct individual *create_population(int N, int ntrain, int ntest)
{
    struct individual *pop = MALLOC(N, sizeof(struct individual));
    int i;

    for (i = 0; i < N; ++i) {
        pop[i].trainy = MALLOC(ntrain, sizeof(double));
        pop[i].testy = MALLOC(ntest, sizeof(double));
    }

    return pop;
}





static void release_population(struct individual *pop, int N)
{
    int i;

    for (i = 0; i < N; ++i) {
        free(pop[i].testy);
        free(pop[i].trainy);
    }

    free(pop);
}




static int initial_generation(struct individual *pop, int N,
                              const double **trainX, const double *traint, int ntrain,
                              const double **testX, const double *testt, int ntest,
                              int min_init_height, int max_init_height,
                              const parsetree_nodeop * const functions, const int * const arity, int nfunc,
                              const parsetree_nodeop * const terminals, int nterm,
                              double (*rnd)(void))
{
    int best, idx, i, j, k, half_pop_size = N / 2;
    int attempt, max_attempts = 10 * N;
    int curr_max_height, next_height_interval, init_height_interval = N / (max_init_height - min_init_height + 1) / 2;
    bool grow, successful;
    double *z = MALLOC(ntrain, sizeof(double));

    struct parsetree **models = MALLOC(N, sizeof(struct parsetree *));
    for (i = 0; i < N; ++i) models[i] = NULL;

    for (k = 0; k < 2; ++k) {
        grow = k != 0;
        curr_max_height = min_init_height;
        next_height_interval = init_height_interval;

        for (i = 0; i < half_pop_size; ++i) {
            if (i >= next_height_interval) {
                next_height_interval += init_height_interval;
                curr_max_height++;
            }

            idx = i + (k * half_pop_size);
            for (attempt = 0; attempt < max_attempts; ++attempt) {
                parsetree_release(models[idx]);
                models[idx] = parsetree_create(0, curr_max_height, grow, functions, arity, nfunc, terminals, nterm, rnd);
                successful = true;
                for (j = 0; (j < idx) && successful; ++j) successful = !parsetree_equal(models[idx], models[j]);
                if (successful) break;
            }

            for (j = 0; j < ntrain; ++j) pop[idx].trainy[j] = parsetree_predict(models[idx], trainX[j]);
            for (j = 0; j < ntest; ++j) pop[idx].testy[j] = parsetree_predict(models[idx], testX[j]);

            pop[idx].fitness = mse(traint, pop[idx].trainy, ntrain);
            pop[idx].test_mse = mse(testt, pop[idx].testy, ntest);
        }
    }

    best = 0;
    for (i = 0; i < N; ++i) {
        if (cmp_double(pop[i].fitness, pop[best].fitness) < 0) best = i;
        parsetree_release(models[i]);
    }
    free(models);
    free(z);

    return best;
}





static int breed_generation(struct individual *pop, struct individual *gen,
                            int N, int tourn_size,
                            const double **trainX, const double *traint, int ntrain,
                            const double **testX, const double *testt, int ntest,
                            double pm, double pc,
                            int max_mut_depth, double mutation_scale,
                            const parsetree_nodeop * const functions, const int * const arity, int nfunc,
                            const parsetree_nodeop * const terminals, int nterm,
                            double (*rnd)(void))
{
    int i, j, best = 0;
    double p, alpha;
    const struct individual *m, *f;
    struct parsetree *r1, *r2;

    for (i = 0; i < N; ++i) {
        m = tournament(pop, N, tourn_size, rnd, NULL);

        p = rnd();
        if (p < pm) {
            r1 = parsetree_create(0, max_mut_depth, true, functions, arity, nfunc, terminals, nterm, rnd);
            r2 = parsetree_create(0, max_mut_depth, true, functions, arity, nfunc, terminals, nterm, rnd);

            for (j = 0; j < ntrain; ++j) gen[i].trainy[j] = m->trainy[j] + mutation_scale * (parsetree_predict(r1, trainX[j]) - parsetree_predict(r2, trainX[j]));
            for (j = 0; j < ntest; ++j) gen[i].testy[j] = m->testy[j] + mutation_scale * (parsetree_predict(r1, testX[j]) - parsetree_predict(r2, testX[j]));

            parsetree_release(r2);
            parsetree_release(r1);
        } else if (p < (pm + pc)) {
            f = tournament(pop, N, tourn_size, rnd, m);

            alpha = rnd();
            for (j = 0; j < ntrain; ++j) gen[i].trainy[j] = alpha * m->trainy[j] + (1 - alpha) * f->trainy[j];
            for (j = 0; j < ntest; ++j) gen[i].testy[j] = alpha * m->testy[j] + (1 - alpha) * f->testy[j];
        } else {
            memcpy(gen[i].trainy, m->trainy, ntrain * sizeof(double));
            memcpy(gen[i].testy, m->testy, ntest * sizeof(double));
        }

        gen[i].fitness = mse(traint, gen[i].trainy, ntrain);
        gen[i].test_mse = mse(testt, gen[i].testy, ntest);

        if (cmp_double(gen[i].fitness, gen[best].fitness) < 0) best = i;
    }

    return best;
}





/* static void report(int g, struct individual *pop, __attribute__((unused)) int N, int best, */
/*                    double msr_train, double msr_test) */
/* { */
/*     double train, test; */

/*     train = pop[best].fitness; */
/*     test  = pop[best].test_mse; */

/*     fprintf(stdout, "%3d %g %g   |   %g %g\n", g, train, test, sqrt(train / msr_train), sqrt(test / msr_test)); */
/*     fflush(stdout); */
/* } */





double *gsgp_evolve(const double **trainX, const double *traint, int ntrain, int nvar,
                    const double **testX, const double *testt, int ntest,
                    int N, int G, int tourn_size,
                    double pc, double pm,
                    int min_init_depth, int max_init_depth,
                    int max_mut_depth, double mutation_scale,
                    const parsetree_nodeop * const functions, const int * const arity, int nfunc,
                    const parsetree_nodeop * const terminals, int nterm,
                    void (*report)(int, const double *, double *, int, double, const double *, double *, int, double, void *),
                    double (*rnd)(void), void *report_args)
{
    static const parsetree_nodeop def_functions[] = { exec_add, exec_sub, exec_mul };
    static const int def_arity[] = { 2, 2, 2 };
    static const int def_nfunc = 3;

    static const parsetree_nodeop def_terminals[] = { exec_erc };
    static const int def_nterm = 1;

    double msr_train = var(traint, ntrain, false), msr_test = var(testt, ntest, false);

    int best, g = 0, i;

    const parsetree_nodeop * int_functions, *int_terminals;
    parsetree_nodeop *var_terminals;
    const int *int_arity;

    struct individual *pop = create_population(N, ntrain, ntest), *gen = create_population(N, ntrain, ntest), *tmp;

    double *ret = MALLOC(ntest, sizeof(double));

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

    best = initial_generation(pop, N,
                              trainX, traint, ntrain,
                              testX, testt, ntest,
                              min_init_depth, max_init_depth,
                              int_functions, int_arity, nfunc, int_terminals, nterm,
                              rnd);

    if (report) report(g, traint, pop[best].trainy, ntrain, msr_train, testt, pop[best].testy, ntest, msr_test, report_args);

    for (g = 1; g <= G; ++g) {
        best = breed_generation(pop, gen, N, tourn_size,
                                trainX, traint, ntrain,
                                testX, testt, ntest,
                                pm, pc,
                                max_mut_depth, mutation_scale,
                                int_functions, int_arity, nfunc, int_terminals, nterm,
                                rnd);

        tmp = pop;
        pop = gen;
        gen = tmp;

        if (report) report(g, traint, pop[best].trainy, ntrain, msr_train, testt, pop[best].testy, ntest, msr_test, report_args);
    }

    memcpy(ret, pop[best].testy, ntest * sizeof(double));

    release_population(gen, N);
    release_population(pop, N);

    free(var_terminals);

    return ret;
}
