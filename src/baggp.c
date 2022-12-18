#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <string.h>

#include "alloc.h"
#include "baselearner.h"
#include "parsetree.h"
#include "util.h"









struct model {
    struct parsetree *m;

    double a;
    double b;
    double f;
};










static struct model **create_population(int N)
{
    int i;
    struct model **pop = MALLOC(N, sizeof(struct model *));

    for (i = 0; i < N; ++i) {
        pop[i] = MALLOC(1, sizeof(struct model));
        pop[i]->a = NAN;
        pop[i]->b = NAN;
        pop[i]->f = INFINITY;
        pop[i]->m = NULL;
    }

    return pop;
}





static void release_population(struct model **pop, int N)
{
    int i;

    for (i = 0; i < N; ++i) {
        parsetree_release(pop[i]->m);
        free(pop[i]);
    }

    free(pop);
}




static void replace_model(struct model *dest, const struct model *src)
{
    parsetree_release(dest->m);
    dest->m = parsetree_copy(src->m);
    dest->a = src->a;
    dest->b = src->b;
    dest->f = src->f;
}





static void evaluate_model(struct model *ind, const double * const *X, const double *t, double *z, int n)
{
    int j;
    for (j = 0; j < n; ++j) z[j] = parsetree_predict(ind->m, X[j]);
    ind->f = fit(t, z, n, NULL, &(ind->a), &(ind->b));
}





static double fitness_model(struct model *ind, const double * const *X, const double *t, double *z, int n)
{
    int j;
    for (j = 0; j < n; ++j) z[j] = ind->a + ind->b * parsetree_predict(ind->m, X[j]);
    return mse(t, z, n);
}





static int initial_generation(struct model **pop, int N,
                              const double * const *trainX, const double *traint, int ntrain,
                              int min_init_height, int max_init_height,
                              const parsetree_nodeop *functions, const int *arity, int nfunc,
                              const parsetree_nodeop *terminals, int nterm,
                              double (*rnd)(void))
{
    int best = 0, idx, i, j, k, half_pop_size = N / 2;
    int attempt, max_attempts = 10 * N;
    int curr_max_height, next_height_interval, init_height_interval = N / (max_init_height - min_init_height + 1) / 2;
    bool grow, successful;
    double *z = MALLOC(ntrain, sizeof(double));

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
                parsetree_release(pop[idx]->m);
                pop[idx]->m = parsetree_create(0, curr_max_height, grow, functions, arity, nfunc, terminals, nterm, rnd);
                successful = true;
                for (j = 0; (j < idx) && successful; ++j) successful = !parsetree_equal(pop[idx]->m, pop[j]->m);
                if (successful) break;
            }
            evaluate_model(pop[idx], trainX, traint, z, ntrain);
            if (cmp_double(pop[idx]->f, pop[best]->f) <= 0) best = idx;
        }
    }

    free(z);

    return best;
}



static int tournament(struct model **pop, int N, int K, double (*rnd)(void), int last)
{
    int pick, i = 1, j = 0, *samples = sample_m_of_n(K + 1, N, rnd);

    do { pick = samples[j++]; } while (pick == last);

    while (i < K) {
        if (samples[j] != last) {
            i++;
            if (cmp_double(pop[samples[j]]->f, pop[pick]->f) <= 0) {
                pick = samples[j];
            }
        }
        j++;
    }

    free(samples);

    return pick;
}



static int breed_generation(struct model **pop, struct model **gen, int N, int K,
                            const double * const *trainX, const double *traint, int ntrain,
                            double pm, double pc,
                            int max_mut_depth, int max_tree_height, int max_tree_size,
                            const parsetree_nodeop *functions, const int *arity, int nfunc,
                            const parsetree_nodeop *terminals, int nterm,
                            double (*rnd)(void))
{
    int best = 0, i, m, f;
    double p, *z = MALLOC(ntrain, sizeof(double));
    bool eval;

    for (i = 1; i < N; ++i) {
        m = tournament(pop, N, K, rnd, -1);

        eval = false;
        p = rnd();
        if (p < (pm + pc)) {
            parsetree_release(gen[i]->m);
            if (p < pm) {
                gen[i]->m = parsetree_mutation(parsetree_copy(pop[m]->m), max_mut_depth, functions, arity, nfunc, terminals, nterm, rnd);
            } else {
                f = tournament(pop, N, K, rnd, m);
                gen[i]->m = parsetree_crossover(parsetree_copy(pop[m]->m), pop[f]->m, rnd);
            }
            eval = (max_tree_size < 0 || parsetree_size(gen[i]->m) <= max_tree_size) && (max_tree_height < 0 || parsetree_height(gen[i]->m) <= max_tree_height);
        }

        if (eval) {
            evaluate_model(gen[i], trainX, traint, z, ntrain);
        } else {
            replace_model(gen[i], pop[m]);
        }
        if (cmp_double(gen[i]->f, gen[best]->f) <= 0) best = i;
    }

    free(z);

    return best;
}










struct baselearner *gp_evolve(
        const double * const *trainX, const double *traint, int ntrain, int nvar,
        const double * const *testX, const double *testt, int ntest,
        int N, int G, int K,
        double pc, double pm,
        int min_init_depth, int max_init_depth,
        int max_mut_depth,
        int max_tree_height, int max_tree_size,
        const parsetree_nodeop * functions, const int * arity, int nfunc,
        const parsetree_nodeop * terminals, int nterm,
        bool verbose,
        double (*rnd)(void))
{
    static const parsetree_nodeop def_functions[] = { exec_add, exec_sub, exec_mul };
    static const int def_arity[] = { 2, 2, 2 };
    static const int def_nfunc = 3;

    static const parsetree_nodeop def_terminals[] = { exec_erc };
    static const int def_nterm = 1;

    int g = 0, best = 0, i;
    struct baselearner *M = baselearner_create();
    struct model **pop, **gen, **tmp;

    double *testz = MALLOC(ntest, sizeof(double));

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

    pop = create_population(N);
    gen = create_population(N);

    best = initial_generation(pop, N,
                              trainX, traint, ntrain,
                              min_init_depth, max_init_depth,
                              int_functions, int_arity, nfunc, int_terminals, nterm,
                              rnd);
    if (verbose) fprintf(stderr, "%3d %g %g\n", g, pop[best]->f, fitness_model(pop[best], testX, testt, testz, ntest));
    for (g = 1; g <= G; ++g) {
        replace_model(gen[0], pop[best]);

        best = breed_generation(pop, gen, N, K,
                                trainX, traint, ntrain,
                                pm, pc, max_mut_depth, max_tree_height, max_tree_size,
                                int_functions, int_arity, nfunc, int_terminals, nterm,
                                rnd);

        tmp = pop;
        pop = gen;
        gen = tmp;

        if (verbose)fprintf(stderr, "%3d %g %g\n", g, pop[best]->f, fitness_model(pop[best], testX, testt, testz, ntest));
    }

    baselearner_set_model(M, parsetree_copy(pop[best]->m), pop[best]->a, pop[best]->b);

    release_population(gen, N);
    release_population(pop, N);

    free(var_terminals);

    free(testz);

    return M;
}
