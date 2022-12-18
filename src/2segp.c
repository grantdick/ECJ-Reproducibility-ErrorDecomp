#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <string.h>

#include "alloc.h"
#include "2segp_ensemble.h"
#include "parsetree.h"
#include "util.h"









struct model {
    struct parsetree *m;

    double *a;
    double *b;
    double *f;
};










static struct model **create_population(int N, int nboot)
{
    int b, i;
    struct model **pop = MALLOC(N, sizeof(struct model *));

    for (i = 0; i < N; ++i) {
        pop[i] = MALLOC(1, sizeof(struct model));
        pop[i]->a = MALLOC(nboot, sizeof(double));
        pop[i]->b = MALLOC(nboot, sizeof(double));
        pop[i]->f = MALLOC(nboot, sizeof(double));
        pop[i]->m = NULL;

        for (b = 0; b < nboot; ++b) {
            pop[i]->a[b] = NAN;
            pop[i]->b[b] = NAN;
            pop[i]->f[b] = INFINITY;
        }
    }

    return pop;
}





static void release_population(struct model **pop, int N)
{
    int i;

    for (i = 0; i < N; ++i) {
        parsetree_release(pop[i]->m);
        free(pop[i]->f);
        free(pop[i]->b);
        free(pop[i]->a);
        free(pop[i]);
    }

    free(pop);
}




static void replace_model(struct model *dest, const struct model *src, int b, int n)
{
    parsetree_release(dest->m);
    dest->m = parsetree_copy(src->m);
    memcpy(dest->f, src->f + b, n * sizeof(double));
    memcpy(dest->a, src->a + b, n * sizeof(double));
    memcpy(dest->b, src->b + b, n * sizeof(double));
}





static void evaluate_model(struct model *ind,
                           const double * const *X, const double *t, double *z, int n,
                           const int * const *B, int nboot,
                           struct model **ens_elite)
{
    int j, b;
    for (j = 0; j < n; ++j) z[j] = parsetree_predict(ind->m, X[j]);
    for (b = 0; b < nboot; ++b) {
        ind->f[b] = fit(t, z, n, B[b], ind->a + b, ind->b + b);
        if (cmp_double(ind->f[b], ens_elite[b]->f[0]) < 0) replace_model(ens_elite[b], ind, b, 1);
    }
}





static void initial_generation(struct model **pop, int N,
                               const int * const *B, int nboot, struct model **ens_elite,
                               const double * const *trainX, const double *traint, int ntrain,
                               int min_init_height, int max_init_height,
                               const parsetree_nodeop *functions, const int *arity, int nfunc,
                               const parsetree_nodeop *terminals, int nterm,
                               double (*rnd)(void))
{
    int idx, i, j, k, half_pop_size = N / 2;
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
            evaluate_model(pop[idx], trainX, traint, z, ntrain, B, nboot, ens_elite);
        }
    }

    free(z);
}





static void breed_generation(struct model **pop, struct model **off, int N,
                             const int * const *B, int nboot, struct model **ens_elite,
                             const double * const *trainX, const double *traint, int ntrain,
                             double pm, double pc,
                             int max_mut_depth, int max_tree_height, int max_tree_size,
                             const parsetree_nodeop *functions, const int *arity, int nfunc,
                             const parsetree_nodeop *terminals, int nterm,
                             double (*rnd)(void))
{
    int i;
    double p, *z = MALLOC(ntrain, sizeof(double));
    bool eval;

    for (i = 0; i < N; ++i) {
        eval = false;
        p = rnd();
        if (p < (pm + pc)) {
            parsetree_release(off[i]->m);
            if (p < pm) {
                off[i]->m = parsetree_mutation(parsetree_copy(pop[i]->m), max_mut_depth, functions, arity, nfunc, terminals, nterm, rnd);
            } else {
                off[i]->m = parsetree_crossover(parsetree_copy(pop[i]->m), pop[(int)(rnd() * N)]->m, rnd);
            }
            eval = (max_tree_size < 0 || parsetree_size(off[i]->m) <= max_tree_size) && (max_tree_height < 0 || parsetree_height(off[i]->m) <= max_tree_height);
        }

        if (eval) {
            evaluate_model(off[i], trainX, traint, z, ntrain, B, nboot, ens_elite);
        } else {
            replace_model(off[i], pop[i], 0, nboot);
        }
    }

    free(z);
}





static int cmp_pop(const void *a_ptr, const void *b_ptr, void *arg)
{
    int i = *((const int *)arg);
    const struct model *a = *((const struct model **)a_ptr), *b = *((const struct model **)b_ptr);

    return cmp_double(a->f[i], b->f[i]);
}
static void selection(const struct model * const *pop, const struct model * const *off, struct model **gen, int N, int nboot,
                      double (*rnd)(void))
{
    int i, j, b;
    int *S = MALLOC(nboot, sizeof(int)), per_s = (int)round(N /(double)nboot);
    int sum_s = per_s * nboot;
    const struct model **sel = MALLOC(2 * N, sizeof(const struct model *));

    memcpy(sel, pop, N * sizeof(struct model *));
    memcpy(sel + N, off, N * sizeof(struct model *));

    for (b = 0; b < nboot; ++b) S[b] = per_s;
    while(sum_s < N) {
        do { i = (int)(rnd() * nboot); } while (S[i] < per_s);
        S[i]++;
        sum_s++;
    }
    while(sum_s > N) {
        do { i = (int)(rnd() * nboot); } while (S[i] == 1);
        S[i]--;
        sum_s--;
    }

    for (i = 0, b = 0; b < nboot; ++b) {
        qsort_r(sel, 2 * N, sizeof(struct model *), cmp_pop, &b);
        for (j = 0; j < S[b]; ++j, ++i) replace_model(gen[i], sel[j], 0, nboot);
    }

    free(sel);
    free(S);
}






static void tournament(const struct model * const *pop, const struct model * const *off, struct model **gen, int N, int nboot, int tournament_size, double (*rnd)(void))
{
    int i, j, k, c, b, pick;
    int *S = MALLOC(nboot, sizeof(int)), per_s = (int)round(N / (double)nboot);
    int sum_s = per_s * nboot;
    const struct model **sel = MALLOC(2 * N, sizeof(const struct model *));

    memcpy(sel, pop, N * sizeof(const struct model *));
    memcpy(sel + N, off, N * sizeof(const struct model *));

    for (b = 0; b < nboot; ++b) S[b] = per_s;
    while(sum_s < N) {
        S[(int)(rnd() * nboot)]++;
        sum_s++;
    }
    while(sum_s > N) {
        S[(int)(rnd() * nboot)]--;
        sum_s--;
    }

    for (i = 0, b = 0; b < nboot; ++b) {
        for (j = 0; j < S[b]; ++j, ++i) {
            pick = (int)(rnd() * 2 * N);
            for (k = 1; k < tournament_size; ++k) {
                c = (int)(rnd() * 2 * N);
                if (cmp_double(sel[c]->f[b], sel[pick]->f[b]) < 0) pick = c;
            }

            replace_model(gen[i], sel[pick], 0, nboot);
        }
    }

    free(sel);
}





static double elite_loss(const struct model * const *elite, int nboot,
                         const double * const *X, const double *t, int n)
{
    int b, j, c;
    double f, y, r, mse = 0;

    for (j = 0; j < n; ++j) {
        y = 0;
        for (c = 0, b = 0; b < nboot; ++b) {
            f = elite[b]->a[0] + elite[b]->b[0] * parsetree_predict(elite[b]->m, X[j]);
            if (sane(f)) {
                y += (f - y) / (++c);
            }
        }
        if (c == 0) y = NAN;

        r = (t[j] - y);
        mse += (r * r - mse) / (j + 1);
    }
    return mse;
}


static int elite_tree_size(const struct model * const *elite, int nboot)
{
    int b, n = 0;

    for (b = 0; b < nboot; ++b) n += parsetree_size(elite[b]->m);

    return n;
}


static int elite_unique(const struct model * const *elite, int nboot)
{
    int b, j, n = nboot;
    bool *dup = MALLOC(nboot, sizeof(bool));

    for (b = 0; b < nboot; ++b) dup[b] = false;

    for (b = 0; b < nboot; ++b) {
        if (dup[b]) continue;
        for (j = b + 1; j < nboot; ++j) {
            if (parsetree_equal(elite[b]->m, elite[j]->m)) dup[j] = true;
        }
    }

    for (b = 0; b < nboot; ++b) if (dup[b]) n--;

    free(dup);

    return n;
}


static void report(int g, const struct model * const *ens_elite, int nboot,
                   const double * const *trainX, const double *traint, int ntrain, double msr_train,
                   const double * const *testX, const double *testt, int ntest, double msr_test)
{
    double mse_train = elite_loss(ens_elite, nboot, trainX, traint, ntrain);
    double mse_test  = elite_loss(ens_elite, nboot, testX, testt, ntest);

    fprintf(stdout, "%3d %g %g %g %g %d %d\n",
            g,
            sqrt(mse_train), sqrt(mse_test),
            sqrt(mse_train / msr_train),
            sqrt(mse_test / msr_test),
            elite_tree_size(ens_elite, nboot),
            elite_unique(ens_elite, nboot));
    fflush(stdout);
}





struct ensemble *twosegp_evolve(
        const double * const *trainX, const double *traint, int ntrain, int nvar,
        const double * const *testX, const double *testt, int ntest,
        int N, int G, int nboot,
        double pc, double pm,
        int min_init_depth, int max_init_depth,
        int max_mut_depth,
        int max_tree_height, int max_tree_size,
        bool verbose,
        const parsetree_nodeop * functions, const int * arity, int nfunc,
        const parsetree_nodeop * terminals, int nterm,
        double (*rnd)(void))
{
    static const parsetree_nodeop def_functions[] = { exec_add, exec_sub, exec_mul, exec_div, exec_log, exec_sqrt };
    static const int def_arity[] = { 2, 2, 2, 2, 1, 1 };
    static const int def_nfunc = 6;

    static const parsetree_nodeop def_terminals[] = { exec_erc };
    static const int def_nterm = 1;

    int g = 0, b, i, j;
    struct ensemble *E = ensemble_create();
    struct model **pop, **off, **gen, **tmp, **ens_elite;
    double msr_train = var(traint, ntrain, false), msr_test = var(testt, ntest, false);
    int **B;
    bool use_tournament = false;

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
        nboot = 1;
        B = MALLOC(nboot, sizeof(int *));
        B[0] = CALLOC(ntrain, sizeof(int));
        for (j = 0; j < ntrain; ++j) B[0][j] = j;
        /* use_tournament = true; */
    } else if (nboot < 0) {
        nboot = abs(nboot);
        B = MALLOC(nboot, sizeof(int *));
        B[0] = CALLOC(nboot * ntrain, sizeof(int));
        for (b = 0; b < nboot; ++b) {
            B[b] = B[0] + b * ntrain;
            for (j = 0; j < ntrain; ++j) B[b][j] = j;
        }
    } else {
        B = MALLOC(nboot, sizeof(int *));
        B[0] = CALLOC(nboot * ntrain, sizeof(int));
        for (b = 0; b < nboot; ++b) {
            B[b] = B[0] + b * ntrain;
            for (j = 0; j < ntrain; ++j) B[b][j] = (int)(rnd() * ntrain);
        }
    }

    pop = create_population(N, nboot);
    off = create_population(N, nboot);
    gen = create_population(N, nboot);
    ens_elite = create_population(nboot, 1);

    initial_generation(pop, N, (const int * const *)B, nboot, ens_elite,
                       trainX, traint, ntrain,
                       min_init_depth, max_init_depth,
                       int_functions, int_arity, nfunc, int_terminals, nterm,
                       rnd);

    if (verbose) report(g, (const struct model * const *)ens_elite, nboot, trainX, traint, ntrain, msr_train, testX, testt, ntest, msr_test);

    for (g = 1; g <= G; ++g) {
        breed_generation(pop, off, N, (const int * const *)B, nboot, ens_elite,
                         trainX, traint, ntrain,
                         pm, pc, max_mut_depth, max_tree_height, max_tree_size,
                         int_functions, int_arity, nfunc, int_terminals, nterm,
                         rnd);

        if (use_tournament) {
            tournament((const struct model * const *)pop, (const struct model * const *)off, gen, N, nboot, 8, rnd);
        } else {
            selection((const struct model * const *)pop, (const struct model * const *)off, gen, N, nboot, rnd);
        }

        tmp = pop;
        pop = gen;
        gen = tmp;
        if (verbose) report(g, (const struct model * const *)ens_elite, nboot, trainX, traint, ntrain, msr_train, testX, testt, ntest, msr_test);
    }

    for (b = 0; b < nboot; ++b) ensemble_add_model(E, ens_elite[b]->m, ens_elite[b]->a[0], ens_elite[b]->b[0]);

    release_population(ens_elite, nboot);
    release_population(gen, N);
    release_population(off, N);
    release_population(pop, N);

    free(B[0]);
    free(B);
    free(var_terminals);

    return E;
}
