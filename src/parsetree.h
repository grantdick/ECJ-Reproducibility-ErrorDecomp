#ifndef PARSETREE_H
#define PARSETREE_H

#ifdef  __cplusplus
extern "C" {
#endif

    #include <stdio.h>

    #include <stdbool.h>

    struct parsetree; /* opaque structure - definition not needed by external users */

    typedef double (*parsetree_nodeop)(const struct parsetree *, const double *);

    double exec_add(const struct parsetree *tree, const double *Xi);
    double exec_sub(const struct parsetree *tree, const double *Xi);
    double exec_mul(const struct parsetree *tree, const double *Xi);
    double exec_div(const struct parsetree *tree, const double *Xi);
    double exec_aqt(const struct parsetree *tree, const double *Xi);
    double exec_sin(const struct parsetree *tree, const double *Xi);
    double exec_cos(const struct parsetree *tree, const double *Xi);
    double exec_exp(const struct parsetree *tree, const double *Xi);
    double exec_sqrt(const struct parsetree *tree, const double *Xi);
    double exec_log(const struct parsetree *tree, const double *Xi);
    double exec_erc(const struct parsetree *tree, const double *Xi);
    double exec_var(const struct parsetree *tree, const double *Xi);

    void parsetree_release(struct parsetree *tree);

    struct parsetree *parsetree_copy(const struct parsetree * tree);

    struct parsetree *parsetree_create(int cur_depth, int max_height, bool grow,
                                       const parsetree_nodeop * functions, const int * arity, int nfunc,
                                       const parsetree_nodeop * terminals, int nterm,
                                       double (*rnd)(void));

    struct parsetree *parsetree_mutation(struct parsetree *tree, int max_mut_height,
                                        const parsetree_nodeop * functions, const int * arity, int nfunc,
                                        const parsetree_nodeop * terminals, int nterm,
                                        double (*rnd)(void));

    struct parsetree *parsetree_crossover(struct parsetree *m, const struct parsetree *f,
                                        double (*rnd)(void));


    double parsetree_predict(const struct parsetree *tree, const double *X);

    bool parsetree_equal(const struct parsetree *a, const struct parsetree *b);

    int parsetree_depth(const struct parsetree *tree);
    int parsetree_height(const struct parsetree *tree);
    int parsetree_size(const struct parsetree *tree);

    void parsetree_print(FILE *out, const struct parsetree *tree);

#ifdef  __cplusplus
}
#endif

#endif
