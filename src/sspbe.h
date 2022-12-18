#ifndef SSPBE_H
#define SSPBE_H

#ifdef  __cplusplus
extern "C" {
#endif

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
                                double (*rnd)(void));

#ifdef  __cplusplus
}
#endif

#endif
