#ifndef SSPBE_H
#define SSPBE_H

#ifdef  __cplusplus
extern "C" {
#endif

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
            double (*rnd)(void));

#ifdef  __cplusplus
}
#endif

#endif
