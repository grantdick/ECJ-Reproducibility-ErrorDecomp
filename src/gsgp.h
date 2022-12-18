#ifndef GSGP_H
#define GSGP_H

#ifdef  __cplusplus
extern "C" {
#endif

    #include "parsetree.h"

    double *gsgp_evolve(const double **trainX, const double *traint, int ntrain, int nvar,
                        const double **testX, const double *testt, int ntest,
                        int N, int G, int tourn_size,
                        double pc, double pm,
                        int min_init_depth, int max_init_depth,
                        int max_mut_depth, double mutation_scale,
                        const parsetree_nodeop * const functions, const int * const arity, int nfunc,
                        const parsetree_nodeop * const terminals, int nterm,
                        void (*report)(int, const double *, double *, int, double, const double *, double *, int, double, void *),
                        double (*rnd)(void), void *report_args);

#ifdef  __cplusplus
}
#endif

#endif
