#ifndef BASELEARNER_H
#define BASELEARNER_H

#ifdef  __cplusplus
extern "C" {
#endif

    #include "parsetree.h"

    struct baselearner; /* opaque structure - definition not needed by external users */

    struct baselearner *baselearner_create();

    void baselearner_release(struct baselearner *l);

    void baselearner_set_model(struct baselearner *l, struct parsetree *model, double intercept, double slope);

    const struct parsetree *baselearner_get_model(const struct baselearner *l);

    struct baselearner *baselearner_copy(const struct baselearner *l);

    double baselearner_transform(const struct baselearner *l, double z);
    
    double baselearner_predict(const struct baselearner *l, const double *Xi);

    void baselearner_print(FILE *out, const struct baselearner *l);

#ifdef  __cplusplus
}
#endif

#endif
