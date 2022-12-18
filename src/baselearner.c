#include <stdlib.h>

#include <math.h>

#include "alloc.h"
#include "baselearner.h"
#include "parsetree.h"
#include "util.h"

struct baselearner
{
    double intercept;
    double slope;

    struct parsetree *model;
};










struct baselearner *baselearner_create()
{
    struct baselearner *l = MALLOC(1, sizeof(struct baselearner));

    l->intercept = 0;
    l->slope = 1;
    l->model = NULL;

    return l;
}





void baselearner_release(struct baselearner *l)
{
    if (l == NULL) return;

    parsetree_release(l->model);

    free(l);
}





void baselearner_set_model(struct baselearner *l, struct parsetree *model, double intercept, double slope)
{
    parsetree_release(l->model);
    l->model = model;

    if (sane(slope) && log10(fabs(slope)) < 7) {
        l->intercept = intercept;
        l->slope = slope;
    } else {
        l->intercept = l->slope = NAN;
    }
}





const struct parsetree *baselearner_get_model(const struct baselearner *l)
{
    return l->model;
}





struct baselearner *baselearner_copy(const struct baselearner *l)
{
    struct baselearner *cpy = baselearner_create();

    cpy->intercept = l->intercept;
    cpy->slope = l->slope;
    cpy->model = parsetree_copy(l->model);

    return cpy;
}





double baselearner_transform(const struct baselearner *l, double z)
{
    return l->intercept + l->slope * z;
}





double baselearner_predict(const struct baselearner *l, const double *Xi)
{
    if (l->model == NULL) return NAN;
    return l->intercept + l->slope * parsetree_predict(l->model, Xi);
}





void baselearner_print(FILE *out, const struct baselearner *l)
{
    fprintf(out, "%f + %f * ", l->intercept, l->slope);
    parsetree_print(out, l->model);
}
