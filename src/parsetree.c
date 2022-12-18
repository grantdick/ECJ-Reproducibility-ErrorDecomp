#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#include "alloc.h"
#include "parsetree.h"
#include "util.h"

struct parsetree {
    parsetree_nodeop op;
    double data;

    int argc;
    struct parsetree **argv;
    struct parsetree *parent;

    int depth;
    int height;
    int size;
};

double exec_add(const struct parsetree *tree, const double *Xi)
{
    return parsetree_predict(tree->argv[0], Xi) + parsetree_predict(tree->argv[1], Xi);
}

double exec_sub(const struct parsetree *tree, const double *Xi)
{
    return parsetree_predict(tree->argv[0], Xi) - parsetree_predict(tree->argv[1], Xi);
}

double exec_mul(const struct parsetree *tree, const double *Xi)
{
    return parsetree_predict(tree->argv[0], Xi) * parsetree_predict(tree->argv[1], Xi);
}

double exec_div(const struct parsetree *tree, const double *Xi)
{
    double b = parsetree_predict(tree->argv[1], Xi);
    return copysign(1, b) * parsetree_predict(tree->argv[0], Xi) / (fabs(b) + 1.0e-6);
}

double exec_aqt(const struct parsetree *tree, const double *Xi)
{
    double b = parsetree_predict(tree->argv[1], Xi);
    return parsetree_predict(tree->argv[0], Xi) / sqrt(1.0 + b * b);
}

double exec_sin(const struct parsetree *tree, const double *Xi)
{
    return sin(parsetree_predict(tree->argv[0], Xi));
}
double exec_cos(const struct parsetree *tree, const double *Xi)
{
    return cos(parsetree_predict(tree->argv[0], Xi));
}

double exec_exp(const struct parsetree *tree, const double *Xi)
{
    return exp(parsetree_predict(tree->argv[0], Xi));
}

double exec_sqrt(const struct parsetree *tree, const double *Xi)
{
    return sqrt(fabs(parsetree_predict(tree->argv[0], Xi)));
}

double exec_log(const struct parsetree *tree, const double *Xi)
{
    return log(fabs(parsetree_predict(tree->argv[0], Xi)) + 1.0e-2);
}


double exec_erc(const struct parsetree *tree, const double *Xi __attribute__((unused)))
{
    return tree->data;
}

double exec_var(const struct parsetree *tree, const double *Xi)
{
    return Xi[(int)(tree->data)];
}










void parsetree_release(struct parsetree *tree)
{
    if (tree == NULL) return;

    while (tree->argc > 0) parsetree_release(tree->argv[--tree->argc]);
    free(tree->argv);
    free(tree);
}





struct parsetree *parsetree_copy(const struct parsetree *tree)
{
    struct parsetree *cpy;

    if (tree == NULL) return NULL;

    cpy = MALLOC(1, sizeof(struct parsetree));
    cpy->parent = NULL;

    cpy->depth      = tree->depth;
    cpy->height     = tree->height;
    cpy->size       = tree->size;

    cpy->op         = tree->op;
    cpy->data       = tree->data;

    cpy->argv       = MALLOC(tree->argc, sizeof(struct parsetree *));
    cpy->argc       = 0;

    while (cpy->argc < tree->argc) {
        cpy->argv[cpy->argc] = parsetree_copy(tree->argv[cpy->argc]);
        cpy->argv[cpy->argc]->parent = cpy;
        cpy->argc++;
    }

    return cpy;
}





static int collect_subtrees(struct parsetree *root, int required_depth, int N, struct parsetree **out)
{
    int i;
    if (root->depth == required_depth) {
        out[N++] = root;
    } else {
        for (i = 0; i < root->argc; ++i) N = collect_subtrees(root->argv[i], required_depth, N, out);
    }
    return N;
}
static struct parsetree *pick_subtree(struct parsetree *root, double (*rnd)(void))
{
    struct parsetree **res = MALLOC(root->size, sizeof(struct parsetree *));
    int required_depth = (int)(rnd() * (root->height + 1));
    int N = collect_subtrees(root, required_depth, 0, res);
    struct parsetree *pick = res[(int)(rnd() * N)];

    free(res);
    return pick;
}


static int remove_child(struct parsetree *parent, struct parsetree *child)
{
    int arg = 0;

    while (parent->argv[arg] != child) arg++;

    parsetree_release(parent->argv[arg]);
    parent->argv[arg] = NULL;

    return arg;
}



static void insert_child(struct parsetree *parent, struct parsetree *child, int arg)
{
    child->parent = parent;
    parent->argv[arg] = child;
}



static void append_child(struct parsetree *parent, struct parsetree *child)
{
    insert_child(parent, child, parent->argc++);
}



static void update_tree_stats(struct parsetree *tree)
{
    int i;

    while (tree) {
        tree->height = 0;
        tree->size = 1;
        for (i = 0; i < tree->argc; ++i) {
            tree->size += tree->argv[i]->size;
            if (tree->height < (tree->argv[i]->height + 1)) tree->height = tree->argv[i]->height + 1;
        }
        tree = tree->parent;
    }
}
static void update_tree_depths(struct parsetree *tree) {
    int i;

    tree->depth = (tree->parent == NULL) ? 0 : tree->parent->depth + 1;
    for (i = 0; i < tree->argc; ++i) {
        tree->argv[i]->parent = tree;
        update_tree_depths(tree->argv[i]);
    }
}





struct parsetree *parsetree_create(int cur_depth, int max_height, bool grow,
                                   const parsetree_nodeop * functions, const int * arity, int nfunc,
                                   const parsetree_nodeop * terminals, int nterm,
                                   double (*rnd)(void))
{
    int i, pick;
    struct parsetree *tree = MALLOC(1, sizeof(struct parsetree));

    tree->depth = cur_depth;
    tree->height = 0;
    tree->size = 1;
    tree->parent = NULL;
    tree->argv = NULL;
    tree->argc = 0;

    if (cur_depth == max_height) {
        pick = nfunc + (int)(rnd() * nterm);
    } else if (grow) {
        pick = (int)(rnd() * (nfunc + nterm));            
    } else {
        pick = (int)(rnd() * nfunc);
    }

    if (pick < nfunc) {
        tree->op = functions[pick];
        tree->data = NAN;
        tree->argv = MALLOC(arity[pick], sizeof(struct parsetree *));
        for (i = 0; i < arity[pick]; ++i) {
            append_child(tree, parsetree_create(cur_depth + 1, max_height, grow, functions, arity, nfunc, terminals, nterm, rnd));
            tree->size += tree->argv[i]->size;
            if (tree->height < (tree->argv[i]->height + 1)) tree->height = tree->argv[i]->height + 1;
        }
    } else {
        pick -= nfunc;
        tree->op = terminals[pick];
        tree->argc = 0;
        tree->data = (tree->op == exec_var) ? pick : 0.001 * (-5000 + (int)(rnd() * 10000));
    }

    return tree;
}





struct parsetree *parsetree_mutation(struct parsetree *tree, int max_mut_height,
                                     const parsetree_nodeop * functions, const int * arity, int nfunc,
                                     const parsetree_nodeop * terminals, int nterm,
                                     double (*rnd)(void))
{
    struct parsetree *mut, *rep, *p;

    mut = parsetree_create(0, max_mut_height, true, functions, arity, nfunc, terminals, nterm, rnd);
    rep = pick_subtree(tree, rnd);

    p = rep->parent;
    if (p == NULL) {
        parsetree_release(tree);
        return mut;
    }

    insert_child(p, mut, remove_child(p, rep));

    update_tree_stats(p);
    update_tree_depths(mut);

    return tree;
}





struct parsetree *parsetree_crossover(struct parsetree *m, const struct parsetree *f,
                                      double (*rnd)(void))
{
    struct parsetree *n1, *n2, *p;

    n1 = pick_subtree(m, rnd);
    n2 = parsetree_copy(pick_subtree((struct parsetree *)f, rnd));

    p = n1->parent;
    if (p == NULL) {
        parsetree_release(m);
        update_tree_depths(n2);
        return n2;
    }

    insert_child(p, n2, remove_child(p, n1));

    update_tree_stats(p);
    update_tree_depths(n2);

    return m;
}





double parsetree_predict(const struct parsetree *tree, const double *X)
{
    return tree->op(tree, X);
}





bool parsetree_equal(const struct parsetree *a, const struct parsetree *b)
{
    int i;
    if (a && b) {
        if (a->op != b->op) return false;
        if (a->size != b->size) return false;
        if (a->height != b->height) return false;
        if (a->argc > 0) {
            for (i = 0; i < a->argc; ++i) if (!parsetree_equal(a->argv[i], b->argv[i])) return false;
        } else {
            if (a->data != b->data) return false;
        }
        return true;
    } else {
        return false;
    }
}





int parsetree_depth(const struct parsetree *tree)
{
    return (tree == NULL) ? -1 : tree->depth;
}





int parsetree_height(const struct parsetree *tree)
{
    return (tree == NULL) ? -1 : tree->height;
}





int parsetree_size(const struct parsetree *tree)
{
    return tree->size;
}





void parsetree_print(FILE *out, const struct parsetree *tree)
{
    if (tree->op == exec_add) {
        fprintf(out, "(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, "+");
        parsetree_print(out, tree->argv[1]);
        fprintf(out, ")");
    } else if (tree->op == exec_sub) {
        fprintf(out, "(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, "-");
        parsetree_print(out, tree->argv[1]);
        fprintf(out, ")");
    } else if (tree->op == exec_mul) {
        fprintf(out, "(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, "*");
        parsetree_print(out, tree->argv[1]);
        fprintf(out, ")");
    } else if (tree->op == exec_div) {
        fprintf(out, "(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, "/");
        parsetree_print(out, tree->argv[1]);
        fprintf(out, ")");
    } else if (tree->op == exec_sin) {
        fprintf(out, "sin(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, ")");
    } else if (tree->op == exec_cos) {
        fprintf(out, "cos(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, ")");
    } else if (tree->op == exec_exp) {
        fprintf(out, "exp(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, ")");
    } else if (tree->op == exec_sqrt) {
        fprintf(out, "sqrt(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, ")");
    } else if (tree->op == exec_log) {
        fprintf(out, "log(");
        parsetree_print(out, tree->argv[0]);
        fprintf(out, ")");
    } else if (tree->op == exec_var) {
        fprintf(out, "x%.0f", tree->data + 1);
    } else {
        fprintf(out, "%f", tree->data);
    }
}
