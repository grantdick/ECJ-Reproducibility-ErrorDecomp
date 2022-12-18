#ifndef UTIL_H
#define UTIL_H

#ifdef  __cplusplus
extern "C" {
#endif

    #include <stdbool.h>

    int *sample_m_of_n(int m, int n, double (*rnd)(void));

    int *sample_p_of_n(double p, int n, double (*rnd)(void));

    bool *sample_indices(int n, double train_frac, double (*rnd)(void));

    bool sane(double d);

    int cmp_double(double a, double b);

    double median(const double *X, int n);

    void cummean(double *X, int n);

    double mean(const double *X, int n);

    double cov(const double *X, const double *Y, int n, bool sample_var);

    double var(const double *X, int n, bool sample_var);

    double fit(const double * const t, const double * const y, int n, const int * const B, double *a_ptr, double *b_ptr);

    double mse(const double * const t, const double * const y, int n);

    double minof(int n, ...);

    double maxof(int n, ...);

#ifdef  __cplusplus
}
#endif

#endif
