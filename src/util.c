#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>

#include <math.h>

#include "alloc.h"

#define BIG_NUMBER  1.0e+15










int *sample_m_of_n(int m, int n, double (*rnd)(void))
{
    int i, j;

    int *sub;

    sub = MALLOC(m + 1, sizeof(int));

    for (j = 0, i = 0; i < n && j < m; ++i) {
        if (rnd() < ((double)(m - j) / (double)(n - i))) sub[j++] = i;
    }

    sub[m] = -1; /* to indicate the end of the sample */

    return sub;
}





int *sample_p_of_n(double p, int n, double (*rnd)(void))
{
    if (p < 0) {
        fprintf(stderr, "%s:%d - WARNING - supplied value of %f for p, trimming to zero and returning NULL\n", __FILE__, __LINE__, p);
        p = 0;
    } else if (p > 1) {
        fprintf(stderr, "%s:%d - WARNING - supplied value of %f for p, truncating to 1 and returning all\n", __FILE__, __LINE__, p);
        p = 1;
    }

    return sample_m_of_n((int)(p * n), n, rnd);
}





bool *sample_indices(int n, double train_frac, double (*rnd)(void))
{
    int *sample = sample_p_of_n(train_frac, n, rnd), i;
    bool *leave_in = MALLOC(n, sizeof(bool));

    for (i = 0; i < n; ++i) leave_in[i] = false;
    for (i = 0; sample[i] >= 0; ++i) leave_in[sample[i]] = true;

    free(sample);

    return leave_in;
}





bool sane(double d)
{
    if (!isfinite(d)) return false;
    if (fabs(d) >= BIG_NUMBER) return false;
    return true;
}





int cmp_double(double a, double b)
{
    if (sane(a) && sane(b)) {
        return (a < b) ? -1 : (a > b) ? 1 : 0;
    } else if (sane(a)) {
        return -1;
    } else if (sane(b)) {
        return  1;
    } else {
        return  0;
    }
}





static double quickselect(double *arr, unsigned long n, unsigned long k)
{
#define SWAP(a,b) { temp=(a); (a)=(b); (b)=temp; }

    unsigned long i, ir, j, l, mid;
    double a, temp;

    l = 0;
    ir = n-1;
    while (true) {
        if (ir <= l+1) {
            if (ir == l+1 && cmp_double(arr[ir], arr[l]) < 0) SWAP(arr[l], arr[ir]);
            return arr[k];
        } else {
            mid = (l+ir) >> 1;
            SWAP(arr[mid], arr[l+1]);
            if (cmp_double(arr[l], arr[ir]) > 0) SWAP(arr[l], arr[ir]);
            if (cmp_double(arr[l+1], arr[ir]) > 0) SWAP(arr[l+1], arr[ir]);
            if (cmp_double(arr[l], arr[l+1]) >0) SWAP(arr[l], arr[l+1]);

            i = l+1;
            j = ir;
            a = arr[l+1];
            for (;;) {
                do i++; while (cmp_double(arr[i], a) < 0);
                do j--; while (cmp_double(arr[j], a) > 0);
                if (j < i) break;
                SWAP(arr[i], arr[j]);
            }
            arr[l+1] = arr[j];
            arr[j] = a;
            if (j >= k) ir = j-1;
            if (j <= k) l = i;
        }
    }
}

double median(const double *X, int n)
{
    double *arr, med;

    if (n == 0) return NAN;
    if (n == 1) return X[0];

    arr = MALLOC(n, sizeof(double));
    memcpy(arr, X, n * sizeof(double));

    if (n % 2 == 0) {
        med = 0.5 * (quickselect(arr, n, n / 2 - 1) + quickselect(arr, n, n / 2));
    } else {
        med = quickselect(arr, n, n / 2);
    }

    free(arr);

    return med;
}





void cummean(double *X, int n)
{
    int i;
    double m = X[0];

    for (i = 1; i < n; ++i) {
        m += (X[i] - m) / (i + 1);
        X[i] = m;
    }
}





double mean(const double *X, int n)
{
    int i;
    double m;

    m = 0;

    for (i = 0; i < n; ++i) m += (X[i] - m) / (i + 1);

    return m;
}





double cov(const double *X, const double *Y, int n, bool sample_var)
{
    int i;
    double delta, mx, my, C;

    mx = my = C = 0;
    for (i = 0; i < n; ++i) {
        delta = X[i] - mx;

        mx += delta / (i + 1);
        my += (Y[i] - my) / (i + 1);
        C += delta * (Y[i] - my);
    }

    return sample_var ? (C / (n - 1)) : (C / n);
}





double var(const double *X, int n, bool sample_var)
{
    return cov(X, X, n, sample_var);
}





double fit(const double * const t, const double * const y, int n, const int * const B, double *a_ptr, double *b_ptr)
{
    int i, j;
    double dt, dy, mt = 0, my = 0, st = 0, sy = 0, C = 0;
    double a, b, r, mse = 0;

    for (j = 0; j < n; ++j) {
        i = B ? B[j] : j;

        dt = t[i] - mt;
        dy = y[i] - my;

        mt += dt / (j + 1);
        my += dy / (j + 1);

        C  += dt * (y[i] - my);
        st += dt * (t[i] - mt);
        sy += dy * (y[i] - my);
    }

    b = (fabs(C) < 1.0e-7 || fabs(sy) < 1.0e-7) ? 0 : C / sy;
    a = mt - b * my;

    for (j = 0; j < n; ++j) {
        i = B ? B[j] : j;

        r = t[i] - (a + b * y[i]);
        mse += (r * r - mse) / (j + 1);
    }

    if (a_ptr) *a_ptr = a;
    if (b_ptr) *b_ptr = b;

    return mse;
}





double mse(const double * const t, const double * const y, int n)
{
    int i;
    double m = 0, r;

    for (i = 0; i < n; ++i) {
        r = t[i] - y[i];
        m += (r * r - m) / (i + 1.0);
    }

    return m;
}





double minof(int n, ...)
{
    int i;

    va_list ap;

    double d, m;

    if (n == 0) return NAN;

    va_start(ap, n);

    m = INFINITY;

    for (i = 0; i < n; ++i) {
        d = va_arg(ap, double);
        if (!sane(d)) return NAN;
        if (d < m) m = d;
    }

    va_end(ap);

    return m;
}





double maxof(int n, ...)
{
    int i;

    va_list ap;

    double d, m;

    if (n == 0) return NAN;

    va_start(ap, n);

    m = -INFINITY;

    for (i = 0; i < n; ++i) {
        d = va_arg(ap, double);
        if (!sane(d)) return NAN;
        if (d > m) m = d;
    }

    va_end(ap);

    return m;
}
