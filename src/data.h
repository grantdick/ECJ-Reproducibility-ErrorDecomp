#ifndef DATA_H
#define DATA_H

#ifdef  __cplusplus
extern "C" {
#endif

    void generate_data(const char *problem,
                       double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                       double ***testX_ptr, double **testt_ptr, int *ntest_ptr,
                       double ***intervals_ptr,
                       int *p_ptr, double (*rnd)(void), double (*rnd_gauss)(double, double));

    void read_data(const char *src,
                   double ***X_ptr, double **t_ptr, double ***intervals_ptr,
                   int *n_ptr, int *p_ptr);

    void release_data(double **X, double *t, double **intervals);

    void configure_data(const char *split_file, int split_idx,
                        double **X, double *t, int n, int p,
                        double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                        double ***testX_ptr, double **testt_ptr, int *ntest_ptr);

    void bootstrap_data(double **X, double *t, int n, int p,
                        double ***bootX_ptr, double **boott_ptr,
                        double (*rnd)(void));

    void split_data(double **X, double *t, bool *leave_in, int n, int p,
                    double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                    double ***testX_ptr, double **testt_ptr, int *ntest_ptr);

    void train_test_split(double **X, double *t, int *train_idx, int n, int p,
                    double ***trainX_ptr, double **traint_ptr, int *ntrain_ptr,
                    double ***testX_ptr, double **testt_ptr, int *ntest_ptr);

    void scale_data(double **trainX, int ntrain, int p,
                    double **testX, int ntest,
                    bool bias_correct,
                    double **mu_ptr, double **s_ptr);

#ifdef  __cplusplus
}
#endif

#endif
