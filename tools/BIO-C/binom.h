/* beta.c */
double fmin2(double x, double y);
double fmax2(double x, double y);
double stirlerr(double n);
double log1p(double x);
double expm1(double x);
double chebyshev_eval(double x, const double *a, const int n);
double gammafn(double x);
double lgammafn(double x);
double lgammacor(double x);
double lbeta(double a, double b);
double pbeta_raw(double x, double pin, double qin, int lower_tail);
double pbeta(double x, double pin, double qin, int lower_tail, int log_p);
double pbinom(double x, double n, double p, int lower_tail, int log_p);
double cumbinom(int k, int n, double p);
int imax2(int x, int y);
double wilcoxtest(double* x, int nx, double* y, int ny, int alter);
double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p);
