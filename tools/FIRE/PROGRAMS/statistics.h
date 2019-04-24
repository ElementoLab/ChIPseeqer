#ifndef STATISTICS_H
#define STATISTICS_H

/* statistics.c */

#define M_SQRT_32     5.656854249492380195206754896838
#define M_1_SQRT_2PI  0.398942280401432677939946059934
#define DBL_EPSILON 2.2204460492503131e-16




void  default_set_seed(int seed);
float default_rand(); 
void  ran1_set_seed(long seed);
float ran1_rand(); 
   
void  quick_set_seed(unsigned long t);
float quick_rand(); 
 
void get_all_permutations(int n, int*** perm, int* nperms);
int eq_perm(int* a, int* b, int n); 

void myqsort (void *a, unsigned long n, unsigned long es, int (*cmp)());
void bubbleSort(int *a, int n);
int* bubbleSortIndex(float *ia, int n);
int* bubbleSortIndexInt(int *ia, int n);
 
void marsa_set_seed(unsigned int i1, unsigned int i2);
void marsa_get_seed(unsigned int *i1, unsigned int *i2);
double marsa_rand(void);

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
double ad_test(float* x, int n);
void writeFloatSet(char* file, float* f, int n); 
//int CmpInt(int* a, int* b);
int CmpInt(const void* a, const void* b);

float median_int (int* d, int n); 
float median (float* d, int n); 
float sum(float* sample, int n);
float pearson_int(int* data1, int* data2, int m); 
float pearson(float* data1, float* data2, int m); 

void varnorm( float* sample, int n, float** sample_r ); 


int CmpFloatRegular(const void* _a, const void* _b);

void count_for_hypergeom(int* v1, int* v2, int n, int* ov, int* s1, int* s2);
void transpose(int** m1, int n, int m, int*** m2);
void transpose_f(float** m1, int n, int m, float*** m2);

float welch_t(float *x1, int n1, float *x2, int n2);
float dof_t(float *x1, int n1, float *x2, int n2);
double sqr(double x);
float average(float *sample, int n);
float stddev_short ( short* sample, int n ); 
float average_short ( short* sample, int n ); 

float weighted_average(float* sample, int n, float* weights);
float stddev(float *sample, int n);
double stddev_dbl ( double* sample, int n );
double average_dbl(double* sample, int n); 
float average_int(int* sample, int n); 



double MannWhitney(double *sample1, int n1, double *sample2, int n2, int t);
int CmpFuncIndividuals(const void *_a, const void *_b);
float euclidean(float* data1, float* data2, int m);

double normal(const double z);
double lcumhyper(int i, int s1, int s2, int N);
double cumhyper(int i, int s1, int s2, int N);
double cumhyper_u(int i, int s1, int s2, int N);

void nrerror(char str[]);
double hypergeom(int i, int s1, int s2, int N);
double factln(int n);
double gammln(double xx);
double factrl(int n);
double cumbino(int k, int N, double p);
double bico(int n, int k);
double lbico(int n, int k);
double lcumbino(int k, int N, double p);
double lbino(int k, int N, double p);
void   getBiPermutations(int n1, int n2, int ***perms1, int ***perms2, int *nc);
void   _recBiPermutations(int ic, int n1, int jc, int n2, int **perm1, int **perm2, int **ref1, int **ref2, int ***perms1, int ***perms2, int *nc);
float *getBootstrappedVector(float *sample, int n);
int *  shuffle(int *v, int n);
int*   shuffleInt(int* v, int n); 

float* shuffleFloat(float* v, int n); 

int    *getIndices(int n);
int    *randPermOfIndices(int n);
void   fitSimpleLinearModel(float *x, float *y, int n, float *b, float *a);
void   fitDirectLinearModel(float* x, float* y, int n, float* b); 

int    maxArrayIndex (int* a, int n); 
int*   bubbleSortIndex(float *a, int n); 
int    CmpFloat(const void* _a, const void* _b);
int    CmpDblRegular(const void* _a, const void* _b);
float* intersect_binary_vector_f(float* v1, float* v2, int n); 
float* union_binary_vector_f    (float* v1, float* v2, int n); 
double median_dbl (double* d, int n);
double quantile_dbl (double* d, int n, float q); 
 

#endif
