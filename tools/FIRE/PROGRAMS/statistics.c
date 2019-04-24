


#define min(a, b)	((a) < (b) ? a : b)


typedef struct _Individual {
  
  float value;
  short dataset; 
  double   rank;

} Individual;


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>


#include "dataio.h"
#include "statistics.h"
#include "information.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

static unsigned long quick_seed = 12345;


static long          ran1_seed = -12345;




void default_set_seed(int seed) {
  //marsa_set_seed(1234, 5678); 
  //quick_set_seed(seed);
  ran1_set_seed(seed);
}


float default_rand() {
  //return marsa_rand();
  //return quick_rand();
  return ran1_rand();
}


void ran1_set_seed(long seed)
{
  if (seed == 0)
    die("Seed must be non-zero\n");

  if (seed < 0)
    ran1_seed = seed;
  else   
    ran1_seed = -seed;
}

float ran1_rand() 
{
  
  int    j;
  long   k;
  static long iy = 0;
  static long iv[32];
  float  temp;

    
  if ((ran1_seed <= 0) || (!iy)) {

    if (-ran1_seed < 1)
      ran1_seed = 1;
    else 
      ran1_seed = -ran1_seed;

    for (j=39; j>=0; j--) {

      k = ran1_seed/127773;
      ran1_seed = 16807 * (ran1_seed - k * 127773) - 2836 * k;
      if (ran1_seed < 0)
	ran1_seed += 2147483647;
      if (j < 32)
	iv[j] = ran1_seed;
    }
    iy = iv[0];
  }

  k = ran1_seed / 127773;
  ran1_seed = 16807 * (ran1_seed - k * 127773) - 2836 * k;

  if (ran1_seed < 0) 
    ran1_seed += 2147483647;

  j     = iy / (1 + ( 2147483647 - 1) / 32);
  iy    = iv[j];
  iv[j] = ran1_seed;

  if ((temp = (1.0/2147483647)*iy) > (1.0 - 1.2e-7)) 
    return (1.0 - 1.2e-7);
  else 
    return temp;
}






void quick_set_seed(unsigned long t) 
{
  quick_seed = t;
}

float quick_rand() 
{
  
  float rand;
  unsigned long temp;
  static unsigned long m1 = 0x3f800000;
  static unsigned long m2 = 0x007fffff;

  quick_seed = 1664525L*quick_seed + 1013904223L;
  temp = m1 | (m2 & quick_seed);
  rand = (float)((*(float *) & temp) - 1.0);
  if ((rand < 0.0) || (rand > 1.0)) {
    die("Problem with quick_rand(): random number < 0.0 or > 1.0.\n");
  }
  //  printf("%f\n", rand);

  return rand;
}




static unsigned int I1=1234, I2=5678;

void marsa_set_seed(unsigned int i1, unsigned int i2)
{
  I1 = i1; I2 = i2;
}

void marsa_get_seed(unsigned int *i1, unsigned int *i2)
{
  *i1 = I1; *i2 = I2;
}


double marsa_rand(void)
{
  double d ;

  I1= 36969*(I1 & 0177777) + (I1>>16);
  I2= 18000*(I2 & 0177777) + (I2>>16);
  
  d = ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10;
    
  //printf("%f\n", d);

  return d; 
}



void get_all_permutations(int n, int*** perm, int* np)
{
  int i, m, idx, j, perm_idx=0;

  int* a;
  int* b;
  int* c;
  int* t;

  int  nperms;

  a = (int*)malloc(n * sizeof(int));
  b = (int*)malloc(n * sizeof(int));
  c = (int*)malloc(n * sizeof(int));
  t = (int*)malloc(n * sizeof(int));

  nperms  = (int) (exp ( factln( n ) ));
  *np = nperms - 1;
  (*perm) = (int**)malloc(nperms * sizeof(int*));

  for (i=0; i<n; i++) {
    a[i] = i;
    b[i] = n-1-i;
  }
  
  while (! eq_perm( a, b, n )) {
    
    m = -1;
    for (m=n-2; m>=0; m--) {
      if (a[m+1] > a[m]) {
	break;
      }
    }
    

    if (m < 0) {
      die("m<0, should not happen\n");
    }

    // find the smallest value 
    // larger than a[m] to the right of a[m]
    idx = -1;
    for (i=m+1; i<n; i++) {
      if (a[i] > a[m]) {
	if (idx == -1) {
	  idx = i;
	} else {	
	  if (a[i] < a[idx]) {
	    idx = i;
	  }	
	}
      }  
    }
    

    for (i=0; i<m; i++) {
      c[i] = a[i];
    }
    
    c[m] = a[idx];
    
    j = 0;
    for (i=m; i<n; i++) {
      if (i == idx) {
	continue;
      } else {
	t[j] = a[i];
	j++;
      }
    }
    
    bubbleSort(t, j);
    
    for (i=m+1,j=0; i<n; i++, j++) {
      c[i] = t[j];
    }
    
    (*perm)[ perm_idx ] = (int*)malloc(n * sizeof(int));
    memcpy((*perm)[ perm_idx ], c, n* sizeof(int));
    perm_idx ++;
   
    memcpy(a, c, n*sizeof(int));
    
  }

  if (perm_idx != nperms-1) {
    printf("Problem ? perm_idx = %d != nbperm-s = %d\n", perm_idx, nperms-1);
  }

}


int eq_perm(int* a, int* b, int n) 
{
  int i;
  for (i=0; i<n; i++) {
    if (a[i] != b[i])
      return 0;
  }
  return 1;
}








void count_for_hypergeom(int* v1, int* v2, int n, int* ov, int* s1, int* s2) {
  int i;

  for (i=0; i<n; i++) {
    if (v1[i] == 1)
      (*s1) ++;
    if (v2[i] == 1)
      (*s2) ++;
    if ((v1[i] == 1) && (v2[i] == 1)) 
      (*ov) ++;
  }

}


void transpose_f(float** m1, int n, int m, float*** m2)
{
  int i, j;

  (*m2) = (float**)malloc(m * sizeof(float*));
  for (i=0; i<m; i++) {
    (*m2)[i] = (float*)malloc(n * sizeof(float));   
  }
  
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
      (*m2)[j][i] = m1[i][j];
    }
  }
  

}


void transpose(int** m1, int n, int m, int*** m2)
{
  int i, j;

  (*m2) = (int**)malloc(m * sizeof(int*));
  for (i=0; i<m; i++) {
    (*m2)[i] = (int*)malloc(n * sizeof(int));   
  }
  
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
      (*m2)[j][i] = m1[i][j];
    }
  }
  

}


float* intersect_binary_vector_f(float* v1, float* v2, int n) 
{
  int i;
  float* v3 = (float*)malloc(n * sizeof(float));
  for (i=0; i<n; i++) {
    v3[i] = ((v1[i]>0.0)&&(v2[i]>0.0)?1.0:0.0);
  }
  return v3;
}



float* union_binary_vector_f(float* v1, float* v2, int n) 
{
  int i;
  float* v3 = (float*)malloc(n * sizeof(float));
  for (i=0; i<n; i++) {
    v3[i] = ((v1[i]>0.0)||(v2[i]>0.0)?1.0:0.0);
  }
  return v3;
}






float welch_t(float* x1, int n1, float* x2, int n2) {
  
  float s1, s2, m1, m2, z;

  m1 = average(x1, n1); //printf("m1=%3.2f\n", m1);
  m2 = average(x2, n2); //printf("m2=%3.2f\n", m2);
    
  s1 = stddev(x1, n1); //printf("s1=%3.2f\n", s1);
  s2 = stddev(x2, n2); //printf("s2=%3.2f\n", s2);

  z = ( m1 - m2 ) / sqrt ( s1 * s1 / n1 + s2 * s2 / n2 );

  return z;

} 


float dof_t(float *x1, int n1, float *x2, int n2) {

  float f;
  float s1;
  float s2;

  s1 = stddev(x1, n1); //printf("s1=%3.2f\n", s1);
  s2 = stddev(x2, n2);

  f = sqr( (s1 * s1 / n1) + (s2 * s2 / n2) ) /
    ( sqr(s1 * s1 / n1) / (n1 - 1) +
      sqr(s2 * s2 / n2) / (n2 - 1) );

  return f;
}


double sqr(double x) 
{
  return x * x;
}



//
//  little helper function
//
int* copy_int_vector(int* d, int n) {
  int* d_copy;
  d_copy = (int*)malloc(n * sizeof(int));
  memcpy(d_copy, d, n * sizeof(int));
  return d_copy;
}


float* copy_float_vector(float* d, int n) {
  float* d_copy;
  d_copy = (float*)malloc(n * sizeof(float));
  memcpy(d_copy, d, n * sizeof(float));
  return d_copy;
}



double* copy_double_vector(double* d, int n) {
  double* d_copy;
  d_copy = (double*)malloc(n * sizeof(double));
  memcpy(d_copy, d, n * sizeof(double));
  return d_copy;
}



//
// 
//
double quantile_dbl (double* d, int n, float q) 
{

  int    nm;
  double* d_copy;
  
  if (n == 1) {
    return d[0];
  } else {
    d_copy = copy_double_vector(d, n);
    qsort((void*)d_copy, n, sizeof(double), CmpDblRegular);
    nm    = (int)(0.5+ n * q);
    return d_copy[nm];
  }
}


//
// 
//
double median_dbl (double* d, int n) 
{

  int    nm;
  double* d_copy;

  if (n == 1) {
    return d[0];
  } else {
    d_copy = copy_double_vector(d, n);
    
    qsort((void*)d_copy, n, sizeof(double), CmpDblRegular);
    
    nm    = n / 2 - 1;  // floor - 1
    
    printf("1:%4.3f, 2:%4.3f 3:%4.3f\n", d_copy[nm], d_copy[nm+1], d_copy[nm+2]); 
    
    if (n % 2 == 1)  {
      return d_copy[nm + 1];
    } else {
      return ( d_copy[nm] + d_copy[nm + 1] ) / 2.0;
    }

  }

}



//
// 
//
float median (float* d, int n) 
{

  int    nm;
  float* d_copy;

  if (n == 1) {
    return d[0];
  } else {
    d_copy = copy_float_vector(d, n);
    
    qsort((void*)d_copy, n, sizeof(float), CmpFloat);
    
    nm    = n / 2 - 1;  // floor - 1
    
    if (n % 2 == 1)  {
      return d_copy[nm + 1];
    } else {
      return ( d_copy[nm] + d_copy[nm + 1] ) / 2.0;
    }

  }

}



//
// 
//
float median_int (int* d, int n) 
{

  int    nm;
  int* d_copy;

  float me ;

  if (n == 1) {
    return d[0];
  } else {
    d_copy = copy_int_vector(d, n);
    
    qsort((void*)d_copy, n, sizeof(int), CmpInt);
    
    nm    = n / 2 - 1;  // floor - 1
        
    if (n % 2 == 1)  {
      me    = (float)(d_copy[nm + 1]);
    } else {
      me    = ( d_copy[nm] + d_copy[nm + 1] ) / 2.0;
    }
    
    free(d_copy);
    return me;
  }

}

int CmpInt(const void* _a, const void* _b)
{

  const int a = * ((const int*) _a);
  const int b = * ((const int*) _b);
  

  if (a<b)
    return -1;
  else if (a == b)
    return 0;
  else
    return 1;
}


double average_dbl(double* sample, int n) 
{

  //
  // calculate the sums
  //

  double  sum      = 0.0;
  int    i;
  int    n_actual = 0;

  for (i=0; i<n; i++) {
    
    if (!isnan(sample[i])) {
      sum   += (double)(sample[i]);
      n_actual ++;
    }

  }
  
  return sum / n_actual;
}


float average_short(short* sample, int n) 
{
  double  sum      = 0.0;
  int    i;
  int    n_actual = 0;
  for (i=0; i<n; i++) {
    sum   += (double)(sample[i]);
    n_actual ++;
  }  
  return (float)(sum / n_actual);
}


float average_int(int* sample, int n) 
{

  double  sum      = 0.0;
  int    i;
  int    n_actual = 0;

  for (i=0; i<n; i++) {
    
    //if (!isnan(sample[i])) {
    sum   += (double)(sample[i]);
    n_actual ++;
    //}

  }
  
  return (float)(sum / n_actual);
}

float average(float* sample, int n) 
{

  //
  // calculate the sums
  //

  double  sum      = 0.0;
  int    i;
  int    n_actual = 0;

  for (i=0; i<n; i++) {
    
    if (!isnan(sample[i])) {
      sum   += (double)(sample[i]);
      n_actual ++;
    }

  }
  
  return (float)(sum / n_actual);
}



float sum(float* sample, int n)
{
  float  sum      = 0.0;
  int    i;

  for (i=0; i<n; i++) {
    
    if (!isnan(sample[i])) {
      sum   += sample[i];
    }

  }
  
  return sum;


}


float weighted_average(float* sample, int n, float* weights) 
{
  
  // calculate the sums
  float sum   = 0.0;
  float sumw  = 0.0;
  
  int    i;

  // normalize the weights
  for (i=0; i<n; i++) {
    sumw  += weights[i];
  }

  for (i=0; i<n; i++) {
    sum   += weights[i] * sample[i] / sumw;
  }
  
  return sum;
}




double stddev_dbl ( double* sample, int n ) 
{
  
  
  // calculate the sums
  double sum   = 0.0;
  double sum_2 = 0.0;
  double std;
  int    i;

  for (i=0; i<n; i++) {
    sum   += sample[i];
    sum_2 += sample[i] * sample[i];
  }
  
  // calculate the standard deviation
  std = sqrt( (sum_2 - sum * sum / n ) / ( n - 1 )); 
  
  return std;
  
}


float stddev ( float* sample, int n ) 
{
  
  
  // calculate the sums
  double sum   = 0.0;
  double sum_2 = 0.0;
  double std;
  int    i;

  for (i=0; i<n; i++) {
    sum   += (double)( sample[i] );
    sum_2 += (double)(sample[i]) * (double)(sample[i] );
  }
  
  // calculate the standard deviation
  std = sqrt( (sum_2 - sum * sum / n ) / ( n - 1 )); 
  
  return (float) std;
  
}



float stddev_short ( short* sample, int n ) 
{
  
  
  // calculate the sums
  double sum   = 0.0;
  double sum_2 = 0.0;
  double std;
  int    i;

  for (i=0; i<n; i++) {
    sum   += (double)( sample[i] );
    sum_2 += (double)(sample[i]) * (double)(sample[i] );
  }
  
  // calculate the standard deviation
  std = sqrt( (sum_2 - sum * sum / n ) / ( n - 1 )); 
  
  return (float) std;
  
}


void varnorm( float* sample, int n, float** sample_r ) 
{

  float  a;
  float  s;
  int    i;
    
  a = average(sample, n);
  s = stddev (sample, n);
  
  *sample_r = (float*)malloc(n * sizeof(float));
  
  for (i=0; i<n; i++) {
    (*sample_r)[i] = ( sample[i] - a ) / s;
  }	

  
}




//
//   get a MW p-value for two populations
//
double MannWhitney(double* sample1, int n1, double* sample2, int n2, int t) {
  
    
  int    n;

  int    i, j, nbties;

  double    W1, W2, sum;

  double U1, U2;
  
  double mu, si, zi;
  
  Individual* indivs; 
  

  
  n = n1 + n2;
    
  //
  //  calculate the U statistics
  //
  
  //
  // put all the individuals in an array 
  //
  
  indivs = (Individual*)malloc(n * sizeof(Individual));
  for (i=0; i<n1; i++) {
    indivs[i].dataset = 1;
    indivs[i].value   = sample1[i];
    indivs[i].rank    = 0.0;
  }
  for (i=0, j=n1; i<n2; i++, j++) {
    indivs[j].dataset = 2;
    indivs[j].value   = sample2[i];
    indivs[j].rank    = 0.0;
  }
  
  for (i=0; i<n; i++) {
    //printf("indivs[%d], value=%3.2f, dataset=%d\n", i, indivs[i].value , indivs[i].dataset);
  }
  //getchar();

  qsort((void*)indivs, n, sizeof(Individual), CmpFuncIndividuals);

  W1 = 0;
  W2 = 0;
  for (i=0; i<n; i++) {
    //printf("indivs[%d], value=%3.2f, dataset=%d\n", i, indivs[i].value , indivs[i].dataset);

    if (indivs[i].rank > 0.0) {
      continue;
    }
    
    //
    // determine all the values ahead of the current one, that have the same value and no rank yet defined
    //


    j   = i+1;
    sum = i+1;
    while ( (indivs[j].rank == 0.0) && (indivs[j].value == indivs[i].value) ) {
      sum += j+1; // add the rank for this indiv
      j ++;
    }

    //
    // if there was a tie
    //
    if (j > (i+1)) {

      nbties = j - i;
      
      sum /= nbties;
      
      for (j=i; j<i+nbties; j++) {
	indivs[j].rank = sum;
      }
    } else {
      indivs[i].rank = (double)(i+1);
    }
  }
  
  //for (i=0; i<n; i++) {
  //  printf("indivs[%d], value=%3.2f, rank = %f, dataset=%d\n", i, indivs[i].value , indivs[i].rank, indivs[i].dataset);
  //}
  
  
  
  for (i=0; i<n; i++) {

    if (indivs[i].dataset == 1)
      W1 += indivs[i].rank;
    else 
      W2 += indivs[i].rank;
  }
  
  //printf("W1 = %f, W2 = %f\n", W1, W2);

  U1 = W1 - n1 * ( n1 + 1) / 2.0;
  U2 = W2 - n2 * ( n2 + 1) / 2.0;
  
  //printf("U1 = %3.2f, U2 = %3.2f\n", U1, U2);
  

  //
  //  calculate mean and standard deviation
  //
  mu = ( n1 * n2 ) / 2;
  si = n1 * n2 * ( n1 + n2 + 1 ) / 12;

  zi = ( U1 - mu ) / sqrt(si);

  //printf("z = %3.2f\n", zi);


  if (t == 0) {
    zi = ( min(U1, U2) - mu ) / sqrt(si);
  } else if (t == 1) {
    zi = ( U2 - mu ) / sqrt(si);
  } else if (t == 2) {
    zi = ( U1 - mu ) / sqrt(si);
  }
  
  return normal(zi);
  
}

//
//   compare two correlation values 
//
int CmpFuncIndividuals(const void* _a, const void* _b)
{
  const Individual* a = (const Individual*) _a;
  const Individual* b = (const Individual*) _b;
  
  if (a->value > b->value) 
    return 1; 
  else if(a->value == b->value) 
    return  0;
  else         
    return -1;
}


float euclidean(float* data1, float* data2, int m) 
{
  
  float sum   = 0.0;
  int   i;

  for (i=0; i<m; i++) {
    
    sum = sum + (data1[i] - data2[i]) * (data1[i] - data2[i]);
    //printf(" sum = %4.3f\n", sum);
  }
  return sqrt(sum);
}



float pearson_int(int* data1, int* data2, int m) 
{
  
  float sum1   = 0.0;
  float sum2   = 0.0;
    
  float sum1_2 = 0.0;
  float sum2_2 = 0.0;

  float avg1   = 0.0;
  float avg2   = 0.0;
    
  float a, b, c;

  int i;

  for (i=0; i<m; i++) {

    sum1   += data1[i];
    sum2   += data2[i];
    
    sum1_2 += data1[i] * data1[i];
    sum2_2 += data2[i] * data2[i];
    
  }
    
  avg1 = sum1 / (float)m;
  avg2 = sum2 / (float)m;

  // calc the Pearson correlation
    
  a = 0.0;
  b = 0.0;
  c = 0.0;

  for (i=0; i<m; i++) {
    a += (data1[i] - avg1) * (data2[i] - avg2);
    b += (data1[i] - avg1) * (data1[i] - avg1);
    c += (data2[i] - avg2) * (data2[i] - avg2);
  }
  
  if (a < DBL_EPSILON)
    return 0.0;
  else    
    return a / sqrt ( b * c );

}

float pearson(float* data1, float* data2, int m) 
{
  
  float sum1   = 0.0;
  float sum2   = 0.0;
    
  float sum1_2 = 0.0;
  float sum2_2 = 0.0;

  float avg1   = 0.0;
  float avg2   = 0.0;
    
  float a, b, c;

  int i;

  for (i=0; i<m; i++) {

    sum1   += data1[i];
    sum2   += data2[i];
    
    sum1_2 += data1[i] * data1[i];
    sum2_2 += data2[i] * data2[i];
    
  }
    
  avg1 = sum1 / m;
  avg2 = sum2 / m;

  // calc the Pearson correlation
    
  a = 0.0;
  b = 0.0;
  c = 0.0;

  for (i=0; i<m; i++) {
    a += (data1[i] - avg1) * (data2[i] - avg2);
    b += (data1[i] - avg1) * (data1[i] - avg1);
    c += (data2[i] - avg2) * (data2[i] - avg2);
  }
  
  //printf("a=%f, b=%f, c=%f, pe=%f\n", a, b, c, a / sqrt ( b * c ));

  if (a < DBL_EPSILON)
    return 0.0;
  else    
    return a / sqrt ( b * c );

}



//
//  Normal approximation
//  returns P(-inf<=x<=z)
//
double normal(const double z) 
{
 
  double b1 = 0.31938153;
  double b2 = -0.356563782;
  double b3 = 1.781477937;
  double b4 = -1.821255978;
  double b5 = 1.330274429;
  double p  = 0.2316419;
  double c2 = 0.3989423;
  double a  = fabs(z);
  double t = 1.0/(1.0+a*p);
  double b = c2*exp((-z)*(z/2.0));
  double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
  
   if (z > 6.0) { return 1.0; }; // this guards against overflow
   if (z < -6.0) { return 0.0; };
  n = 1.0-b*n;

  
  if ( z < 0.0 ) 
    n = 1.0 - n;
  
  
  return n;
}


double lcumhyper(int i, int s1, int s2, int N) {
  
  return log(cumhyper(i, s1, s2, N));

}


double cumhyper(int i, int s1, int s2, int N) {
  
  int min = (s1<s2?s1:s2);
  double prod = 0.0;
  int a;
  double tmp = 0.0;

  for (a=i; a<=min; a++) {
    tmp = (double)hypergeom(a, s1, s2, N);
    prod += (double)hypergeom(a, s1, s2, N);
  }
  
  return prod;
}


double cumhyper_u(int i, int s1, int s2, int N) {
  

  double prod = 0.0;
  int a;
  double tmp = 0.0;

  for (a=0; a<=i; a++) {
    tmp = (double)hypergeom(a, s1, s2, N);
    prod += (double)hypergeom(a, s1, s2, N);
  }
  
  return prod;
}



void nrerror(char str[]) {
  printf("%s\n", str);
}


double hypergeom(int i, int s1, int s2, int N) {
  double factln(int n);
    
  return exp(factln(s1) + factln(N - s1) + factln(s2) + factln(N - s2) - 
	     factln(i) - factln(N) - factln(s1 - i) - factln(s2 - i) 
	     - factln(N - s1 - s2 + i));
    
}



double factln(int n) {
  double gammln(double xx);
  void nrerror(char error_text[]);
  static double a[101];                                        
  if (n < 0) nrerror("Negative factorial in routine factln");
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));                                    
  else return gammln(n+1.0);                                  
}



double gammln(double xx) {
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}



double factrl(int n)  { 
  double gammln(double xx); 
  void nrerror(char error_text[]); 
  static int ntop=4; 
  static double a[33]={1.0,1.0,2.0,6.0,24.0}; 
  int j; 
  
  if (n < 0) 
    nrerror("Negative factorial in routine factrl"); 
  if (n > 32) 
    return exp(gammln(n+1.0)); 
  while (ntop<n) { 
    j=ntop++; 
    a[ntop]=a[j]*ntop; 
  } 
  
  return a[n]; 

}


double cumbino(int k, int N, double p) 
{

  int x;
  double pb;
  
  pb = 0;
  for (x=k; x<=N; x++) {
    pb += bico(N, x) * pow(p, x) * pow(1.0-p, N-x);
    
  }
  
  return pb;
}


double bico(int n, int k)                      
     //Returns the binomial coefficient nk as a doubleing-point number.
{
     double factln(int n);

     return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
     //The floor function cleans up roundoff error for smaller values of n and k.
}


double lbico(int n, int k)                      
     //Returns the binomial coefficient nk as a doubleing-point number.
{
     double factln(int n);

     return factln(n)-factln(k)-factln(n-k);
     //The floor function cleans up roundoff error for smaller values of n and k.
}


double lcumbino(int k, int N, double p) 
{
  int x;
  double pb;
  
  pb = 0;
  for (x=k; x<=N; x++) {
    pb += exp ( lbico(N, x) + log(pow(p, x)) + log(pow(1.0-p, N-x)));
  }
  
  return pb;
}


double lbino(int k, int N, double p) 
{
  double pb;    
  pb = exp ( lbico(N, k) + log(pow(p, k)) + log(pow(1.0-p, N-k)));    
  return pb;
}



void  getBiPermutations(int n1, int n2, int*** perms1, int*** perms2, int* nc) 
{
  
  int  i,j,n;
  int* ref1;
  int* ref2;
  
  // 
  //   fill in the reference vectors
  //
  ref1 = (int*)malloc(n1 * sizeof(int));
  ref2 = (int*)malloc(n2 * sizeof(int));

  for (i=0; i<n1; i++) {
    ref1[i] = i;
  }

  for (i=n1,j=0; i<n1+n2; i++,j++) {
    ref2[j] = i;
  }



  //
  // allocate memory for the shuffles
  //
  n = bico((n1+n2), n1);
  (*perms1) = (int**)malloc(n * sizeof(int*));
  (*perms2) = (int**)malloc(n * sizeof(int*));
    
  for (i=0; i<n2; i++) 
    _recBiPermutations(0, n1, i, n2, &ref1, &ref2, &ref1, &ref2, perms1, perms2, nc);
  
    
  
  
}



void _recBiPermutations(int ic, int n1, int jc, int n2, int** perm1, int** perm2, int** ref1, int** ref2, int*** perms1, int*** perms2, int* nc) 
{
  
  //printf("swap s1[%d] and s2[%d]\n");
  int  i, j;

  int* myperm1;
  int* myperm2;

  //printf("ic=%d, jc=%d\n", ic, jc);

  if (jc == n2) {

    // copy the permutation
    
    

    
    return;

  } else {
    
    for (i=ic; i<n1; i++) {

      //printf("swap s1[%d]=%d and s2[%d]=%d\n", i,  (*perm1)[i], jc, (*ref2)[jc]);
      
      // get copies of the current permutations
      myperm1 = (int*)malloc(n1 * sizeof(int));
      myperm2 = (int*)malloc(n2 * sizeof(int));
      //memcpy(myperm1, *perm1, sizeof(n1 * sizeof(int)));
      for (j=0; j<n1; j++) {
	myperm1[j] = (*perm1)[j];
      }
      for (j=0; j<n2; j++) {
	myperm2[j] = (*perm2)[j];
      }

      myperm1[i]  = (*ref2)[jc];
      myperm2[jc] = (*ref1)[i];

      /**
      printf("FINAL ");
    
      for (j=0; j<n1; j++) {
	printf("%d ", myperm1[j]);
      }
      printf("/ ");
      
      for (j=0; j<n2; j++) {
	printf("%d ", myperm2[j]);
      }
      printf("\n\n");
      **/

      (*perms1)[*nc] = (int*)malloc(n1 * sizeof(int));
      for (j=0; j<n1; j++) {
	(*perms1)[*nc][j] = myperm1[j];
      }
      (*perms2)[*nc] = (int*)malloc(n2 * sizeof(int));
      for (j=0; j<n2; j++) {
	(*perms2)[*nc][j] = myperm2[j];
      }
      (*nc)++;

	    
      for (j=jc+1; j<n2; j++) {
	_recBiPermutations(i+1, n1, j, n2, &myperm1, &myperm2, ref1, ref2, perms1, perms2, nc);
      }
    }

  }
  
}





float* getBootstrappedVector(float* sample, int n) 
{
  int    i;
  float* bt = (float*)malloc(n * sizeof(float));
  
  for (i=0; i<n; i++) {
    bt[i] =  sample[ (int) ((double)n*rand()/(RAND_MAX+1.0)) ];
  }

  return bt;


}


// 
int* shuffle(int* v, int n) 
{
  
  int* vnew;
  int  r, i;

  vnew = (int*)malloc(n * sizeof(int));

  for (i=0; i<n; i++) {
    //r = (int)((double)i*rand()/(RAND_MAX+1.0));
    r       = (int)((double)i*default_rand());
    //if ((r < 0) || (r >= n)) {
    //  die("Something went wrong with the default_rand() random number generator\n");
    //}
    vnew[i] = vnew[r]; 
    vnew[r] = v   [i]; 
  }

  return vnew;
}


int* shuffleInt(int* v, int n) 
{
  
  int* vnew;
  int  r, i;

  vnew = (int*)malloc(n * sizeof(int));

  for (i=0; i<n; i++) {
    r       = (int)((double)i*default_rand());
    //if ((r < 0) || (r >= n)) {
    //  die("Something went wrong with the default_rand() random number generator\n");
    //}
    vnew[i] = vnew[r]; 
    vnew[r] = v   [i]; 
  }

  return vnew;
}


float* shuffleFloat(float* v, int n) 
{
  
  float* vnew;
  int  r, i;

  vnew = (float*)malloc(n * sizeof(float));

  for (i=0; i<n; i++) {
    //r = (int)((double)i*rand()/(RAND_MAX+1.0));
    r       = (int)((double)i*default_rand());

    vnew[i] = vnew[r]; 
    vnew[r] = v   [i]; 
  }

  return vnew;
}





int* getIndices(int n) {
  int* vnew;
  int  i;

  vnew = (int*)malloc(n * sizeof(int));
  
  for (i=0; i<n; i++) {
    vnew[i] = i;
  }
  
  return vnew;
}


//
//  get a random permutation of n indices
//
int* randPermOfIndices(int n) 
{
  int* v1;
  int* v2;

  v1 = getIndices(n);
  
  v2 = shuffle(v1, n);
  
  free(v1);

  return v2;
  
}


void fitDirectLinearModel(float* x, float* y, int n, float* b) 
{
  
  //f = Sum_i (y_i - bx_i)2
  //df/db = 2 * (y_i - bx_i) ( -x_i ) = Sum_i - 2.y_i.x_i + 2b^2x_i^2 = 0 = Sum_i [ -y_i + bx_i ] = 0 =  Sum_i [ -y_i ] + b Sum_i [ x_i ] = 0 

  // df/dy = sum 2 * (y - bx) ( 1 ) = sum y - b sum x
  float sxy = 0.0;
  float sxx = 0.0;
  int   i;

  for (i=0; i<n; i++) {
    sxy += x[i] * y[i];
    sxx += x[i] * x[i];
  }

  *b = sxy/sxx;
    
}


void fitSimpleLinearModel(float* x, float* y, int n, float* b, float* a) 
{

  float ssxx;
  float ssxy;
  float ssyy;

  float sx2 = 0.0;
  float sy2 = 0.0;
  float sxy = 0.0;

  float xav = 0.0;
  float yav = 0.0;
  
  int    i;

  for (i=0; i<n; i++) {

    //printf("DATA %f, %f\n", x[i], y[i]); 

    sx2 += x[i] * x[i];
    sy2 += y[i] * y[i];
    sxy += x[i] * y[i];
  }

  //printf("sx2=%3.2f, sy2=%3.2f, sxy=%3.2f\n", sx2, sy2, sxy); 
  

  xav = average(x, n);
  yav = average(y, n);
  
  ssxx = sx2 - n * xav * xav;
  ssyy = sy2 - n * yav * yav;
  ssxy = sxy - n * xav * yav;

  //printf("ssxx=%3.2f, ssyy=%3.2f, ssxy=%3.2f\n", ssxx, ssyy, ssxy); 
  
  *b = ssxy / ssxx;
  *a = yav - (*b) * xav;
 
  return;
}


int maxArrayIndex (int* a, int n) 
{
  int i;
  int amax = -100000;  // dangerous
  int imax;
  
  
  for (i=0; i<n; i++) {
    if (a[i] > amax) {
      amax = a[i];
      imax = i;
    } 
  }
  
  return i;
}


void bubbleSort(int *a, int n)
{
  int i, j;
  float ftmp;


  for (i=0; i<n-1; i++) {
    for (j=0; j<n-1-i; j++)
      if (a[j+1] < a[j]) {  /* compare the two neighbors */

        ftmp = a[j];         /* swap a[j] and a[j+1]      */
        a[j] = a[j+1];
        a[j+1] = ftmp;

      }
  }

}



//
//  get order of integer
//
int* bubbleSortIndexInt(int *ia, int n)
{
  int i, j;

  int* a;
  int* b;


  // for storing temporary values
  int itmp;
  int ftmp;

  a = (int*)malloc(n * sizeof(int));
  b = (int*)malloc(n * sizeof(int));

  for (i=0; i<n; i++) {
    a[i] = ia[i];
    b[i] = i;
  }

  for (i=0; i<n-1; i++) {
    for (j=0; j<n-1-i; j++)
      if (a[j+1] < a[j]) {  /* compare the two neighbors */

        ftmp = a[j];         /* swap a[j] and a[j+1]      */
        a[j] = a[j+1];
        a[j+1] = ftmp;

        itmp = b[j];         /* swap a[j] and a[j+1]      */
        b[j] = b[j+1];
        b[j+1] = itmp;

      }
  }

  free(a);

  return b;

}


//
//  get order of floats
//
int* bubbleSortIndex(float *ia, int n) 
{
  int i, j;

  float* a;
  int* b;


  // for storing temporary values
  int   itmp;
  float ftmp;
  
  a = (float*)malloc(n * sizeof(float));
  b = (int*)malloc(n * sizeof(int));
  
  for (i=0; i<n; i++) {
    a[i] = ia[i];
    b[i] = i;
  }

  for (i=0; i<n-1; i++) {
    for (j=0; j<n-1-i; j++)
      if (a[j+1] < a[j]) {  /* compare the two neighbors */
        
        ftmp = a[j];         /* swap a[j] and a[j+1]      */
        a[j] = a[j+1];
        a[j+1] = ftmp;

        itmp = b[j];         /* swap a[j] and a[j+1]      */
        b[j] = b[j+1];
        b[j+1] = itmp;

      }
  }

  free(a);

  return b;

}



int CmpFloat(const void* _a, const void* _b)
{
  const float a = * ((const float*) _a);
  const float b = * ((const float*) _b);
  
  if (a > b) 
    return 1; 
  else if(a == b) 
    return  0;
  else         
    return -1;
}


int CmpFloatRegular(const void* _a, const void* _b)
{
  const float a = * ((const float*) _a);
  const float b = * ((const float*) _b);
  
  if (a < b) 
    return -1; 
  else if(a == b) 
    return  0;
  else         
    return 1;
}



int CmpDblRegular(const void* _a, const void* _b)
{
  const double a = * ((const double*) _a);
  const double b = * ((const double*) _b);
  
  if (a < b) 
    return -1; 
  else if(a == b) 
    return  0;
  else         
    return 1;
}





static inline char	*med3 (char *, char *, char *, int (*)());
static inline void	 swapfunc(char *, char *, int, int);


/*
 * Qsort routine from Bentley & McIlroy's "Engineering a Sort Function".
 */
#define swapcode(TYPE, parmi, parmj, n) { 		\
	long i = (n) / sizeof (TYPE); 			\
	register TYPE *pi = (TYPE *) (parmi); 		\
	register TYPE *pj = (TYPE *) (parmj); 		\
	do { 						\
		register TYPE	t = *pi;		\
		*pi++ = *pj;				\
		*pj++ = t;				\
        } while (--i > 0);				\
}

#define SWAPINIT(a, es) swaptype = ((char *)a - (char *)0) % sizeof(long) || \
	es % sizeof(long) ? 2 : es == sizeof(long)? 0 : 1;

static inline void swapfunc( char *a, char *b, int n, int swaptype)
{
  if (swaptype <= 1) 
    swapcode(long, a, b, n)
  else
    swapcode(char, a, b, n)
}

#define swap(a, b)					\
  if (swaptype == 0) {					\
    long t = *(long *)(a);				\
    *(long *)(a) = *(long *)(b);			\
    *(long *)(b) = t;					\
  } else						\
    swapfunc(a, b, es, swaptype)

#define vecswap(a, b, n) 	if ((n) > 0) swapfunc(a, b, n, swaptype)

static inline char *med3 (char *a, char *b, char *c, int (*cmp)())
{
  return cmp(a, b) < 0 ?
    (cmp(b, c) < 0 ? b : (cmp(a, c) < 0 ? c : a ))
    :(cmp(b, c) > 0 ? b : (cmp(a, c) < 0 ? a : c ));
}

void myqsort (void *a, unsigned long n, unsigned long es, int (*cmp)())
{
  char *pa, *pb, *pc, *pd, *pl, *pm, *pn;
  int d, r, swaptype, swap_cnt;
  
 loop:	SWAPINIT(a, es);
  swap_cnt = 0;
  if (n < 7) {
    for (pm = (char *) a + es; pm < (char *) a + n * es; pm += es)
      for (pl = pm; pl > (char *) a && cmp(pl - es, pl) > 0;
	   pl -= es)
	swap(pl, pl - es);
    return;
  }
  pm = (char *) a + (n / 2) * es;
  if (n > 7) {
    pl = a;
    pn = (char *) a + (n - 1) * es;
    if (n > 40) {
      d = (n / 8) * es;
      pl = med3(pl, pl + d, pl + 2 * d, cmp);
      pm = med3(pm - d, pm, pm + d, cmp);
      pn = med3(pn - 2 * d, pn - d, pn, cmp);
    }
    pm = med3(pl, pm, pn, cmp);
  }
  swap(a, pm);
  pa = pb = (char *) a + es;
  
  pc = pd = (char *) a + (n - 1) * es;
  for (;;) {
    while (pb <= pc && (r = cmp(pb, a)) <= 0) {
      if (r == 0) {
	swap_cnt = 1;
	swap(pa, pb);
	pa += es;
      }
      pb += es;
    }
    while (pb <= pc && (r = cmp(pc, a)) >= 0) {
      if (r == 0) {
	swap_cnt = 1;
	swap(pc, pd);
	pd -= es;
      }
      pc -= es;
    }
    if (pb > pc)
      break;
    swap(pb, pc);
    swap_cnt = 1;
    pb += es;
    pc -= es;
  }
  if (swap_cnt == 0) {  /* Switch to insertion sort */
    for (pm = (char *) a + es; pm < (char *) a + n * es; pm += es)
      for (pl = pm; pl > (char *) a && cmp(pl - es, pl) > 0; 
	   pl -= es)
	swap(pl, pl - es);
    return;
  }
  
  pn = (char *) a + n * es;
  r = min(pa - (char *)a, pb - pa);
  vecswap(a, pb - r, r);
  r = min(pd - pc, pn - pd - es);
  vecswap(pb, pn - r, r);
  if ((r = pb - pa) > es)
    qsort(a, r / es, es, cmp);
  if ((r = pd - pc) > es) { 
    /* Iterate rather than recurse to save stack space */
    a = pn - r;
    n = r / es;
    goto loop;
  }
  /*		qsort(pn - r, r / es, es, cmp);*/
}
