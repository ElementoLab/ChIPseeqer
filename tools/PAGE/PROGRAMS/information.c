#include "information.h"
#include "statistics.h"

#define _GNU_SOURCE
#ifndef CPROTO
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#endif

#ifdef USEMPATROL
#include <mpatrol.h>
#endif


#define REPEATS 1




int isDiagonalDominant(int* A_q, int* B_q, int DA, int DB, int n) 
{
  
  int** tmpCounts;
  int   i, j;
  int   c ;

  tmpCounts   = ConstructCountMatrix (A_q, B_q, n, DA, DB);
  //showCountMatrix(tmpCounts, DA, DB, 0, 0, 0);

  for (i=0; i<DA; i++) {
    c = tmpCounts[i][i];
    for (j=0; j<DB; j++) {
      if (i == j)
	continue;
      if (tmpCounts[i][j] > c)
	return 0;
    }
  }
  free(tmpCounts);
  return 1;

}


//
//  calculate the increase in information I(A,E;B) - I(A;B) >= 0
//
double InformationIncrease(int* A_q, int DA, int* B_q, int DB, int* E_q, int DE, int n, double* I1, double* I2) 
{
  
  int*   AE_q;
  int    DAE;
  double mI1;
  double mI2;

  DAE  = DA * DB;
  AE_q = combineQuantizedVectors(A_q,  E_q, n, DA, DE);
  mI1   = CalculateMIbasic       (AE_q, B_q, n, DAE, DB);
  mI2   = CalculateMIbasic        (A_q, B_q, n, DA , DB);

  *I1 = mI1;
  *I2 = mI2;

  return mI1 - mI2;
  
}


float FindInfo(int* qv1, int* qv2, int FullSample, int D1, int D2, float* SampleSizesPerc, int N, int* Trials, int Repeats, float* std) {
  
  int  i;
  int  r;
  int  t;
  int  n;
  int  tmpSampleSize;

  int* SampleSizes; 
  float* SampleSizesInv;
  int* perm;
  int* tmpqv1;
  int* tmpqv2;

  float** Ivals;
  float*  Extrap_I;

  int**  tmpCounts;

  float  a,b;
  
  float* meanIvals;
  float  mi;

  int    trialMax;


  if ((D1 == 1) || (D2 == 1)) {
    *std = 0.0;
    return 0.0;
  }

  
  Extrap_I    = (float*)malloc(Repeats * sizeof(float));
  SampleSizes = (int*)  malloc(N * sizeof(int));
  SampleSizesInv = (float*)  malloc(N * sizeof(float));

  meanIvals   = (float*)malloc(N * sizeof(float));
  Ivals       = (float**)malloc(N * sizeof(float*));
  
  // max number of trials

  trialMax    =  maxArrayIndex(Trials, N);
  //for (i=0; i<N; i++) {
  //  Ivals[i] = (float*)malloc(N * sizeof(float));
  //}

  //printf("D1=%d, D2=%d\n", D1, D2);

  for (i=0; i<N; i++) {
    SampleSizes[i] = ceil(FullSample * SampleSizesPerc[i]);
    SampleSizesInv[i] = 1.0 / SampleSizes[i];
    Ivals[i] = (float*)malloc(Trials[i] * sizeof(float));
  }
  
  
  for (r=0; r<Repeats; r++) {
    
    for (n=0; n<N; n++) {

      //printf("sample size n = %d\n", n);

      tmpSampleSize = SampleSizes[n];
      
      tmpqv1        = (int*)malloc(tmpSampleSize * sizeof(int));
      tmpqv2        = (int*)malloc(tmpSampleSize * sizeof(int));
      

      for (t=0; t<Trials[n]; t++) {
	
	//printf("Trial t = %d\n", t);
	
	if (tmpSampleSize < FullSample) {
	  
	  // create a shuffled list of indices
	  perm   = randPermOfIndices(FullSample);
	  
	  // create the subvectors
	  for (i=0; i<tmpSampleSize; i++) {
	    tmpqv1[i] = qv1[ perm[i] ];
	    tmpqv2[i] = qv2[ perm[i] ];
	  }
	  
	  free(perm);

	  // 
	  
	} else {
	  
	  // create the subvectors
	  for (i=0; i<tmpSampleSize; i++) {
	    tmpqv1[i] = qv1[i];
	    tmpqv2[i] = qv2[i];
	  }
	  
	}
	
	//printf("Constructing count matrix\n");
	tmpCounts   = ConstructCountMatrix (tmpqv1, tmpqv2, tmpSampleSize, D1, D2);

	//showCountMatrix(tmpCounts, D1, D2); getchar();

	//printf("Local MI\n");
	Ivals[n][t] = FindInfoLocalMI(tmpCounts, D1, D2);

	//printf("Intermediate MI[%d][%d] = %3.2f\n", n, t, Ivals[n][t]); 
	//printf("freeing matrix\n");
	freeCountMatrix(tmpCounts, D1, D2);
	free(tmpCounts);
	
      }

      free(tmpqv1);
      free(tmpqv2);

      
      // mean over all trials
      meanIvals[n] = 0.0;
      for (t=0; t<Trials[n]; t++) {
	meanIvals[n] += Ivals[n][t];
      }
      meanIvals[n] /= Trials[n];
      

    }
    //
    //  END sample size loop
    //
    
    fitSimpleLinearModel( SampleSizesInv, meanIvals, N, &b, &a);

    //printf("Output of linear model, b=%3.2f, a=%3.2f\n", b, a);
    
    Extrap_I[r] = a;

    //printf("done with repeat %d\n", r);
   
    
  }
  //
  //  END repeat loop
  //

  // return the average over all the repeats
  mi = 0.0;
  for (i=0; i<Repeats; i++) { 
    //printf("mi[r=%d]=%f\n", i, Extrap_I[i]);
    mi += Extrap_I[i];
  }
  mi = mi/Repeats;

  if (Repeats < 5) 
    *std = stddev(Ivals[0], Trials[0]);
  else {
    *std = stddev(Extrap_I, Repeats); 
  }



  free(Extrap_I);
  free(SampleSizes);
  free(SampleSizesInv);
  free(meanIvals);
  for (i=0; i<N; i++) {
    free(Ivals[i]);
  }
  free(Ivals);
  
  

  return mi;
}




int** Construct3DCountMatrix(int* qv1, int* qv2, int* qv3, int m, int D1, int D2, int D3) {
  
  int   i, j, k;
  int*** counts;

  counts = (int***)malloc(D1 * sizeof(int**));
  if (counts == 0) {
    printf("Construct3DCountMatrix: not enough memory\n");
    exit(0);
  }

  for (i=0; i<D1; i++) {
    counts[i] = (int**)malloc(D2 * sizeof(int*));
    if (counts[i] == 0) {
      printf("Construct3DCountMatrix: not enough memory\n");
      exit(0);
    }

    for (j=0; j<D2; j++) {
      counts[i][j] = (int*)malloc(D3 * sizeof(int));
      if (counts[i][j] == 0) {
	printf("Construct3DCountMatrix: not enough memory\n");
	exit(0);
      }

      for (k=0; k<D3; k++) {
	counts[i][j][k] = 0;
      }
      
    }
    
  }

  for (i=0; i<m; i++) {
    if ((qv1[i] >= D1) || (qv2[i] >= D2) || (qv3[i] >= D3)) {
      printf("Construct3DCountMatrix: something is wrong with the quantization %d >= D1=%d or %d >= D2=%d\n", qv1[i], D1, qv2[i], D2);
      exit(0);
    }
    
    counts[ qv1[i] ] [ qv2[i] ] [ qv3[i] ] ++; 
  }
  
  return counts;
}



int** ConstructCountMatrix(int* qv1, int* qv2, int m, int D1, int D2) {
  
  int   i, j;
  int** counts;

  counts = (int**)malloc(D1 * sizeof(int*));
  if (counts == 0) {
    printf("ConstructCountMatrix: not enough memory\n");
    exit(0);
  }
  for (i=0; i<D1; i++) {
    counts[i] = (int*)malloc(D2 * sizeof(int));
    if (counts[i] == 0) {
      printf("ConstructCountMatrix: not enough memory\n");
      exit(0);
    }
    for (j=0; j<D2; j++) {
      counts[i][j] = 0;
    }
  }

  for (i=0; i<m; i++) {
    if ((qv1[i] >= D1) || (qv2[i] >= D2)) {
      printf("something is wrong with the quantization %d >= D1=%d or %d >= D2=%d\n", qv1[i], D1, qv2[i], D2);
      exit(0);
    }
    counts[ qv1[i] ] [ qv2[i] ] ++; //printf("%d\n", counts[ qv1[i] ] [ qv2[i] ]);
  }
  
  return counts;
}


void showCountMatrix(int** counts, int D1, int D2, char** states1, char** states2, char* prefix) 
{
  
  int   i, j;
  
  if (states2 != NULL) {
    if (prefix != NULL) {
      printf("%s", prefix);
    }
    for (j=0; j<D2; j++) {
      printf("\t%s", states2[j]);
    }
    printf("\n");
  }

  for (i=0; i<D1; i++) {
    
    printf("S%d ", i);

    if (prefix != NULL) {
      printf("%s", prefix);
    }
    if (states1 != NULL) {
      printf("%s", states1[i]);
      printf("\t");
    }

    for (j=0; j<D2; j++) {
      printf("%d\t", counts[i][j]);
    }

    printf("\n");
  }
  
}


void freeCountMatrix(int** counts, int D1, int D2) {
  int i;

  for (i=0; i<D1; i++) {
    free(counts[i]);
  }
  //free(counts);
}


//
//  KL on direct discrete linearized distributions
//
float KL_dist(float* P1, float* P2, int D) 
{
  float e  = 2.204e-16;
  float KL = 0.0;
  int   i;

  for (i=0; i<D; i++) {
    if ( P1[i] >= e )
      KL += P1[i] * log ( P1[i]/P2[i]) / log(2.0);
  }

  return KL;
 
}

//
//  entropy
//
float entropy_from_probas(float* P, int n) 
{
  
  int i;
  float e = 0.0;

  for (i=0; i<n; i++) {
    //printf("%5.4f\n", P[i]);
    if (P[i] > 0.0) {
      e += - P[i] * log( P[i] ) / log(2.0);
    }
  }
  
  return e;
}


//
//    Jensen-Shannon divergence
//
float JS_dist(float* P1, float* P2, int D) 
{

  float  js;
  float* PS;
  int    i;

  PS = (float*)malloc(D * sizeof(float));
  
  for (i=0; i<D; i++) {
    PS[i] = 0.5 * P1[i] + 0.5 * P2[i];
  }

  js = 0.5 * KL_dist(P1, PS, D) + 0.5 * KL_dist(P2, PS, D);

  free(PS);
  
  return js;
  
}


float CalculateCondMIbasic( int* data1,  int* data2, int* data3, int m, int D1, int D2, int D3) 
{
  
  float   mi;
  int i, j;
  int***   tmpCounts;
 
  tmpCounts   = Construct3DCountMatrix (data1, data2, data3, m, D1, D2, D3);
  mi          = FindCondMI(tmpCounts, D1, D2, D3);

  for (i=0; i<D1; i++) {
    for (j=0; j<D2; j++) {
      free(tmpCounts[i][j]);
    }
    free(tmpCounts[i]);
  }
  free(tmpCounts);
  
  return mi;

}

//
//
//
double FindCondMI(int*** counts, int D1, int D2, int D3) 
{

  double CMI   = 0.0;
  double pz, pyz, pxz, pxyz; 

  int   x, y, z;

  int   n = 0;
  double    e = 0.0; //2.204e-16;

  int*  counts_z;
  int** counts_xz;
  int** counts_yz;


  for (x=0; x<D1; x++) {
    for (y=0; y<D2; y++) {
      for (z=0; z<D3; z++) {
	n += counts[x][y][z];
      }
    }
  }

  
  //
  // marginal distribution p(z)
  //
  counts_z = (int*)malloc(D3 * sizeof(int));
  for (z=0; z<D3; z++) {
    counts_z[z] = 0;
    for (x=0; x<D1; x++) {
      for (y=0; y<D2; y++) {
	counts_z[z] += counts[x][y][z];
      }
    }
  }

  //
  // marginal distribution p(x,z)
  //  
  counts_xz = (int**)malloc(D1 * sizeof(int*));
  for (x=0; x<D1; x++) {
    counts_xz[x] = (int*)malloc(D3 * sizeof(int));
  }
  for (x=0; x<D1; x++) {
    for (z=0; z<D3; z++) {
      counts_xz[x][z] = 0;
      for (y=0; y<D2; y++) {
	counts_xz[x][z] += counts[x][y][z];
      }
    }
  }

  
  //
  // marginal distribution p(y,z)
  //
  counts_yz = (int**)malloc(D2 * sizeof(int*));
  for (y=0; y<D2; y++) {
    counts_yz[y] = (int*)malloc(D3 * sizeof(int));
  }
  for (y=0; y<D2; y++) {
    for (z=0; z<D3; z++) {
      counts_yz[y][z] = 0;
      for (x=0; x<D1; x++) {
	counts_yz[y][z] += counts[x][y][z];
      }
    }
  }

  
  // 
  // conditional mi				       
  // 
  CMI = 0.0;
  
  for (x=0; x<D1; x++) {
    for (y=0; y<D2; y++) {
      for (z=0; z<D3; z++) {
	
	if ((counts[x][y][z] > 0) && (counts_z[z] > 0)) {

	  pz   = (e + (double)counts_z  [z]       ) / (double)n;
	  pxz  = (e + (double)counts_xz [x][z]    ) / (double)n;
	  pyz  = (e + (double)counts_yz [y][z]    ) / (double)n;
	  pxyz = (e + (double)counts    [x][y][z] ) / (double)n;
	  
	  CMI += pxyz * log ( pxyz * pz / ( pxz * pyz ) ) / log( 2.0 ); 
	  
	}

	
      }
    }
  }

  free(counts_z);
  

  for (y=0; y<D2; y++) {
    free(counts_yz[y]);
  }
  free(counts_yz);


  for (x=0; x<D1; x++) {
    free(counts_xz[x]);
  }
  free(counts_xz);


  return CMI;

}

//
//  calculate mi from count matrix
//
double FindInfoLocalMI(int** counts, int D1, int D2) {
  
  double MI = 0.0;
  int    x,y;
  int    n;
  double  px, pxy, py;
  int*   counts_x;
  int*   counts_y;
  //double    e = 2.204e-16;
  double      e = 0.0;

  //showCountMatrix(counts, D1, D2, NULL, NULL, NULL);

  
  n = 0;
  for (x=0; x<D1; x++) {
    for (y=0; y<D2; y++) {
      n += counts[x][y];
    }
  }

  if (n == 0) {
    return 0.0;
  }
  

  counts_x = (int*)calloc(D1, sizeof(int));
  counts_y = (int*)calloc(D2, sizeof(int));

  for (x=0; x<D1; x++) {
    counts_x[x] = 0;
    for (y=0; y<D2; y++) {
      counts_x[x] += counts [x][y];
      
    }
  }

  for (y=0; y<D2; y++) {
    counts_y[y] = 0;
    for (x=0; x<D1; x++) {
      counts_y[y] += counts [x][y];
    }
  }

  for (x=0; x<D1; x++) {
    
    for (y=0; y<D2; y++) {
      if (counts [x][y] == 0) {
	//MI += 0.0;
      } else {
	pxy = (e + (double)counts [x][y]) / (double)n;
      	px  = (e + (double)counts_x[ x ]) / (double)n;
	py  = (e + (double)counts_y[ y ]) / (double)n;
	MI += pxy * log ( pxy / (px * py) ) / log(2.0) ; 
	//printf("x=%d,y=%d, MI += pxy(%3.2f) * log( pxy(%3.2f) / (px(%3.2f) * py(%3.2f)) ) = %5.4f\n", x, y, pxy, pxy, px, py, MI);

      }

    }
  }

  free(counts_x);
  free(counts_y);

  
  //printf("MI=%5.4f\n", MI);

  return MI;
}

void InformationGain_D2(int** counts, int D1, int D2, MI_Contribs** mico, int sort) 
{
  
  //float         MI = 0.0;
  int           x,y;
  int           n;
  float         px, pxy, py;
  int*          counts_x;
  int*          counts_y;
  float         e = 2.204e-16;
  //float         mii;
  MI_Contribs*  my_mico;
  int           nm = 0;
  extern int    Cmp_MI_Contribs();
  float         px_bar_y;
  float         ig;
  float         co, best_co;
  int           best_x;

  my_mico = (MI_Contribs*) malloc( D2 * sizeof ( MI_Contribs ));

  n = 0;
  for (x=0; x<D1; x++) {
    for (y=0; y<D2; y++) {
      n += counts[x][y];
    }
  }

  if (n == 0) {
    return;
  }
  

  counts_x = (int*)calloc(D1, sizeof(int));
  counts_y = (int*)calloc(D2, sizeof(int));

  for (x=0; x<D1; x++) {
    counts_x[x] = 0;
    for (y=0; y<D2; y++) {
      counts_x[x] += counts [x][y];
      
    }
  }

  for (y=0; y<D2; y++) {
    counts_y[y] = 0;
    for (x=0; x<D1; x++) {
      counts_y[y] += counts [x][y];
    }
  }

  for (y=0; y<D2; y++) {

    ig = 0.0;
    py = (e + (float)counts_y[ y ]) / (float)n;
        
    
    best_co = -100;
    best_x  = -100;

    for (x=0; x<D1; x++) {

      px  = (e + (float)counts_x[ x ]) / (float)n;
      pxy = (e + (float)counts [x][y]) / (float)n;

      if (counts [x][y] == 0) {
	px_bar_y = 0.0;
	co = 0.0;
      } else {
	px_bar_y = pxy / py; 	// p(x,y) = p(y)p(x|y)
	co = px_bar_y * log ( px_bar_y / px ) / log(2.0) ; 
      }
      
      if (co > best_co) {
	best_co = co;
	best_x  = x;
      }

      ig += co;
    }

    ig = ig * py;

    
    my_mico[ nm ].i  =  y;
    my_mico[ nm ].j  =  best_x;
    my_mico[ nm ].co =  ig;
    nm ++;
  }
  
  if (sort == 1)
    qsort((void*)my_mico, D2, sizeof(MI_Contribs), Cmp_MI_Contribs);
  
  *mico = my_mico;

  free(counts_x);
  free(counts_y);

}

void StudyContribution_FindInfoLocalMI(int** counts, int D1, int D2, MI_Contribs** mico) 
{
  
  float MI = 0.0;
  int    x,y;
  int           n;
  float         px, pxy, py;
  int*          counts_x;
  int*          counts_y;
  float         e = 2.204e-16;
  float         mii;
  MI_Contribs*  my_mico;
  int           nm = 0;
  extern int Cmp_MI_Contribs();

  my_mico = (MI_Contribs*) malloc( D1 * D2 * sizeof ( MI_Contribs ));

  n = 0;
  for (x=0; x<D1; x++) {
    for (y=0; y<D2; y++) {
      n += counts[x][y];
    }
  }

  if (n == 0) {
    return;
  }
  

  counts_x = (int*)calloc(D1, sizeof(int));
  counts_y = (int*)calloc(D2, sizeof(int));

  for (x=0; x<D1; x++) {
    counts_x[x] = 0;
    for (y=0; y<D2; y++) {
      counts_x[x] += counts [x][y];
      
    }
  }

  for (y=0; y<D2; y++) {
    counts_y[y] = 0;
    for (x=0; x<D1; x++) {
      counts_y[y] += counts [x][y];
    }
  }

  for (x=0; x<D1; x++) {
    
    for (y=0; y<D2; y++) {

      if (counts [x][y] == 0) {
	mii = 0;
      } else {
	pxy = (e + (float)counts [x][y]) / (float)n;
      	px  = (e + (float)counts_x[ x ]) / (float)n;
	py  = (e + (float)counts_y[ y ]) / (float)n;
	mii = pxy * log ( pxy / (px * py) ) / log(2.0) ; 
	MI += mii;
      }

      my_mico[ nm ].i  = x;
      my_mico[ nm ].j  = y;
      my_mico[ nm ].co = mii;
      nm ++;
    }
  }

  qsort((void*)my_mico, D1*D2, sizeof(MI_Contribs), Cmp_MI_Contribs);
  
  *mico = my_mico;

  free(counts_x);
  free(counts_y);

}


int Cmp_MI_Contribs(const void* _a, const void* _b)
{
  const MI_Contribs* a = (const MI_Contribs*) _a;
  const MI_Contribs* b = (const MI_Contribs*) _b;
  
  if (a->co < b->co) 
    return 1; 
  else if(a->co == b->co) 
    return  0;
  else         
    return -1;
}


int** allocateQuantizedMatrix(int n) 
{
  
  int** qv = (int**)calloc(n, sizeof(int*));

  return qv;

}

float* sortFloatVector(float* v, int n)
{
  int    i;
  float* vc;
  vc    = (float*)malloc(n * sizeof(float));
  for (i=0; i<n; i++) {
    vc[i] = v[i];
  }
  qsort((void*)vc, n, sizeof(float), CmpDbl);
  
  return vc;

  
}


int* TopQuantize(float* v, int n, int ne, float* th) 
{
  int*   qv;
  int    i;
  float* vc;
  float   t;
  
  qv    = (int*)calloc(n, sizeof(int));

  vc    = (float*)malloc(n * sizeof(float));
  for (i=0; i<n; i++) {
    vc[i] = v[i];
  }

  qsort((void*)vc, n, sizeof(float), CmpDbl);

  t     = vc[n - ne]; //printf("th=%3.2f\n", t);
  *th   = t;

  for (i=0; i<n; i++) {
    if (v[i] >= t) {
      qv[i] = 1;
    } else {
      qv[i] = 0;
    }

  }

  free(vc);
  
  return qv;

}

int* ThresholdQuantizeInt(int* v, int n, int t) 
{
  int*   qv;
  int    i;

  qv    = (int*)malloc(n * sizeof(int));
  
  for (i=0; i<n; i++) {

    if (v[i] == INT_MAX) 
      qv[i] = INT_MAX;
    else {
      if (v[i] >= t) {
	qv[i] = 1;
      } else {
	qv[i] = 0;
      }
    }

  }
  
  return qv;

}

int* ThresholdQuantize(float* v, int n, float t, int sup) 
{
  int*   qv;
  int    i;

  qv    = (int*)malloc(n * sizeof(int));
  
  for (i=0; i<n; i++) {
    
    if (sup == 1) {
      if (v[i] >= t) {
	qv[i] = 1;
      } else {
	qv[i] = 0;
      }
    } else {
      if (v[i] <= t) {
	qv[i] = 1;
      } else {
	qv[i] = 0;
      }
    }

  }
  
  return qv;

}

int* QuantizeUsingThresholds(float* v, int n, int D, float** bins) 
{
  
  int    i, d, j;
  float* vc;
  int*   qv;
  
  int binsize   = (int)ceilf( n / (float)D);
  int n_nomissing = 0;

  int e = 0.0001;
  
  int verbose = 0;

  float* binup;

  //
  // copy the vector first
  //
  vc    = (float*)calloc(n,  sizeof(float));
  for (i=0; i<n; i++) {
    if (!isnan(v[i])) {
      vc[n_nomissing] = v[i];
      n_nomissing++;
    }
  }
  
  //memcpy(vc, v, n*sizeof(float));

  binup = (float*)malloc((D+1) * sizeof(float));
  qv    = (int*)calloc(n,  sizeof(int));

  //
  // sort the vector of floats 
  //

  qsort((void*)vc, n, sizeof(float), CmpDbl);

  for (i=0; i<=D; i++) {
    //printf("i=%d, i*binsize=%d\n", i, i*binsize);
    if (i * binsize < n)
      binup[i] = vc[i * binsize]; // binup[0] = vc[0], binup[1] = vc [ 1 * binsize ], etc 
    else {
      // printf("last value is %4.3f\n", vc[n-1]);
      binup[i] = vc[n-1]; // last value in the sorted array
      break;
    }
  }
  binup[0] -= e;
    
  //
  // now actually do the binning
  //
  for (i=0; i<n; i++) {
    
    if (isnan(v[i])) {
      qv[i] = INT_MAX;
    } else {

      //
      //  find the bin
      //
      
      //printf("v[i=%d] = %3.2f\n", i, v[i]);

      for (j=0; j<D; j++) {
	
	if (verbose == 1)
	  printf(" is v[i] between %3.2f and %3.2f\n", binup[j], binup[j+1]);
	
	if ((v[i] > binup[j]) && (v[i] <= binup[j+1])) {
	  
	  if (verbose == 1)
	    printf("  yes\n"); 
	  
	  qv[i] = j;
	  break;
	} else {
	  if (verbose == 1)
	    printf("no\n");
	}
	
      }
    }
    
  }

  binup[0] += e;

  if (bins != NULL)
    *bins = binup;
  else
    free(binup);

  free(vc);

  return qv;
}



int* Quantize(float* v, int n, int D, float** binup) 
{
  
  int    i, d, j;
  float* vc;
  int*   qv;
  FloatAndIndex* fi;

  int binsize;
  int toj;
  
  binsize = (int)nearbyintf ( n / (float)D);
  fi      = (FloatAndIndex*)malloc( n * sizeof( FloatAndIndex ) );
  for (i=0; i<n; i++) {
    fi[i].i = i;
    fi[i].v = v[i];
  }
  
  qsort((void*)fi, n, sizeof(FloatAndIndex), CmpFI);

  //for (i=0; i<n; i++) {
  //  printf("%d\t%f\n", fi[i].i, fi[i].v);
  //}
  

  qv      = (int*)calloc(n,  sizeof(int));
  //for (i=0; i<n; i++) {
  //  qv[i] = -1;
  //}
  
  if (binup != 0)
    *binup   = (float*)malloc((D+1) * sizeof(float));
  
  for (i=0; i<D; i++) {
    j        = i * binsize;
    if (binup != 0)
      (*binup)[i] = fi[j].v;
    
    toj = (i+1)*binsize;
    
    // very last bin
    if (i == D-1) {
      toj = n;
    }
    while ((j < toj) && (j < n)) {
      qv[ fi[j].i ] = i;
      j ++;
    }
    
  }

  if (binup != 0)
    (*binup)[D] = fi[ n-1 ].v;

  free(fi);

  return qv;
}


void QuantizePointerized(float* v, int n, int D, int** qv, float** bins) 
{
  
  int    i, d, j;
  float* vc;

  
  int binsize   = (int)ceilf( n / (float)D);
  int n_nomissing = 0;

  float e = 0.01;
  
  int verbose = 1;

  float* binup;

  //
  // copy the vector first
  //
  vc    = (float*)malloc(n * sizeof(float));
  for (i=0; i<n; i++) {
    if (!isnan(v[i])) {
      vc[n_nomissing] = v[i];
      n_nomissing++;
    }
  }
  
  //memcpy(vc, v, n*sizeof(float));
  
  binup = (float*)calloc((D+1), sizeof(float));
  *qv   = (int*)  malloc(n     * sizeof(int));

  //
  // sort the vector of floats 
  //

  qsort((void*)vc, n, sizeof(float), CmpDbl);

  for (i=0; i<n; i++) {
    printf("% 3.2f", vc[i]);
  }
  printf("\n");

  for (i=0; i<=D; i++) {
    printf("i=%d, i*binsize=%d\n", i, i*binsize);
    if (i * binsize < n) {
      binup[i] = vc[i * binsize]; 
      printf("binup[i=%d] = %3.2f (vc[%d])\n", i, vc[i * binsize], i*binsize);
    }
    // binup[0] = vc[0], binup[1] = vc [ 1 * binsize ], etc 
    else {
      printf("last value is %4.3f\n", vc[n-1]);
      binup[i] = vc[n-1]; // last value in the sorted array
      break;
    }
  }

  binup[0] -= e;
    
  //
  // now actually do the binning
  //
  for (i=0; i<n; i++) {
    
    if (isnan(v[i])) {
      (*qv)[i] = INT_MAX;
    } else {

      //
      //  find the bin
      //
      
      if (verbose == 1)  printf("v[i=%d] = %3.2f\n", i, v[i]);

      for (j=0; j<D; j++) {
	
	if (verbose == 1) printf(" is v[i]=%3.2f between %3.2f and %3.2f\n", v[i],binup[j], binup[j+1]);
	
	if ((v[i] > binup[j]) && (v[i] <= binup[j+1])) {
	  
	  if (verbose) printf("  yes, qv[i] is %d\n", j); 
	  
	  (*qv)[i] = j;
	  break;
	} else {
	  if (verbose == 1) printf("no\n");
	}
	
      }
    }
    
  }

  binup[0] += 3;

  if (bins != NULL)
    *bins = binup;
  else
    free(binup);

  free(vc);

  return; // qv;
}




float CalculateMI(int* data1,  int* data2, int m, int D1, int D2, int Repeats) 
{
  //int     Repeats = REPEATS;
  int     FullSample;
  float   SampleSizesPerc[4] = {0.7, 0.7875, 0.9, 1.0};
  int     Trials [4]         = { 21,     16,  12,   1};
  float   mi;
  float   std;
  int     N = 4;
  
  int     i;
  int**   tmpCounts;

  /*
  tmpCounts   = ConstructCountMatrix (data1, data2, m, D1, D2);
  mi = FindInfoLocalMI(tmpCounts, D1, D2);
  showCountMatrix(tmpCounts, D1, D2, NULL, NULL, NULL);
  for (i=0; i<D1; i++) {
    free(tmpCounts[i]);
  }
  free(tmpCounts);
  */
  mi = FindInfo(data1, data2, m, D1, D2, SampleSizesPerc, N, Trials, Repeats, &std); 

  return mi;
}



float CalculateMIbasic(int* data1,  int* data2, int m, int D1, int D2) 
{
  float   mi;
  int i;
  int**   tmpCounts;
 
  tmpCounts   = ConstructCountMatrix (data1, data2, m, D1, D2);
  mi = FindInfoLocalMI(tmpCounts, D1, D2);
  for (i=0; i<D1; i++) {
    free(tmpCounts[i]);
  }
  free(tmpCounts);
  return mi;
}


float CalculateNormalizedMIbasic(int* data1,  int* data2, int m, int D1, int D2) 
{
  float   mi;
  int     i;
  int**   tmpCounts;
  float   en;

  tmpCounts   = ConstructCountMatrix (data1, data2, m, D1, D2);
  mi = FindInfoLocalMI(tmpCounts, D1, D2);

  en = entropy(data1, m, D1); 

  for (i=0; i<D1; i++) {
    free(tmpCounts[i]);
  }
  free(tmpCounts);
  return mi / en;
}






//
//  calculate and display a count matrix 
//
void CalculateAndShowCountMatrix(int* data1,  int* data2, int m, int D1, int D2, char** states1, char** states2, char* prefix) 
{

  int     i;
  int**   tmpCounts;
  float   mymi;
 
  tmpCounts   = ConstructCountMatrix (data1, data2, m, D1, D2);
  showCountMatrix(tmpCounts, D1, D2, states1, states2, prefix);

  mymi = FindInfoLocalMI(tmpCounts, D1, D2);
  printf("mi=%3.2f\n", mymi);

  freeCountMatrix(tmpCounts, D1, D2);
  free(tmpCounts);
  
}



int* combineQuantizedVectors(int* data1,  int* data2, int m, int D1, int D2) 
{
  
  int* newdata;
  int  i, j;

  newdata = (int*)malloc(m * sizeof(int));
  if (newdata == 0) {
    printf("combineQuantizedVectors : not enough memory\n");
    exit(0);
  }
  for (i=0; i<m; i++) {
    if ((data1[ i ] == INT_MAX) || (data2[ i ] == INT_MAX))
      newdata [ i ] = INT_MAX;
    else
      newdata [ i ] = data1[ i ]  + data2[ i ] * D1;
  }


  return newdata;
  
}


int* combineBinaryVectorsAND(int* data1,  int* data2, int m) 
{
  
  int* newdata;
  int  i, j;

  newdata = (int*)malloc(m * sizeof(int));

  for (i=0; i<m; i++) {
    if ((data1[ i ] == 1) && (data2[ i ] == 1)) {
      newdata [ i ] = 1;
    } else {
      newdata [ i ] = 0;
    }
  
  }

  return newdata;
  
}



void  removeMissingValuesFloatFloat(float* data1, float* data2, 
				    int m, 
				    float** data1_nomissing, 
				    float** data2_nomissing,
				    int* m_nomissing) 
{
  int i;

  *data1_nomissing = (float*)malloc(m * sizeof(float));
  *data2_nomissing = (float*)malloc(m * sizeof(float));
  *m_nomissing = 0;

  for (i=0; i<m; i++) {
    if ((!isnan(data1[i])) && (!isnan(data2[i]))) {
      (*data1_nomissing)[*m_nomissing] = data1[i];
      (*data2_nomissing)[*m_nomissing] = data2[i];
      (*m_nomissing) ++;
    }

  }
}

//
//  remove missing data from a vector, and output a mask of the removed entries
//
void removeMissingValuesFloat(float* data1, int m, float** data1_nomissing, int* m_nomissing, char** mask) 
{
  int i;
  *data1_nomissing = (float*)malloc(m * sizeof(float));
  *mask            = (char*) malloc(m * sizeof(char));

  *m_nomissing = 0;
  for (i=0; i<m; i++) {
    //printf("%f\n", data1[i]);
    if (!isinf(data1[i])) {
      (*data1_nomissing)[*m_nomissing] = data1[i];
      (*m_nomissing) ++;
      (*mask)[i]     = 1;
    } else {
      (*mask)[i]     = 0;
    }

  }

}

//
//  get a vector whose original entries have been masked 
//
void  getMaskedVector(int* data1, int m, char* mask, int** data1_nomissing) 
{
  int i;
  
  int m_nomissing = 0;
  
  *data1_nomissing = (int*)malloc(m * sizeof(int));
  for (i=0; i<m; i++) {
    if (mask[i] == 1) {
      (*data1_nomissing)[m_nomissing] = data1[i];
      m_nomissing ++;
    }
  }
}



void  removeMissingValuesIntFloat(int* data1, float* data2, 
				    int m, 
				    int** data1_nomissing, 
				    float** data2_nomissing,
				    int* m_nomissing) 
{
  int i;

  *data1_nomissing = (int*)malloc(m * sizeof(int));
  *data2_nomissing = (float*)malloc(m * sizeof(float));
  *m_nomissing = 0;

  for (i=0; i<m; i++) {
    if ((data1[i] != INT_MAX) && (!isinf(data2[i]))) {
      (*data1_nomissing)[*m_nomissing] = data1[i];
      (*data2_nomissing)[*m_nomissing] = data2[i];
      (*m_nomissing) ++;
    }

  }
}



void  removeMissingValuesIntInt(int* data1, int* data2, 
				    int m, 
				    int** data1_nomissing, 
				    int** data2_nomissing,
				    int* m_nomissing) 
{
  int i;

  *data1_nomissing = (int*)malloc(m * sizeof(int));
  *data2_nomissing = (int*)malloc(m * sizeof(int));
  
  if ((*data1_nomissing == 0) || (*data1_nomissing == 0)) {
    printf("removeMissingValuesIntInt : not enough memory\n");
    exit(0);
  }
  
  *m_nomissing = 0;

  for (i=0; i<m; i++) {
    if ((data1[i] != INT_MAX) && (data2[i] != INT_MAX)) {
      (*data1_nomissing)[*m_nomissing] = data1[i];
      (*data2_nomissing)[*m_nomissing] = data2[i];
      (*m_nomissing) ++;
    } //else {
      //printf("GUESS WHAT ..(%d)\n", i);
    //}

  }
}


float CalculateMIFromUniquePreQuantizedVectors(int* quantized_vector1, int D1, int* quantized_vector2, int D2, int m, int Repeats, int remove_missing_values, int correct_sample_size) 
{
  int** qvs1;
  int** qvs2;
  int*  vqd1;
  int*  vqd2;
  float mi;
  int   gogo1, gogo2;
    
  qvs1    = (int**)malloc(1 * sizeof(int*));
  qvs1[0] = quantized_vector1;

  qvs2    = (int**)malloc(1 * sizeof(int*));
  qvs2[0] = quantized_vector2;

  vqd1    = (int*)malloc(1 * sizeof(int));
  vqd1[0] = D1;

  vqd2    = (int*)malloc(1 * sizeof(int));
  vqd2[0] = D2;

  mi      = CalculateMIFromPreQuantizedVectors(qvs1, vqd1, 1, qvs2, vqd2, 1, m, &gogo1, &gogo2, Repeats, remove_missing_values, correct_sample_size);

  free(qvs1); free(qvs2); free(vqd1); free(vqd2);
  
  return mi;


}



//
//  use a pre-quantization
//
float CalculateMIFromPreQuantizedVectors(int** quantized_vectors1, int* vqd1, int nqd1, int** quantized_vectors2, int* vqd2, int nqd2, int m, int* OD1, int* OD2, int Repeats, int remove_missing_values, int correct_sample_size) 
{
  int*    qv1;
  int*    qv2;
  //int     Repeats = REPEATS;
  int     D1;
  int     D2;
  int     FullSample;
  float   SampleSizesPerc[4] = {0.7, 0.7875, 0.9, 1.0};
  int     Trials [4]         = { 21,     16,  12,   1};
  int     D;
  int     i;
  float   mi;
  float   std = -1;
  int     maxnbbins = 5;
  int     nbbins = 0;
  float*  mymuts;
  float*  mystds;
  int     N = 4;

  float   max_mi;
  
  int     ID1, ID2;
  int     my_OD1, my_OD2;

  int*    quantized_vectors1_nomissing;
  int*    quantized_vectors2_nomissing;
  int     m_nomissing;

  //int correct_sample_size = 0;
  int**   tmpCounts;

  my_OD1 = -1;
  my_OD2 = -1;
  
  // printf("nqd1*nqd2=%d\n", nqd1 * nqd2);

  mymuts   = (float*)calloc(nqd1 * nqd2, sizeof(float));
  mystds   = (float*)calloc(nqd1 * nqd2, sizeof(float));


  for (ID1=0; ID1<nqd1; ID1++) {
        
    for (ID2=0; ID2<nqd2; ID2++) {

      if (remove_missing_values == 1) {
	//printf("ID1 = %d\n", ID1);
	//showIntVector(quantized_vectors1[ID1], m);

	removeMissingValuesIntInt(quantized_vectors1[ID1], quantized_vectors2[ID2], m, &quantized_vectors1_nomissing, &quantized_vectors2_nomissing, &m_nomissing);
	//showIntVector(quantized_vectors1_nomissing, m_nomissing);
      } else {
	quantized_vectors1_nomissing = quantized_vectors1[ID1];
	quantized_vectors2_nomissing = quantized_vectors2[ID2];
	m_nomissing                  = m;
      }
	  
      if (correct_sample_size == 1) {
	mi = FindInfo(quantized_vectors1_nomissing, quantized_vectors2_nomissing, m_nomissing, vqd1[ID1], vqd2[ID2], SampleSizesPerc, N, Trials, Repeats, &std); 


      } else {
	tmpCounts   = ConstructCountMatrix (quantized_vectors1_nomissing, quantized_vectors2_nomissing, m_nomissing, vqd1[ID1], vqd2[ID2]);
	
	//showCountMatrix(tmpCounts, vqd1[ID1], vqd2[ID2], NULL, NULL, NULL); //getchar();
	

	mi = FindInfoLocalMI(tmpCounts, vqd1[ID1], vqd2[ID2]);
	for (i=0; i<vqd1[ID1]; i++) {
	  free(tmpCounts[i]);
	}
	free(tmpCounts);
      }
      
      if (nbbins == 0) {
	max_mi = mi;
	my_OD1   = vqd1[ID1];
	my_OD2   = vqd2[ID2]; 
      } else {
	for (i=0; i<nbbins; i++) {
	  if (mi - 0.5*std > mymuts[i] + 0.5 * mystds[i]) {
	    max_mi = mi;
	    my_OD1   = vqd1[ID1];
	    my_OD2   = vqd2[ID2]; 
	  }
	}
      }
      
      // printf("nbbins=%d\n", nbbins);
      
      mystds[ nbbins ] = std;
      mymuts[ nbbins ] = mi;
      
      nbbins ++;
    
    
      if (remove_missing_values == 1) {
	free(quantized_vectors1_nomissing);
	free(quantized_vectors2_nomissing);
      }

      // printf("mi(D=%d)=%3.2f, std=%f\n", D, mi, std);
      
      //getchar();
    }

  }

  *OD1 = my_OD1;
  *OD2 = my_OD2;

  free(mymuts);
  free(mystds);

  return max_mi;
}





float QuantizeAndCalculateMI(float* data1, float* data2, int m, int* quantized_data1, int QD1, int* quantized_data2, int QD2) 
{
  int*    qv1;
  int*    qv2;
  int     Repeats = REPEATS;
  int     D1;
  int     D2;
  int     FullSample;
  float   SampleSizesPerc[4] = {0.7, 0.7875, 0.9, 1.0};
  int     Trials [4]         = { 21,     16,  12,   1};
  int     D;
  int     i;
  float   mi;
  float   std;
  int     maxnbbins = 5;
  int     nbbins = 0;
  float*  muts;
  float*  stds;
  int     N = 4;

  float   max_mi;
  
  int     bins1[4]      = {2, 3, 4, 5};
  int     maxnbbins1 = 4;
  int     bins2[4]      = {2, 3, 4, 5};
  int     maxnbbins2 = 4;
  int     ID1, ID2;

  float* bins;
  
  if (quantized_data1 != 0) {
    bins1[0]   = QD1;
    maxnbbins1 = 1;
  }
  
  if (quantized_data2 != 0) {
    bins2[0]   = QD2;
    maxnbbins2 = 1;
  }

  
  muts = (float*)calloc(maxnbbins1 * maxnbbins2, sizeof(float));
  stds = (float*)calloc(maxnbbins1 * maxnbbins2, sizeof(float));


  for (ID1=0; ID1<maxnbbins1; ID1++) {
    
    //printf("Quantize for D = %d\n", D);
    
    if (quantized_data1 == 0) {
      D1  = bins1[ID1];
      qv1 = Quantize(data1, m, D1, NULL); 
    } else { 
      qv1 = quantized_data1;  
      D1  = QD1;
    }

    for (ID2=0; ID2<maxnbbins2; ID2++) {

      if (quantized_data2 == 0) {
	D2  = bins2[ID2];
	qv2 = Quantize(data2, m, D2, NULL);
      } else {
	qv2 = quantized_data2;
	D2  = QD2;
      }

      //for (i=0; i<m; i++) {
      //  printf("qv1[%d]=%d\n", i, qv1[i]);
      //  printf("qv2[%d]=%d\n", i, qv2[i]);
      //}
            
      mi = FindInfo(qv1, qv2, m, D1, D2, SampleSizesPerc, N, Trials, Repeats, &std); 
      
      //printf("mi for %d and %d bins = %f\n", D1, D2, mi);
      
      if (nbbins == 0) {
	max_mi = mi;
      } else {
	for (i=0; i<nbbins; i++) {
	  if (mi - 0.5*std > muts[i] + 0.5 * stds[i])
	    max_mi = mi;
	}
      }
      
      stds[ nbbins ] = std;
      muts[ nbbins ] = mi;
      
      nbbins ++;
      
      // printf("mi(D=%d)=%3.2f, std=%f\n", D, mi, std);
      
      //getchar();
    }

  }
  
  return max_mi;
}


int CmpDbl(const void* _a, const void* _b) {

  const float a = *((const float*) _a);
  const float b = *((const float*) _b);

  if (a > b) 
    return 1;
  else if (a == b) {
    return 0;
  } else {
    return -1;
  }
}





int* logicCombine(int* data1,  int* data2, int m, int D1, int D2, char* type) 
{
 
  int* newdata;
  int  i, j;

  newdata = (int*)malloc(m * sizeof(int));
  if (newdata == 0) {
    printf("combineQuantizedVectors : not enough memory\n");
    exit(0);
  }
  for (i=0; i<m; i++) {
    
    
    if ((data1[ i ] == INT_MAX) || (data2[ i ] == INT_MAX))
      newdata [ i ] = INT_MAX;
    else {

      if (strcmp(type, "AND") == 0) {
	newdata [ i ] = data1[ i ]  & data2[ i ];
      } else 
	if (strcmp(type, "XOR") == 0) {
	  newdata [ i ] = data1[ i ]  ^ data2[ i ];
	} else 
	  if (strcmp(type, "OR") == 0) {
	    newdata [ i ] = data1[ i ] | data2[ i ];
	  } else 
	    if (strcmp(type, "NOT") == 0) {
	      newdata [ i ] = data1[ i ] & !data2[ i ]; // | (!data1[ i ] & data2[ i ]));
	    }
    }
    
  }

  return newdata;
  
}


float entropy(int* data, int m, int D) 
{
  
  int bins[D];
  int i;
  float h = 0.0;
  float    e = 2.204e-16;
  
  for (i=0; i<D; i++) {
    bins[ i ] = 0; 
  }

  for (i=0; i<m; i++) {
    bins[ data[i] ] ++; 
  }

  for (i=0; i<D; i++) {
    //printf("h+= %3.2f * %3.2f(log) / log(2.0)\n", (float)bins[i]/(float)m, log((float)bins[i]/(float)m));
    h += ( (float)bins[i] + e )/(float)m * log( ((float)bins[i] + e ) / (float)m) / log(2.0); 
  }

  return - h;

}




int CmpFI(const void* _a, const void* _b)
{
  const FloatAndIndex* a = (const FloatAndIndex*) _a;
  const FloatAndIndex* b = (const FloatAndIndex*) _b;
  
  if (a->v < b->v) 
    return -1; 
  else if(a->v == b->v) 
    return  0;
  else         
    return  1;
}

