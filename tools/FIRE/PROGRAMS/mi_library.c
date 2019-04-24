

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "information.h"
#include "statistics.h"
#include "dataio.h"
#include "mi_library.h"





void encode_kmer(char* kmer, int kmersize, int add5, int add3, char** d_ch, char** ch, int nbchars, int** encoded_kmer, char** encoded_re)
{
  int k, l;
  
  for (k=0; k<add5; k++) {
    (*encoded_kmer)[k] = 0;
    strcat(*encoded_re, ".");
  }
  
  for (k=0; k<kmersize; k++) {
    for (l=0; l<nbchars; l++) {       
      if (kmer[k] == d_ch[l][0]) {
	(*encoded_kmer)[add5+k] = l;
	strcat(*encoded_re, ch[ l ]);
	break;
      }     
    }
  }
  
  for (k=0; k<add3; k++) {
    (*encoded_kmer)[add5+kmersize+k] = 0;
    strcat(*encoded_re, ".");
  }

}


void bubbleSortFI(FloatAndIndex* fi, int n)
{
  int i, j;

  // for storing temporary values
  int   itmp;
  float ftmp;


  for (i=0; i<n-1; i++) {
    for (j=0; j<n-1-i; j++)
      if (fi[j+1].v < fi[j].v) {  /* compare the two neighbors */

        ftmp = fi[j].v;         /* swap a[j] and a[j+1]      */
        fi[j].v = fi[j+1].v;
        fi[j+1].v = ftmp;

        itmp = fi[j].i;         /* swap a[j] and a[j+1]      */
        fi[j].i = fi[j+1].i;
        fi[j+1].i = itmp;

      }
  }

}

void add_small_values_to_identical_floats(float* E, int n) {
  
  //float* Ec;
  FloatAndIndex* fi;
  float  eps = 2.204e-16;
  float  d;
  float  d_prevdiff;
  float  d_nextdiff;
  int    i;
  float* nextdiff;

  if (n <= 1)
    return;
  
  fi      = (FloatAndIndex*)malloc( n * sizeof( FloatAndIndex ) );
  for (i=0; i<n; i++) {
    fi[i].i = i;
    fi[i].v = E[i];
  }
  //myqsort((void*)fi, n, sizeof(FloatAndIndex), CmpFI);
  
  bubbleSortFI(fi, n);
  
  // create a vector of next differences
  nextdiff      = (float*)malloc(n * sizeof(float));
  d_nextdiff    = 0.1;
  nextdiff[n-1] = d_nextdiff;
  for (i=n-2; i>=0; i--) {
    d = fabs(fi[i+1].v - fi[i].v);
    if (d > eps) {
      d_nextdiff = d; // update nextdiff
    }    
    nextdiff[i] = d_nextdiff;
  }

  E[ fi[0].i ] = fi[0].v;
  d_prevdiff = 0.1;
  //printf("%f\t%f\t%f\t%f\n", d_prevdiff, nextdiff[0], fi[0].v, Ec[0]);
	
  for (i=1; i<n; i++) {
    d = fabs(fi[i-1].v - fi[i].v);
    
    //printf("d=%f\n", d);
    
    if (d < eps) {
      E[ fi[i].i  ] = fi[i].v + min(0.1, min( nextdiff[i],  d_prevdiff )) * default_rand();
      //printf("%d\t%f\t%f\t%f\t%f\n", fi[i].i, d_prevdiff, nextdiff[i], fi[i].v, E[ fi[i].i  ]);
    } else {
      E[ fi[i].i  ] = fi[i].v;
      d_prevdiff = d; 
      //printf("%d\t%f\t%f\t%f\t%f\n", fi[i].i, d_prevdiff, nextdiff[i], fi[i].v, E[ fi[i].i  ]);

    }

    
  }
  
  free(nextdiff);
  free(fi);
  
}


float calc_gc_content(char* seq, int l)
{
  int  gc = 0; 
  int  l_p = 0;
  int  i;
  
  for (i=0; i<l; i++) {
    if ((seq[i] == 'C') || (seq[i] == 'G')) 
      gc ++;
    if ((seq[i] == 'C') || (seq[i] == 'G') || (seq[i] == 'T') || (seq[i] == 'A') ) 
      l_p ++;
  }

  if (l_p > 0) 
    return (float)gc / l_p;
  else
    return 0.0;
}


float calc_CpG_content(char* seq, int l)
{
  int  gc = 0; 
  int  l_p = 0;
  int  i;
  
  for (i=0; i<l-1; i++) {
    if ((seq[i] == 'C') && (seq[i+1] == 'G')) 
      gc ++;
    if (((seq[i] == 'C') || (seq[i] == 'G') || (seq[i] == 'T') || (seq[i] == 'A') ) && 
	((seq[i+1] == 'C') || (seq[i+1] == 'G') || (seq[i+1] == 'T') || (seq[i+1] == 'A') ))
      l_p ++;
  }

  if (l_p > 0) 
    return (float)gc / l_p;
  else
    return 0.0;
}


void get_Params(int argc, char** argv, Params* p) 
{
  
  p->expfile         = get_parameter(argc, argv, "-expfile");

  if (exist_parameter(argc, argv, "-kmerfile")) 
    p->kmerfile        = get_parameter(argc, argv, "-kmerfile");
  else if (exist_parameter(argc, argv, "-motiffile")) 
    p->kmerfile        = get_parameter(argc, argv, "-motiffile");
  
  p->fastafile       = get_parameter(argc, argv, "-fastafile");
  p->kmersize        = atoi(get_parameter(argc, argv, "-k"));
  
  
  if (exist_parameter(argc, argv, "-quantized")) 
    p->quantized = atoi(get_parameter(argc, argv, "-quantized"));
   
  if (exist_parameter(argc, argv, "-shuffle")) 
    p->shuffle = atoi(get_parameter(argc, argv, "-shuffle"));

  if (exist_parameter(argc, argv, "-shuffle_rank")) 
    p->shuffle_rank = atoi(get_parameter(argc, argv, "-shuffle_rank"));
   
  if (exist_parameter(argc, argv, "-rna")) 
    p->singlestrand    = atoi(get_parameter(argc, argv, "-rna"));   

  if (exist_parameter(argc, argv, "-verbose")) 
    p->verbose         = atoi(get_parameter(argc, argv, "-verbose"));   
 
  if (exist_parameter(argc, argv, "-mbins")) 
    p->mbins         = atoi(get_parameter(argc, argv, "-mbins"));   

  if (exist_parameter(argc, argv, "-mbins_dist")) 
    p->mbins_dist         = atoi(get_parameter(argc, argv, "-mbins_dist"));   
  
  if (exist_parameter(argc, argv, "-outfile")) 
    p->outfile         = get_parameter(argc, argv, "-outfile");   
  
  if (exist_parameter(argc, argv, "-report")) 
    p->report         = atoi(get_parameter(argc, argv, "-report"));
  

}




//
//  simple max rank test
//
int max_rank_test_cond(double score, int* M_q, int mbins, int* E_q, int ebins, int* A_q, int abins, int nborfs, int shuffle, double* max_cmi_shu, int* val, double* z) 
{
  int       j;
  int*      E_q_shu;
  double    mymi;
  double*   a_mi_shu;
  double    mi_shu_avg;
  double    mi_shu_std;
  double    c;
  // double    cmi;

  a_mi_shu = (double*)malloc(shuffle * sizeof(double));
  
  for (j=0; j<shuffle; j++) {
    E_q_shu = shuffleInt(E_q, nborfs); 
    mymi    = CalculateCondMIbasic(M_q, E_q_shu, A_q, nborfs, mbins, ebins, abins); 
    a_mi_shu[j] = mymi;       
    free(E_q_shu);
  }
  
  qsort((void*)a_mi_shu, shuffle, sizeof(double), CmpDblRegular);

  c  = a_mi_shu[ shuffle - 1  ];
  
  *max_cmi_shu = c;
  
  //if (do_val == 1) {
    j = shuffle - 1;
    while ((j >= 0) && (score < a_mi_shu[j])) {   // go while the shuffled score is higher than the real one
      j --;
    }
    *val = shuffle - j - 1;
  //}

//if (do_z == 1) {
    mi_shu_avg = average_dbl(a_mi_shu, shuffle);
    mi_shu_std = stddev_dbl (a_mi_shu, shuffle);

    *z = ( score - mi_shu_avg ) / mi_shu_std;
//}

  free(a_mi_shu);

  if (score > c) 
    return 1;
  else
    return 0;
}




//
//  simple max rank test
//
int max_rank_test(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int do_val, int* val, int do_z, double* z) 
{
  int      j;
  int*     E_q_shu;
  double    mymi;
  double*   a_mi_shu;
  double    mi_shu_avg;
  double    mi_shu_std;
  double    c;

  a_mi_shu = (double*)malloc(shuffle * sizeof(double));
  
  for (j=0; j<shuffle; j++) {
    E_q_shu = shuffleInt(E_q, nborfs); 
    mymi = CalculateMIbasic(M_q, E_q_shu, nborfs, mbins, ebins);
    a_mi_shu[j] = mymi;       
    free(E_q_shu);
  }
  
  qsort((void*)a_mi_shu, shuffle, sizeof(double), CmpDblRegular);

  c  = a_mi_shu[ shuffle - 1 ];
  
  if (do_val == 1) {

    if ((c < 1e-10) && (score < 1e-10))  // shortcut
      *val = shuffle; 
    else if (score <= a_mi_shu[0]) {  // other shortcut
      *val = shuffle;
    } else {
      j = shuffle - 1;
      while ((j >= 0) && (score <= a_mi_shu[j])) {   // go while the shuffled score is higher than the real one
	j --;
      }
      *val = shuffle - j - 1;
    }

  }

  if (do_z == 1) {
    mi_shu_avg = average_dbl(a_mi_shu, shuffle);
    mi_shu_std = stddev_dbl (a_mi_shu, shuffle);

    *z = ( score - mi_shu_avg ) / mi_shu_std;

    if (isnan(*z)) {
      
      //printf("c=%e,s=%e, %s,%s\n", c, score, (c>score?"Y":"N"),(c<score?"Y":"N"));
      
      *z = 0;
      //printf("c=%5.4f, std=%5.4f\n", c, mi_shu_std);
    }
  }

  free(a_mi_shu);

  //printf("%f\t%f\n", score, c);

  if (score > c) 
    return 1;
  else
    return 0;
}


int jacknife_max_rank_test(int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int jn, int jn_f, int jn_t, int do_val, int* val) 
{

  int*  M_q_cross;
  int*  E_q_cross;
  int   pass = 0;
  int   l;
  int   nborfs_retained;
  float r;
  int   j;
  double mymi;

  M_q_cross = (int*)malloc(nborfs * sizeof(int));
  E_q_cross = (int*)malloc(nborfs * sizeof(int));

  for (l=0; l<jn; l++) {
    //printf("l=%d\n", l);
    // keep random 90% of the data
    nborfs_retained = 0;
    for (j=0; j<nborfs; j++) {
      //r = rand()/(RAND_MAX+1.0); 
      r       = default_rand();
      
      if (r > 1.0/jn_f) {
	M_q_cross[ nborfs_retained ] = M_q[j];
	E_q_cross[ nborfs_retained ] = E_q[j];
	nborfs_retained ++;
      }
    }
    
    //    printf("nborfs_retained = %d\n", nborfs_retained);

    // eval mi
    mymi      = CalculateMIbasic(M_q_cross, E_q_cross, nborfs_retained, mbins, ebins);
    pass     += max_rank_test(mymi, M_q_cross, mbins, E_q_cross, ebins, nborfs_retained, shuffle, 0, 0, 0, 0); 

  }
  
  //printf("pass=%d\n", pass);
  
  free(M_q_cross);
  free(E_q_cross);

  if (do_val == 1) {
    *val = pass;
  }

  if ( pass >= jn_t )
    return 1;
  else 
    return 0;

}
  

void get_rank_and_zscore(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int* rank, double* zscore) 
{
  int      j;
  int*     M_q_shu;
  double    mymi;
  double*   a_mi_shu;
  double    mi_shu_avg;
  double    mi_shu_std;
  double    z;
  int cnt_higher = 0;

  a_mi_shu = (double*)malloc(shuffle * sizeof(double));
  
  for (j=0; j<shuffle; j++) {
    M_q_shu = shuffleInt(M_q, nborfs); 
    mymi = CalculateMIbasic(M_q_shu, E_q, nborfs, mbins, ebins);
    a_mi_shu[j] = mymi;       
    if (mymi >= score) {
      cnt_higher++;
    }
    free(M_q_shu);
  }
       
  mi_shu_avg = average_dbl(a_mi_shu, shuffle);
  mi_shu_std = stddev_dbl (a_mi_shu, shuffle);
  
  z = ( score - mi_shu_avg ) / mi_shu_std;
  
  *zscore = z;
  *rank   = cnt_higher;

  free(a_mi_shu);

  return;
  
}
 


double get_zscore_and_rank_value(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int shuffle_rank, double* c) 
{
  int      j;
  int*     M_q_shu;
  double    mymi;
  double*   a_mi_shu;
  double    mi_shu_avg;
  double    mi_shu_std;
  double    z;
  
  a_mi_shu = (double*)malloc(shuffle * sizeof(double));
  
  for (j=0; j<shuffle; j++) {
    M_q_shu = shuffleInt(M_q, nborfs); 
    mymi = CalculateMIbasic(M_q_shu, E_q, nborfs, mbins, ebins);
    a_mi_shu[j] = mymi;       
    free(M_q_shu);
  }
       
  mi_shu_avg = average_dbl(a_mi_shu, shuffle);
  mi_shu_std = stddev_dbl (a_mi_shu, shuffle);
  
  z = ( score - mi_shu_avg ) / mi_shu_std;
  
  qsort((void*)a_mi_shu, shuffle, sizeof(double), CmpDblRegular);

  *c  = a_mi_shu[ shuffle_rank - 1  ];

  free(a_mi_shu);

  return z;
  
}


double get_zscore_and_do_interval(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int nbkmers, double* c) 
{
  int      j;
  int*     M_q_shu;
  double    mymi;
  double*   a_mi_shu;
  double    mi_shu_avg;
  double    mi_shu_std;
  double    z;
  int      ci;

  a_mi_shu = (double*)malloc(shuffle * sizeof(double));
  
  for (j=0; j<shuffle; j++) {
    M_q_shu = shuffleInt(M_q, nborfs); 
    mymi = CalculateMIbasic(M_q_shu, E_q, nborfs, mbins, ebins);
    a_mi_shu[j] = mymi;       
    free(M_q_shu);
  }
       
  mi_shu_avg = average_dbl(a_mi_shu, shuffle);
  mi_shu_std = stddev_dbl (a_mi_shu, shuffle);
  
  z = ( score - mi_shu_avg ) / mi_shu_std;

  
  qsort((void*)a_mi_shu, shuffle, sizeof(double), CmpDblRegular);
  
  ci = (int)(0.5 + (1 - 0.05 / nbkmers) * shuffle ) - 1;
  *c  = a_mi_shu[ ci ];

  free(a_mi_shu);

  return z;
  
}


double get_zscore(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle) 
{
  int      j;
  int*     M_q_shu;
  double    mymi;
  double*   a_mi_shu;
  double    mi_shu_avg;
  double    mi_shu_std;
  double    z;

  a_mi_shu = (double*)malloc(shuffle * sizeof(double));
  
  for (j=0; j<shuffle; j++) {
    M_q_shu = shuffleInt(M_q, nborfs); 
    mymi = CalculateMIbasic(M_q_shu, E_q, nborfs, mbins, ebins);
    a_mi_shu[j] = mymi;     
    //printf("%f\n", mymi);
    free(M_q_shu);
  }
       
  mi_shu_avg = average_dbl(a_mi_shu, shuffle);
  mi_shu_std = stddev_dbl (a_mi_shu, shuffle);
  
  //printf("mi_shu_avg=%f, mi_shu_std=%f\n", mi_shu_avg, mi_shu_std);

  z = ( score - mi_shu_avg ) / mi_shu_std;

  free(a_mi_shu);

  return z;
  
}


void add_to_report(FILE* fpr, char* kmer, int kmersize, int gap, int* M_q, int mbins, int* E_q, int ebins, float* E_q_bins, int nborfs)
{

  char*    stmp;
  int      j;
  int**    tmpCounts;
  MI_Contribs* mico;
  int      im; //, ie;
  double    mi;
  int*     mdensity;
  int*     esize;

  

  stmp = (char*)calloc( 100, sizeof(char));
  
  if (gap > 0) {
    strncpy(stmp, kmer, kmersize/2);  
    stmp[kmersize/2] = '\0';
    for (j=0; j<gap; j++) 
      strcat(stmp, "N");
    strcat(stmp, kmer + kmersize/2);
  } else {
    strcat(stmp, kmer);
    kmersize = strlen(kmer);
    stmp[kmersize] = '\0';
  }
  
  tmpCounts = ConstructCountMatrix (M_q, E_q, nborfs, mbins, ebins);
  mi        = FindInfoLocalMI(tmpCounts, mbins, ebins);

  InformationGain_D2(tmpCounts, mbins, ebins, &mico, 0); 

  mdensity = (int*)calloc( ebins, sizeof(int));
  esize    = (int*)calloc( ebins, sizeof(int));

  for (j=0; j<nborfs; j++) {
    if (M_q[j] == 1)
      mdensity[ E_q[j] ] ++;
    esize[ E_q[j] ] ++;
  }

  for (im=0; im<ebins; im++) {
    fprintf(fpr, "%s\t%d\t%d\t%d", stmp, mico[im].i, (int)(0.5 + 100 * mico[im].co / mi ), mico[im].j);
    
    if (E_q_bins != 0)
      fprintf(fpr, "\t%4.3f\t%4.3f", E_q_bins[im], E_q_bins[im+1]); 
    
    fprintf(fpr, "\t%4.3f\t%d", (mdensity[im]==0?0:mdensity[im] / (float)(esize[im])), esize[im]);
    
    fprintf(fpr, "\n");
  }

  //  fprintf(fpr, "\n");
  free(mico);
  for (im=0; im<mbins; im++)
    free(tmpCounts[im]); 
  free(tmpCounts);
  free(stmp);
  
  free(mdensity);
  free(esize);

}




int CmpShort(short* a, short* b)
{
  if (*a<*b)
    return -1;
  else if (*a == *b) 
    return 0;
  else 
    return 1;
}

void quantize_M_counts( short* M, int nborfs, int* mbins, int** M_q) 
{
  int j;
  int max_count = 0;
  int* counts;
  int  cnt;
  extern int CmpShort();
  short* MC;

  //
  //  sort the count profile
  //
  MC = (short*)malloc( nborfs * sizeof( short ));
  memcpy(MC, M, nborfs * sizeof( short ));

  qsort(MC, nborfs, sizeof(short), CmpShort);
  max_count = MC[nborfs-1];
 
  //printf("maxcount = %d\n", max_count);

  counts = (int*)calloc( max_count+1,  sizeof(int));
  for (j=0; j<=max_count; j++) {
    counts[ j ] = -1;
  }

  cnt = 0;
  for (j=0; j<nborfs; j++) {
    if ( counts[ MC[ j ] ] == -1 ) {
      counts[ MC[ j ] ] = cnt;
      //printf("CNT=%d, b=%d\n", MC[j], cnt);
      cnt ++;
      
    }
  }

  //
  // finds the maximum number of counts
  //
  *M_q = (int*)malloc(nborfs * sizeof(int));
  for (j=0; j<nborfs; j++) {
    (*M_q)[j] = counts[ M[j] ];
    //if ((*M_q)[j] == 137) {
    //  printf("j=%d for %d\n", j, M[j]);
    //}
    //printf("for %d, I will have %d\n", M[j], counts[ M[j] ]);
  }

  free(counts);
  free(MC);

  *mbins = cnt;

}
  

//
//
//
void quantize_M_zero_eqpop( float* M, int nborfs, int mbins, int** M_q) 
{
  int j;
  int max_count = 0;
  int* counts;
  //int  cnt;
  extern int CmpShort();
  //short* MC;
  int    s;
  int    sum;
  int    tsum = 0;
  float  ratio;
  int    k;
  float  ratprev;
  int    sumprev;

  int*   st;
  int*   en;

  *M_q = (int*)malloc(nborfs * sizeof(int));
  
  if (mbins == 2) {
    // simple case
    for (j=0; j<nborfs; j++) {

      if (M[j]>0.5) {
	(*M_q)[j] = 1; 
      } else {
	(*M_q)[j] = 0;
      }
    }

  } else {
    // more complicated case
    st = (int*)malloc( mbins * sizeof( int ));
    en = (int*)malloc( mbins * sizeof( int ));
    
    st[0] = 0; en[0] = 1;
    
    
    max_count = -1; 
    for (j=0; j<nborfs; j++) {
      if (M[j] > max_count) {
	max_count = (float)M[j];
      }
    }
    
    
    counts = (int*)calloc( max_count+1,  sizeof(int));
    for (j=0; j<nborfs; j++) {
      counts[ (int)M[j] ] ++;
    }
    
    if (counts[0] > 0) 
      s = 1;
    else 
      s = 0;
    
    tsum = 0;
    
    for (j=s; j<=max_count; j++) {
      tsum  += counts[ j ];
    }
    
    k = 1;
    j = s;
    while (k <= (mbins-1)) {
      
      sum   = 0;
      ratio = -1;
      j     = s;
      ratprev = -1;
      
      while ((ratio < 1.0/(mbins-1)) && (j<=max_count)) {
	
	sumprev = sum;
	ratprev = ratio;

	sum   += counts[ j ];
	ratio  = sum / (float)tsum;
	
	j ++;
      } 
      
      //
      // in theory here, should look at one count before (in some cases)
      //
      
      st[k] = s; en[k] = j;
      s = j;
      
      // try one step more
      k ++;
    }
    
    
    for (j=0; j<nborfs; j++) {
      for (k=0; k<mbins; k++) {
	if (((int)M[j] >= st[k]) && ((int)M[j] < en[k])) {
	  (*M_q)[j] = k;
	  break;
	}
      }
    }
    
    free(st);
    free(en);
  }


}




void quantize_E(float* E, int nborfs, int quantized, int* ebins, int** E_q, float** E_q_bins) 
{
  int k;
  
  if (quantized == 0) {

     for (k=0; k<nborfs; k++) {
       if (isnan(E[k]))
	 E[k] = 0.0;
     }
     
     //
     //  quantize
     //
     *E_q = Quantize(E, nborfs, *ebins, E_q_bins);

   } else {

    *ebins = -1;

    // finds the maximum count
    for (k=0; k<nborfs; k++) {
      if ((int)(E[k]) > *ebins) {
	*ebins = (int)(E[k]);
      }
    }
    *ebins = *ebins + 1;
    
    *E_q = (int*)malloc(nborfs * sizeof(int));
    for (k=0; k<nborfs; k++) {
      (*E_q)[k] = (int)(E[k]);
    }
    
    E_q_bins  = 0;

  }
  
  
}




void quantize_M( float* M, int nborfs, int mbinary, int mnkmax, int* mbins, int** M_q) 
{
  int j;
  int Mi;

  if (mbinary == 1) {
    *M_q = (int*)malloc(nborfs * sizeof(int));
    for (j=0; j<nborfs; j++) {
      if (M[j]>0.5) {
	(*M_q)[j] = 1; 
      } else {
	(*M_q)[j] = 0;
      }
      //printf("%d\n", (*M_q)[j]);

    }
    *mbins = 2;
    
  } else if (mnkmax > 0) {

    *M_q = (int*)malloc(nborfs * sizeof(int));

    for (j=0; j<nborfs; j++) {
      
      Mi = (int)(M[j]);
      
      if (Mi > mnkmax) {
	(*M_q)[j] = mnkmax; 
      } else {
	(*M_q)[j] = Mi;
      } 

      //printf("%f\t%d\n", (M[j]), (*M_q)[j]);

    }
    *mbins = mnkmax+1;
    
  } else {

    *M_q = Quantize(M, nborfs, *mbins, NULL);
    die("this type of quantization is not allowed now\n");
  }

}
  



void readKmers (char* kmer, char* kmerfile, int max_nb_kmers, int inc_nb_kmers, char*** kmers, float** oldmis, int* nbkmers, int* kmersize) 
{
  char* buff;
  FILE* fp;
  char* s;
  int   len = 10000;
  //int   i;

  buff = (char*)malloc(1000 * sizeof(char));

  if (strlen(kmer) > 0) {
    *nbkmers = 1;
    *kmers  = (char**)malloc(*nbkmers * sizeof(char*));  
    *oldmis = (float*)malloc(max_nb_kmers * sizeof(float));
    (*kmers)[0] = strdup(kmer);
    
  } else {
    
    //
    // START READING THE LIST OF KMERS
    //
    *kmers  = (char**)malloc(max_nb_kmers * sizeof(char*));
    *oldmis = (float*)malloc(max_nb_kmers * sizeof(float));
    
    buff = (char*)malloc(10000 * sizeof(char));
    fp = fopen(kmerfile, "r");
    if (!fp) {
      printf("could not open kmer/seed file %s\n", kmerfile);
    }
    
    *nbkmers  = 0;
    *kmersize = 0;
    
    while (!feof(fp)) {

      fgets(buff, len, fp);

      if (feof(fp))
	break;
      
      s = mystrtok(buff, '\t');
      (*kmers)[*nbkmers] = (char*)calloc(100, sizeof(char));   
      strcat((*kmers)[*nbkmers], s);    
      free(s);
      
      s = mystrtok(0, '\t');
      (*oldmis)[ *nbkmers ] = atof( s );
      free(s);
      
      (*nbkmers)++;

      // allocate more memory if we need it
      if (*nbkmers == max_nb_kmers) {
	die("Olivier: allocate more space for the kmers\n");
	*kmers = allocate_more_kmers(*kmers, *nbkmers, max_nb_kmers, inc_nb_kmers);
	max_nb_kmers += inc_nb_kmers;
      }
    }
  
    

    // motif size becomes the length of the first k-mer in the list
    if (*nbkmers > 0) 
      *kmersize = strlen((*kmers)[0]);    
    else 
      *kmersize = -1;

    fclose(fp);     
    //
    // END READING THE LIST OF KMERS
    //
    
  }

}     



void readKmers_general_special_optim (char* kmerfile, int max_nb_kmers, int inc_nb_kmers, char*** motifs, char*** seeds, int* nbkmers, int* kmersize) 
{
  char* buff;
  FILE* fp;
  char* s;
  int   len = 10000;
  int   i;

  buff = (char*)malloc(1000 * sizeof(char));
    
  //
  // START READING THE LIST OF KMERS
  //
  *motifs  = (char**)malloc(max_nb_kmers * sizeof(char*));
  *seeds   = (char**)malloc(max_nb_kmers * sizeof(char*));
  
  
  buff = (char*)malloc(10000 * sizeof(char));
  fp = fopen(kmerfile, "r");
  if (!fp) {
    printf("could not open kmer/seed/motif file %s\n", kmerfile);
    exit(0);
  }
  
  *nbkmers  = 0;
  *kmersize = 0;
  
  while (!feof(fp)) {
    
    fgets(buff, len, fp);
    
    if (feof(fp))
      break;

    chomp(buff);

    s = mystrtok(buff, '\t');
    (*motifs)[*nbkmers] = (char*)calloc(100, sizeof(char));   
    strcat((*motifs)[*nbkmers], s);    
    free(s);
          
    for (i=1; i<=3; i++) {
      s = mystrtok(0, '\t'); free(s);
    }

    s = mystrtok(0, '\t');
    (*seeds)[*nbkmers] = (char*)calloc(100, sizeof(char));   
    strcat((*seeds)[*nbkmers], s);    
    free(s);
    
    (*nbkmers)++;

    // allocate more memory if we need it
    if (*nbkmers == max_nb_kmers) {
      die("Olivier: allocate more space for the seeds and motifs\n");
      
      *seeds  = allocate_more_kmers(*seeds,  *nbkmers, max_nb_kmers, inc_nb_kmers);
      *motifs = allocate_more_kmers(*motifs, *nbkmers, max_nb_kmers, inc_nb_kmers);
      
      max_nb_kmers += inc_nb_kmers;
    }
  }
  
  // motif size becomes the length of the first k-mer in the list
  if (*nbkmers > 0) {
    *kmersize = strlen((*seeds)[0]);    
  }

  fclose(fp);     
  //
  // END READING THE LIST OF KMERS
  //
  
}     



void readKmers_general (char* kmerfile, int max_nb_kmers, int inc_nb_kmers, char*** kmers, int* nbkmers, int* kmersize) 
{
  char* buff;
  FILE* fp;
  char* s;
  int   len = 10000;
  //int   i;

  buff = (char*)malloc(1000 * sizeof(char));
    
  //
  // START READING THE LIST OF KMERS
  //
  *kmers  = (char**)malloc(max_nb_kmers * sizeof(char*));
  
  
  buff = (char*)malloc(10000 * sizeof(char));
  fp = fopen(kmerfile, "r");
  if (!fp) {
    printf("could not open kmer/seed/motif file %s\n", kmerfile);
    exit(0);
  }
  
  *nbkmers  = 0;
  *kmersize = 0;
  
  while (!feof(fp)) {
    
    fgets(buff, len, fp);
    
    if (feof(fp))
      break;

    chomp(buff);

    s = mystrtok(buff, '\t');
    (*kmers)[*nbkmers] = (char*)calloc(200, sizeof(char));   
    strcat((*kmers)[*nbkmers], s);    
    free(s);
          
    (*nbkmers)++;
    
    // allocate more memory if we need it
    if (*nbkmers == max_nb_kmers) {
      die("Olivier: allocate more space for the kmers\n");
      *kmers = allocate_more_kmers(*kmers, *nbkmers, max_nb_kmers, inc_nb_kmers);
      max_nb_kmers += inc_nb_kmers;
    }
  }
  
  // motif size becomes the length of the first k-mer in the list
  *kmersize = strlen((*kmers)[0]);    
  fclose(fp);     
  //
  // END READING THE LIST OF KMERS
  //
  
}     








//
//  function to allocate more kmers
//
char** allocate_more_kmers(char** kmers, int n, int max_nb_kmers, int inc_nb_kmers) 
{

  char** more_kmers;
  int    i;

  more_kmers = (char**)malloc( (max_nb_kmers + inc_nb_kmers) * sizeof(kmers));
  for (i=0; i<n; i++) {
    more_kmers[i] = strdup(kmers[i]); 
  }
  
  return more_kmers;

}


