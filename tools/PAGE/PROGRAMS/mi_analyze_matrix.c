//
// calculate pairwise MI, z-scores
//


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <search.h>

#include "dataio.h"
#include "statistics.h"
#include "information.h"
#include "hashtable.h"
#include "mi_library.h"

#ifndef NAN
static const double NAN = (0.0 / 0.0);
#endif

typedef struct _MI {
  float zscore;
  float mi;
  int   rank;
  float pe;
} MI;

int main(int argc, char** argv) {
  
  char*   datafile = argv[1];
  float** data;
  char**  rownames;
  char**  colnames;
  int     n, m;
  
  int     i, j;

  MI**    mu;

  int     ebins1 = 10;
  int     ebins2 = 10;
  
  int*    E1_q;
  int*    E2_q;
  float*  E1_q_bins;
  float*  E2_q_bins;
  
  int     verbose = 1;
  int     k ;
  if (argc == 1) {
    printf("Usage : mi_analyze_matrix -datafile FILE -expnames FILE -quantized 0/1 [ -ebins INT ]\n");
    exit(0);
  }
  
  default_set_seed(12345);
  
  //						
  // read in profiles
  //
  readFloatTable(datafile, &m, &n, &data, &rownames, &colnames, 0, 1);

  if (verbose == 1) {
    printf("read matrix %d x %d\n", n, m);
  }

  mu = (MI**)malloc(m * sizeof(MI*));
  for (i=0; i<m; i++) {
    mu[i] = (MI*)malloc(m * sizeof(MI));
  }
  
  //
  // calc mi for all pairs
  //
  for (i=0; i<m-1; i++) {
    


    float* E1 = f_matrix_column(data, i, n);
    
   // getchar();
	

    for (j=i+1; j<m; j++) {
      
      float* E2 = f_matrix_column(data, j, n);

      float* E1_nom;
      float* E2_nom;
      int    n_nom;

      removeMissingValuesFloatFloat(E1, E2, n, &E1_nom, &E2_nom, &n_nom);


      quantize_E(E1_nom, n_nom, (colnames[i][0] == 'D'?1:0), &ebins1, &E1_q, &E1_q_bins);	
      quantize_E(E2_nom, n_nom, (colnames[j][0] == 'D'?1:0), &ebins2, &E2_q, &E2_q_bins);

      if (verbose == 1) {
	for (k=0; k<n_nom; k++) {
	  //printf("%d\t%d\n", E1_q[k], E2_q[k]);
	}
      }

      // correlate
      float mi = CalculateMIbasic(E1_q, E2_q, n_nom, ebins1, ebins2);      
      float pe = pearson(E1_nom, E2_nom, n_nom);

      // test
      int   shuffle = 10000;
      double zscore;
      int   rank = bounded_rank_test(mi, E1_q, ebins1, E2_q, ebins2, n_nom, shuffle, shuffle/10, 0, &zscore) ;
      
      mu[i][j].mi     = mi;
      mu[i][j].rank   = rank;
      mu[i][j].zscore = zscore;
      mu[i][j].pe     = pe;

      if (verbose == 1) {
	printf("%s\t%s\t%f\t%d\t%f\t%f\n", colnames[i], colnames[j], mi, rank, zscore, pe);
      }
      mu[i][j] = mu[j][i];

    }
    
  }
  return 0;

}
