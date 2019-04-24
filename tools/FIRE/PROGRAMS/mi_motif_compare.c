#define _GNU_SOURCE
#include <search.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>
#include "statistics.h"

#ifdef USEMPATROL
#include <mpatrol.h>
#endif

#include "dataio.h"
#include "statistics.h"
#include "information.h"
#include "mi_library.h"
#include "regexp.h"
#include "sequences.h"

typedef struct _MIARRAY {
 int    index;      // where is it found
 double score;
} MIArray;

int    CmpFunc(const void* _a, const void* _b) ;
//double minCondInfoNormalized(int *A_go_indices, short** go_profile, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, int shuffle, double minr, int v, int* midx, char** motifs) ;
int test_minCondInfoNormalized(int *A_go_indices, short** go_profile, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, int shuffle, double minr, int v, int* midx, char** motifs, int idxmain, float* minratio);

float  compare_regexp_motifs(char* m1, char* m2); 

int main(int argc, char ** argv){
  
  FILE*   f;
 char*   outfile;
 int     i,j;
 
//expfile
 char*   expfile;
 int     m, n;
 char**  rownames;
 char**  colnames;
 float** data;
 float*  E;
 int*    E_q;
 float*  E_q_bins;
 
//motiffile
 char*   motiffile;
 int     singlestrand;
 int     nbmotifs;
 int     nbcols;
 char*** motiftable;
 char**  motifs;
 
 float*  M;
 int*    M_q;
 double  motif_mi;
 float   minr    = 5.0; // min imi for accepting a motif
 float   minratio ;
 int     midx;
 
//hash
 ENTRY   e;
 ENTRY   *ep;
 
//fastafile
 int     nborfs;
 char*   fastafile;
 int     nextSequence_started = 0;
 int     nextSequence_ended   = 0;
 char*   nextSequence_currentLine=0; 
 FILE    *fp;
 int     size;
 char*   seq;
 char*   name;
 short** gene_motif;
 int     motif_np;
 
 int      true_idx = 0;
 int      idx;
 char*    realname;
 int      count_seq;
 char*    count_seq_index;
 
 int      quantized   = 0;
 int      ebins       = 0;
 int      mbins       = 2;
 float    divbins     = 50.0;
 MIArray* mi_array ;
 
 int      shuffle     = 10000;
 int*     idx_gene_name;
 int*     A_motif_indices  ;
 int      nb_cat           = 0;
 int      pass = 0;

 if (argc == 1) {
    printf("Usage : mi_motif_compare -expfile FILE -motiffile FILE -fastafile FILE -rna INT -quantized 1/0 -outfile FILE \n");
    exit(0);
  }
  
  expfile         = get_parameter(argc, argv, "-expfile");
  motiffile       = get_parameter(argc, argv, "-motiffile");
  fastafile       = get_parameter(argc, argv, "-fastafile");
  outfile         = get_parameter(argc, argv, "-outfile");
  quantized       = atoi(get_parameter(argc, argv, "-quantized"));

  initialize_nt();
  
  if (exist_parameter(argc, argv, "-rna"))     
    singlestrand    = atoi(get_parameter(argc, argv, "-rna"));

  if (exist_parameter(argc, argv, "-minr"))     
    minr            = atof(get_parameter(argc, argv, "-minr"));

  readFloatTable(expfile, &m, &n, &data, &rownames, &colnames, 0, 1);
  printf("%s loaded: %d by %d\n", expfile, n, m) ;
  
  E              = (float* ) malloc( n * sizeof(float));
  idx_gene_name  = (int*)malloc( n * sizeof(int));
  gene_motif     = (short**) malloc( n * sizeof(short*));

  readStringTable(motiffile, &motiftable,  &nbmotifs, &nbcols) ;
  printf("%s loaded: %d by %d\n\n", motiffile, nbmotifs, nbcols) ;

  if (nbmotifs == 0) {
    die("Please a non-empty motif file.\n");
  }
  motifs = (char**) malloc (nbmotifs * sizeof(char*)) ;
  for(i=0 ; i<nbmotifs ; i++){
    int len = strlen(motiftable[i][0]);
    motifs[i]=(char*) malloc ((len+1) * sizeof(char)) ;
    memcpy(motifs[i],motiftable[i][0],len);
    motifs[i][len]='\0';
    //printf("%d/%d. Motif: %s\n", i+1, nbmotifs, motifs[i]);
  }
  //
  //  create a hash out of the gene names
  //
  hcreate(100000);
  for (i=0; i<n; i++) {
    e.key = strdup(rownames[i]); 
    e.data = (char*)i;
    hsearch(e, ENTER);
  }

  //
  //  read in the sequences that are also in the microarray
  //

  fp = fopen(fastafile, "r");
  if (!fp) {
    printf("cannot open %s ..\n", fastafile);
    exit(0);
  }

  nextSequence_currentLine = (char*)malloc(200000 * sizeof(char));
  nborfs  = -1;
  realname = (char*)calloc(100, sizeof(char));

  count_seq       = -1;
  count_seq_index = (char*)calloc(100000,  sizeof(char));
  printf("\nReading Sequences:") ;
  while ( (seq = nextSequence(fp, &name, &size, &nextSequence_started, &nextSequence_ended, nextSequence_currentLine)) ) {   
     
    count_seq ++; if (count_seq == 100000) die("pb with count_seq > 100000\n");
    //printf(".") ;
    fflush(stdout) ;
    //
    // cut out stuff after the first space
    //
    i = 0;
    while ((name[i] != ' ') && (i >= 0) && (i<strlen(name))) i++;
    strncpy(realname, name, i);
    realname[i] = '\0';
    
    // printf("%s\n", realname);

    //
    // accept the orf if the expression data also has it
    // 
    e.key = realname;
    ep = hsearch(e, FIND);
    if (!ep) {
      free(seq);
      free(name);
      count_seq_index[ count_seq ] = 0;
      continue;   
    } 
     
    count_seq_index[ count_seq ] = 1;
     
    idx = (int)ep->data;       

    // copy expression data into a new array
    E[ true_idx ]     = data[idx][0];
    idx_gene_name[ true_idx ] = idx;  
     
    //						
    // build a kmer profile for this sequence .. use regexp or not ? YES, because k-mers are preselected
    //     
    gene_motif    [ true_idx ] = (short*)calloc(nbmotifs, sizeof(short));
       
    for (i=0; i<nbmotifs; i++) {
      motif_np = 0;
      motif_np =  re_matches(motifs[i], seq, singlestrand); 
      //printf("%s vs %s: %d\n", motifs[i], rownames[idx_gene_name[true_idx]], motif_np);
      //getchar() ;
      gene_motif[ true_idx ][ i ] = (short)motif_np;    
    }

    true_idx ++;
     
    free(seq);
    free(name);
  }
  fclose(fp);
  printf("Completed\n") ;
  fflush(stdout);

  nborfs = true_idx;

  if ((quantized == 0) && (ebins == 0)) {
    ebins = (int)(0.5 + (float)nborfs / ( divbins * mbins ));
  }

  if (quantized == 0) {
    // add a little random number to each value in the E vector
    add_small_values_to_identical_floats(E, nborfs);    
  }

  if (outfile == 0) {
    fp = stdout;
  } else {
    fp = fopen(outfile, "w");
    printf ("\nOutput file (%s) loaded...\n", outfile) ;
    if (fp == 0) {
      printf("cannot open outfile: %s\n", outfile);
    }
  }
   
  //
  //  PROCESS EXPRESSION DATA
  //
  quantize_E(E, nborfs, quantized, &ebins, &E_q, &E_q_bins); 

  //						
  //  PROCESS ALL MOTIFS[/SEED PAIRS]
  //
  mi_array = (MIArray*) malloc (nborfs*sizeof(MIArray)) ;
  for (i=0; i<nbmotifs; i++) {
    M = stof_matrix_column(gene_motif, i, nborfs);
    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

    motif_mi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins);
    printf("MI: E_q vs %s = %f\n", motifs[i], motif_mi) ;
    fflush(stdout);
    mi_array[i].score = motif_mi ;
    mi_array[i].index = i ;
    free(M);
    free(M_q);
  }

  // sort motifs based on MI
  qsort((void*)mi_array, nbmotifs, sizeof(MIArray), CmpFunc);

  // will store accepted motifs (indices)
  A_motif_indices  = (int*)malloc(nbmotifs * sizeof(int));

  f = fopen(outfile, "w") ;
  
  // iterate through motifs
  for (i=0; i<nbmotifs; i++) {

    j = mi_array[i].index;
    
    // get raw profile
    M = stof_matrix_column(gene_motif, j, nborfs);     

    // quantize
    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

    // test using cond info
    if (i != 0) {
      printf("Accept or reject motif %s (index = %d)?\n", motifs[mi_array[i].index], mi_array[i].index);
      minratio = -1;
      pass = test_minCondInfoNormalized(A_motif_indices, gene_motif, mbins, nb_cat, M_q, mbins, E_q, ebins, nborfs, shuffle, minr, 1, &midx, motifs, j, &minratio);
      if (pass == 0){
	continue;
      }
    }
    
    A_motif_indices[ nb_cat ] = mi_array[i].index;
    nb_cat++ ;
    free(M_q);
    printf("Passed: motif %s, mi = %4.3f\n", motifs[mi_array[i].index], mi_array[i].score) ;
    fprintf(f,"%s\n", motifs[mi_array[i].index]) ;
   }

  return 0;
}

int CmpFunc(const void* _a, const void* _b)
{
 const MIArray* a = (const MIArray*) _a;
 const MIArray* b = (const MIArray*) _b;

 if (a->score < b->score)
   return 1;
 else if(a->score == b->score)
   return  0;
 else
   return -1;
}


//
// for each motif already in list, calculate R and p
//   pass = 0
//   if ((R > minr) || (p <= 0.25))
//         pass = 1;

//
int test_minCondInfoNormalized(int *A_go_indices, short** go_profile, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, int shuffle, double minr, int v, int* midx, char** motifs, int idxmain, float* minratio)
{
  int    i, j;
  int*   AB_q;
  int    DAE, DAB;
  double mi_ab_e, mi_a_e, mi_a_b;
  double cmi;
  
  double ratio;
  int*   M_q;
  float* M;
  int    pass = 0;
  float  pe;

  *minratio = 1e6;

  DAE    = DA * DE;
  DAB    = DA * DB;
  *midx  = -1;

  for (i=0; i<nA; i++) {  // nA is the number of already kept GO cats

    j   = A_go_indices[i];  // get index of existing category
    M = stof_matrix_column(go_profile, j, n);
    quantize_M_zero_eqpop(M, n, DA, &M_q);


    // calculate conditional information I(Cnew;E|Cexist)  = I(A;E|B)  =  I(A,B;E) - I(A;E)
    AB_q        = combineQuantizedVectors(M_q, B_q, n, DA, DB);
    mi_ab_e     = CalculateMIbasic       (AB_q,   E_q, n, DAB, DE);
    mi_a_e      = CalculateMIbasic       (M_q, E_q, n, DA, DE);
    cmi         = mi_ab_e - mi_a_e;
   
    // calculate MI I(Cnew;Cexist)
    mi_a_b      = CalculateMIbasic       (M_q, B_q, n, DA, DB);

    if (cmi < 1e-16)
      cmi = 0.0;
    else if (mi_a_b < 1e-16)
      mi_a_b = 1e-16;

    ratio       = cmi / mi_a_b;
    pe          = compare_regexp_motifs(motifs[idxmain], motifs[j]); 

    printf("\t%s (index=%d): %f, p=%f\n", motifs[j], j, ratio, pe);

    if (i == 0)
      *minratio = ratio;
    else {

      if (ratio < *minratio) {
	*minratio = ratio;
      }
    }

    free(M_q);
    free(AB_q);

    if ((ratio > minr) || (pe < 0.50)) {
      
      pass  = 1;
      
    } else {
      *midx = i;
      pass  = 0;
      return pass;
    }




  }

  return pass;

}


float compare_regexp_motifs(char* m1, char* m2) 
{

  float* lwm1;
  float* lwm2;
  int    lw1;
  int    lw2;
  float  p;

  getLinearFloatWMfromRegexp(m1, &lwm1, &lw1);
  getLinearFloatWMfromRegexp(m2, &lwm2, &lw2);
 
  p = compare_motifs(lwm1, lw1, lwm2, lw2, 0, 4); 

  return p;
}
