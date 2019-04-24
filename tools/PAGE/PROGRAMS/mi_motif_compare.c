#define _GNU_SOURCE
#include <search.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>
#include <statistics.h>

#ifdef USEMPATROL
#include <mpatrol.h>
#endif

#include "dataio.h"
#include "statistics.h"
#include "information.h"
#include "mi_library.h"

typedef struct _MIARRAY {
 int    index;      // where is it found
 double score;
} MIArray;

int CmpFunc(const void* _a, const void* _b) ;

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

int*     A_motif_indices  ;
int      nb_cat           = 0;

int main(int argc, char ** argv){
  
  if (argc == 1) {
    printf("Usage : mi_motif_compare -expfile FILE -motiffile FILE -fastafile FILE -rna INT -quantized 1/0 -outfile FILE \n");
    exit(0);
  }
  
  expfile         = get_parameter(argc, argv, "-expfile");
  motiffile       = get_parameter(argc, argv, "-motiffile");
  fastafile       = get_parameter(argc, argv, "-fastafile");
  outfile         = get_parameter(argc, argv, "-outfile");
  singlestrand    = atoi(get_parameter(argc, argv, "-rna"));
  quantized       = atoi(get_parameter(argc, argv, "-quantized"));

  readFloatTable(expfile, &m, &n, &data, &rownames, &colnames, 0, 1);
  E = (float* ) malloc( n * sizeof(float));


  readStringTable(motiffile, &motiftable, &nbcols, &nbmotifs)
  if (nbmotifs == 0) {
    die("Please a non-empty motif file.\n");
  }
  motifs = (char**) malloc (nbmotifs * sizeof(char*)) ;
  for(i=0 ; i<nbmotifs ; i++){
    strcpy(motifs[i],motiftable[i][0]);
    printf("Motif: %s", motifs[i]);
    getchar();
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
  while ( (seq = nextSequence(fp, &name, &size, &nextSequence_started, &nextSequence_ended, nextSequence_currentLine)) ) {   
     
    count_seq ++; if (count_seq == 100000) die("pb with count_seq > 100000\n");
     
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
      gene_motif[ true_idx ][ i ] = (short)motif_np;    
    }

    true_idx ++;
     
    free(seq);
    free(name);
  }
  fclose(fp);

  nborfs = true_idx;

  if ((quantized == 0) && (ebins == 0)) {
    ebins = (int)(0.5 + (float)nborfs / ( divbins * mbins ));
    if (verbose == 1)
     printf("number of bins for expression data %d\n", ebins);
  }

  if (quantized == 0) {
    // add a little random number to each value in the E vector
    add_small_values_to_identical_floats(E, nborfs);    
  }

  if (outfile == 0) {
    fp = stdout;
  } else {
    fp = fopen(outfile, "w");
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
  mi_array = (MIArray*) malloc (sizeof(MIArray)) ;
  for (i=0; i<nbmotifs; i++) {
    
    //
    // START: MOTIF
    //
    M = stof_matrix_column(gene_motif, i, nborfs);     
    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

    motif_mi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins);
    mi_array[i].score = motif_mi ;
    mi_array[i].index = i ;
  }

  qsort((void*)mi_array, nbmotifs, sizeof(MIArray), CmpFunc);
  A_motif_indices  = (int*)malloc(nbmotifs * sizeof(int));

  f = fopen(outfile, "w") ;
  for (i=0; i<nbmotifs; i++) {
    M = stof_matrix_column(gene_motif, i, nborfs);     
    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

    if (i != 0){
      minratio = minCondInfoNormalized(A_motif_indices, gene_motif, mbins, nb_cat, M_q, mbins, E_q, ebins, nborfs, shuffle, minr, 1, &midx);
      if (minratio<minr){
	continue;
      }
    }
    
    A_motif_indices[ nb_cat ] = mi_array[i].index;
    nb_cat++ ;
    free(M_q);
    printf("Passed: motif %s, mi = %4.3f\n", motifs[i], mi_array[i].score) ;
    fprintf(f,"%s\n", motifs[i]) ;
    //printf("Passed: cat %d, mi = %4.3f, pvalue = %4.3f, %s\n", go_array[i].index, go_array[i].score, pass, go_desc[ go_array[nk].index ] );
   }
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
