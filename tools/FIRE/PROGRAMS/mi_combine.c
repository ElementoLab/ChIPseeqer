// takes as input a set of ORFs => cluster ID
//  for each member of this dataset, calculate the k-closest neighbor (Pearson)
// 



#define _GNU_SOURCE
#include <search.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>


#include "dataio.h"
#include "regexp.h"
#include "statistics.h"
#include "prefix.h"
#include "information.h"
#include "mi_library.h"

#define MIN_KMER_COUNT 20

#define DIVBINS 50.0
#define DIVBINS_DIST 10.0               // for position bias 

#define PE_T 0.00

int CmpFunc(const void* _a, const void* _b);
void readKmers_this (char* kmerfile, int max_nb_kmers, int inc_nb_kmers, char*** kmers, int* nbkmers, int* kmersize, int** dna_rna); 

int  testMotifClustering(float* M1_cnt, float* M2_cnt, int*** dist, int a, int b, int* gene_length, int n, float* d_avg);

int  testMotifPairPosition(float* M1_cnt, float* M2_cnt, int*** dist, int a, int b, 
			   int n,
			   int mbins_dist, int shuffle, int* rank, double* z, float* pe);


int  testMotifColocWithSelf(float* E, int quantized, float* M_cnt, 
			   unsigned short*** dist, int a, int n,
			    int shuffle, int divbins, int* rank, double* z, float* p);


int  testMotifsInteraction(float* M1_cnt, float* M2_cnt, int*** dist, int a, int b, 
			   int n,
			   int mbins_dist, int shuffle, int* rank, double* z, float* p);

int  testMotifPairDistance(float* E, int quantized, float* M1_cnt, float* M2_cnt, 
			   unsigned short*** dist, int a, int b, int n,
			   int shuffle, int divbins, int* rank, double* z, float* p);

int  testMotifPairRelativeOrientation(float* E, int quantized, float* M1_cnt, float* M2_cnt, 
				      int*** dist, int*** ori, int a, int b, int n,
				      int shuffle, int divbins, int* rank, double* z, float* p);

int  testMotifPairOrder(float* E, int quantized, float* M1_cnt, float* M2_cnt, 
			int*** dist, int a, int b, int n,
			int shuffle, int divbins, int* rank, double* z, float* p);

int main(int argc, char ** argv)
{
     

  // GENERAL
  int i, j;
  FILE  *fp = 0, *fpo = 0;  ENTRY e;
  ENTRY *ep;
  
  // general parameters
  int    shuffle = 10000;

  // kmer stufff
  char**   kmers;
  char**   kmers_rna = 0;
  char**   kmers_dna = 0;
  
  

  int      kmersize;
  int      nbkmers = 0;
  
  int      nbkmers_rna = 0;
  int      nbkmers_dna = 0;
  
  char*    kmerfile_rna = 0;
  char*    kmerfile_dna = 0;
  char*    kmerfile = 0;

  int      max_nb_kmers = 10000;
  int      inc_nb_kmers = 10000;
   
  // sequence stuff
  int      nborfs = 0;
    
  char*    realname;
  int      size;
  char*    seq;
  char*    name;
  int      nextSequence_started = 0;
  int      nextSequence_ended   = 0;
  char*    nextSequence_currentLine=0; 

  char*    fastafile_dna;
  char*    fastafile_rna;

  int      true_idx = 0;
  int      idx;
  int      l;
  
  // motif stuff
  unsigned short*        kmer_pos;
  int**       gene_kmer_cnt;
  int**       gene_kmer_cnt_masked;
  int         mbins = 2;
  int         no, np;

  unsigned short***      gene_kmer_pos;
  
  //
  // expression stuff
  //
  char*   expfile;
  int     m, n;
  char**  rownames;
  char**  colnames;
  float** data;
  float*  E;
  float  *Ei;
  int     quantized = 0;

  int     shuffle_rank;

  char*   outfile = 0;

  //int**   tmpCounts;

  int     outprm = 0;
  int     seed   = 12345;
    
  extern   int CmpInt();
  char*    outmatrix = 0;
 
 
  int     a, b, c;
  float*  Mcnt_a, *Mcnt_b, *Mcnt_c;
  float*  Mcnt_a_f, *Mcnt_b_f, *Mcnt_c_f;
 
  int*    Mcnt_a_q, *Mcnt_b_q, *Mcnt_c_q;
  double  mi_cnt_ab, mi_cnt_bc;
  int     rank_cnt_ab, rank_cnt_bc; 
  double  z_cnt_ab, z_cnt_bc;
  float   pe_ab, pe_bc;
  int     cnt = 0;
  int     good = 0;
  int*    oldE_newE;
  int*    dna_rna;
  int*    taken;
  int     poscor = 0;
  int**   passed;
  char*   outfullmimatrix = 0;
  FILE*   fpfo;
  int     doallstats = 0;  
  char**  mi_calculated;     
  int*    gene_length_dna;
  int*    gene_length_rna;

  float   minz = -100000.0;
  int     cluster_rank;
  int     mbins_interval_dna = 200;
  int     mbins_interval_rna = 200;
  double  avg_len_dna        = 0;
  double  avg_len_rna        = 0;

  int     rank;
  double  z;  
  float   pe;
  float   divbins     = 50.0;

  float** Edata_tmp;
  float** Edata;
  
  int     verbose = 0;
  int     max_num_sites = 100;

  if (argc == 1) {
    printf("Usage : mi_combine -summaryfile FILE -expfile FILE -quantized INT -fastafile_dna FILE -poscor INT -outmimatrix FILE -outfullmimatrix FILE -shuffle INT -doallstats INT -minz FLOAT -mbins_interval_dna INT -mbins_interval_rna INT\n");
    exit(0);
  }

  expfile         = get_parameter(argc, argv, "-expfile");

  if (exist_parameter(argc, argv, "-summaryfile_dna")) 
    kmerfile_dna    = get_parameter(argc, argv, "-summaryfile_dna");
  
  if (exist_parameter(argc, argv, "-summaryfile_rna")) 
    kmerfile_rna    = get_parameter(argc, argv, "-summaryfile_rna");
  
  if (exist_parameter(argc, argv, "-summaryfile")) 
    kmerfile        = get_parameter(argc, argv, "-summaryfile");
   
  if (exist_parameter(argc, argv, "-outfile")) 
    outfile         = get_parameter(argc, argv, "-outfile");
  
  if (exist_parameter(argc, argv, "-outmatrix")) 
    outmatrix       = get_parameter(argc, argv, "-outmatrix");

  if (exist_parameter(argc, argv, "-outmimatrix")) 
    outmatrix       = get_parameter(argc, argv, "-outmimatrix");
  
  if (exist_parameter(argc, argv, "-outfullmimatrix")) 
    outfullmimatrix = get_parameter(argc, argv, "-outfullmimatrix");
  
  fastafile_dna     = get_parameter(argc, argv, "-fastafile_dna");
  
  fastafile_rna     = get_parameter(argc, argv, "-fastafile_rna");
  
  if (exist_parameter(argc, argv, "-quantized")) 
    quantized       = atoi(get_parameter(argc, argv, "-quantized"));
  
  if (exist_parameter(argc, argv, "-doallstats")) 
    doallstats      = atoi(get_parameter(argc, argv, "-doallstats"));
 
  if (exist_parameter(argc, argv, "-verbose")) 
    verbose      = atoi(get_parameter(argc, argv, "-verbose"));
  
  if (exist_parameter(argc, argv, "-poscor")) 
    poscor          = atoi(get_parameter(argc, argv, "-poscor"));
  
  if (poscor == 1)
    printf("Using positive correlation.\n");

  if (exist_parameter(argc, argv, "-seed")) 
    seed            = atoi(get_parameter(argc, argv, "-seed"));
  default_set_seed(seed);
  
  if (exist_parameter(argc, argv, "-outprm")) 
    outprm          = atoi(get_parameter(argc, argv, "-outprm"));
  
  if (exist_parameter(argc, argv, "-shuffle")) 
    shuffle          = atoi(get_parameter(argc, argv, "-shuffle"));
  
  if (exist_parameter(argc, argv, "-minz")) 
    minz             = atof(get_parameter(argc, argv, "-minz"));

  if (exist_parameter(argc, argv, "-mbins_interval_dna")) 
    mbins_interval_dna = atoi(get_parameter(argc, argv, "-mbins_interval_dna"));
 
  if (exist_parameter(argc, argv, "-mbins_interval_rna")) 
    mbins_interval_rna = atoi(get_parameter(argc, argv, "-mbins_interval_rna"));

  if (exist_parameter(argc, argv, "-divbins")) 
    divbins = atof(get_parameter(argc, argv, "-divbins"));

  if (shuffle_rank == 0) {
    shuffle_rank    = shuffle;  // default
  }
  
  if (kmerfile_dna != 0)
    readKmers_general (kmerfile_dna, max_nb_kmers, inc_nb_kmers, &kmers_dna, &nbkmers_dna, &kmersize); 
  
  if (kmerfile_rna != 0)
    readKmers_general (kmerfile_rna, max_nb_kmers, inc_nb_kmers, &kmers_rna, &nbkmers_rna, &kmersize); 
  
  if (kmerfile != 0) {
    readKmers_this (kmerfile, max_nb_kmers, inc_nb_kmers, &kmers, &nbkmers, &kmersize, &dna_rna); 
    if (nbkmers == 0) {
      die("No motifs in summary file.\n");
    }
  }
   
   if (kmerfile == 0) {
     //
     // pool together kmers
     //
     nbkmers = nbkmers_dna + nbkmers_rna;
     if (nbkmers == 0) {
       die("No motifs in summary file(s).\n");
     }
     dna_rna = (int*)calloc(nbkmers,  sizeof(int));
     kmers   = (char**)malloc(nbkmers * sizeof(char*)); 
     j = 0;
     for (i=0; i<nbkmers_dna; i++) {
       kmers[ j ] = strdup(kmers_dna[i]);
       dna_rna[j] = 0;
       j ++;
     }
     for (i=0; i<nbkmers_rna; i++) {
       kmers[ j ] = strdup(kmers_rna[i]);
       dna_rna[j] = 1;
       j ++;
     }
   } else {

     nbkmers_dna = 0;
     nbkmers_rna = 0;
     for (i=0; i<nbkmers; i++) {
       if (dna_rna[i] == 0)
	 nbkmers_dna ++;
       if (dna_rna[i] == 1)
	 nbkmers_rna ++;
     }

   }

   
   mi_calculated = (char**)calloc( nbkmers, sizeof(char**));
   for (a=0; a<nbkmers; a++) {
     mi_calculated[a] = (char*)calloc( nbkmers, sizeof(char));
   }

   if (doallstats == 1)
     printf("Doing all stats.\n");

   printf("Summary file has %d DNA motifs and %d RNA motifs.\n", nbkmers_dna, nbkmers_rna);

     

   //
   //  read in expression data
   //
   readFloatTable(expfile, &m, &n, &data, &rownames, &colnames, 0, 1);

   
   E                 = (float* ) malloc( n * sizeof(float));      
   Edata_tmp         = (float**) malloc( n * sizeof(float*));

   gene_kmer_cnt     = (int**)   malloc( n * sizeof(int*));
   gene_kmer_pos     = (unsigned short***)  malloc( n * sizeof(unsigned short**));

   gene_length_dna   = (int*)    malloc( n * sizeof(int));
   gene_length_rna   = (int*)    malloc( n * sizeof(int));

   oldE_newE         = (int*)    malloc( n * sizeof(int));
   for (i=0; i<n; i++) {
     oldE_newE[i] = -1;
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
   //
   //  DNA DNA DNA DNA read in the DNA sequences that are also in the microarray
   //
   //
   nextSequence_currentLine = (char*)malloc(200000 * sizeof(char));
	

   if (nbkmers_dna > 0) {
     
     printf("Read in DNA sequences and build motif profiles ...");

     fp = fopen(fastafile_dna, "r");
     if (!fp) {
       printf("cannot open %s ..\n", fastafile_dna);
       exit(0);
     }
     
     nborfs      = 0;
     avg_len_dna = 0;

     realname = (char*)calloc(100, sizeof(char));
     while ( (seq = nextSequence(fp, &name, &size, &nextSequence_started, &nextSequence_ended, nextSequence_currentLine)) ) {   
       
       //
       // cut out stuff after the first space
       //
       i = 0;
       while ((name[i] != ' ') && (i >= 0) && (i<strlen(name))) i++;
       strncpy(realname, name, i);
       realname[i] = '\0';
       
       
       //
       // accept the orf if the expression data also has it
       // 
       e.key = realname;
       ep = hsearch(e, FIND);
       
       if (!ep) {
	 free(seq);
	 free(name);
	 continue;
       } 
       
       // get the index for the gene
       idx = (int)ep->data;       
       
       // update expression vector
       E[ true_idx ]             = data[idx][0];
       Edata_tmp[ true_idx ]     = data[idx];
     
       oldE_newE[ idx ] = true_idx;
       
       gene_kmer_cnt   [ true_idx ] = (int*)  calloc(nbkmers,  sizeof(int));
       gene_kmer_pos   [ true_idx ] = (unsigned short**) malloc(nbkmers * sizeof(unsigned short*));

       gene_length_dna [ true_idx ] = strlen(seq);

       l  = strlen(seq);
       avg_len_dna += l;

       // create position profile
       for (i=0; i<nbkmers; i++) {
	 
	 if (dna_rna[i] == 1)
	   continue;
	 
	 no = 0; np = 0;
	 
	 kmer_pos = (unsigned short*)malloc(max_num_sites * sizeof(unsigned short));

	 if (kmer_pos == 0) {
	   die("sorry can't allocate mem for kmer_pos (DNA).\n");
	 }
	 
	 findSites_micombine(kmers[i], seq, kmer_pos, &np, 0, 0, max_num_sites, 0, 1, 0) ;
	 
	 gene_kmer_cnt[ true_idx ][ i ]    = np;	 
	 
	 //if (np >= 2)
	 //  printf("Found %d copies in %s\n", np, name);
	 
	 realloc(kmer_pos, np * sizeof(unsigned short)); 

	 gene_kmer_pos[ true_idx ][ i ]    = kmer_pos;

	 
       }
       
       true_idx ++;
       
       free(seq);
       free(name);
       
     }
     
     printf("Done.\n");

     nborfs = true_idx;
   
     fclose(fp);
     
     avg_len_dna /= nborfs;
     
   }

   //
   //  END reading DNA DNA DNA DNA sequences
   //
   //


   //
   //  START reading RNA RNA RNA RNA sequences
   //

   if (nbkmers_rna > 0) {
     

     printf("Read in RNA sequences and build motif profiles ...");
          
     //
     //  START: read in the *** RNA *** sequences that are also in the microarray
     //
     true_idx = 0;
     
     fp = fopen(fastafile_rna, "r");
     if (!fp) {
       printf("cannot open %s ..\n", fastafile_rna);
       exit(0);
     }
     
     // reset
     nextSequence_started = 0;
     nextSequence_ended   = 0;
     
     avg_len_rna = 0;

     realname = (char*)calloc(100, sizeof(char));
     while ( (seq = nextSequence(fp, &name, &size, &nextSequence_started, &nextSequence_ended, nextSequence_currentLine)) ) {   
       
       //
       // cut out stuff after the first space
       //
       i = 0;
       while ((name[i] != ' ') && (i >= 0) && (i<strlen(name))) i++;
       strncpy(realname, name, i);
       realname[i] = '\0';
       
       
       //
       // accept the orf if the expression data also has it
       // 
       e.key = realname;
       ep = hsearch(e, FIND);
       
       if (!ep) {
	 free(seq);
	 free(name);
	 continue;
       } 
       
 

       // get the index for the gene
       idx = (int)ep->data;       

       gene_length_rna [ true_idx ] = strlen(seq);
       
       // get the expression index
       if (nbkmers_dna == 0) {
	 E[ true_idx ]             = data[idx][0];
	 Edata_tmp[ true_idx ]     = data[idx];

	 gene_kmer_cnt[ true_idx ] = (int*)  calloc(nbkmers, sizeof(int));
	 gene_kmer_pos[ true_idx ] = (unsigned short**) malloc(nbkmers * sizeof(unsigned short*));

       } else {
	 true_idx                  = oldE_newE[ idx ];

       }
 
       // added here that the DNA fasta file must have it too
       // if, for some reason, true_idx = -1, then we cannot consider this sequence
       if (true_idx == -1) {
	 free(seq);
	 free(name);
	 continue;  
       }
       

       l           = strlen(seq);
       avg_len_rna += l;

       //
       // create position profile
       //
       for (i=0; i<nbkmers; i++) {
	 
	 if (dna_rna[i] == 0)
	   continue;
	 
	 np       = 0;      
	 kmer_pos = (unsigned short*)malloc(max_num_sites * sizeof(unsigned short));
	 
	 findSites_micombine(kmers[i], seq, kmer_pos, &np, 0, 0, max_num_sites, 1, 1, 0) ;
	 
	 gene_kmer_cnt[ true_idx ][ i ] = np;


	 realloc(kmer_pos, np * sizeof(unsigned short));
	 gene_kmer_pos[ true_idx ][ i ] = kmer_pos;

       }
      
       if (nbkmers_dna == 0)
	 true_idx ++;
 
       avg_len_rna /= nborfs;

       free(seq);
       free(name);
       
     }
     
     fclose(fp);
     
     if (nbkmers_dna == 0)
       nborfs = true_idx;
     
     printf("Done.\n");

   }

  
   //						
   // here figure out wich motifs are palindromic
   //
     
   if (outfile != 0) {
     fp = fopen(outfile, "w"); 
   }

   if (outmatrix != 0) {
     fpo = fopen(outmatrix, "w");
     if (!fpo) {
       die("Cannot open outmatrix.\n");
     }
   }
   
  
   //
   //  apply mask in gene_kmer_cnt
   //
   
   //
   //  we set to zero certain motif counts, if gene is not an over-rep cluster
   //
   gene_kmer_cnt_masked = (int**)malloc( nborfs * sizeof(int*) );
   for (i=0; i<nborfs; i++) {
     gene_kmer_cnt_masked[i] = (int*)malloc( nbkmers * sizeof(int));
     for (a=0; a<nbkmers; a++) {
       if (Edata_tmp[i][a] == 0)
	 gene_kmer_cnt_masked[i][a] = 0;
       else
	 gene_kmer_cnt_masked[i][a] = gene_kmer_cnt[i][a];
     }
   }
   
  //
  // transpose E matrix
  // 
  transpose_f(Edata_tmp, nborfs, m, &Edata);
  
  
  if (outfullmimatrix != 0) {
     
     printf("Outputing MI matrix ... ");
     
     fpfo = fopen(outfullmimatrix, "w");
     if (!fpfo) {
       die("Cannot open outfullmimatrix.\n");
     }
   
     
     for (a=0; a<nbkmers-1; a++) {
       
       Mcnt_a = itof_matrix_column(gene_kmer_cnt_masked, a, nborfs);
       quantize_M_zero_eqpop(Mcnt_a, nborfs, mbins, &Mcnt_a_q); 
       
       for (b=a+1; b<nbkmers; b++) {
	 
	 Mcnt_b = itof_matrix_column(gene_kmer_cnt_masked, b, nborfs);
	 quantize_M_zero_eqpop(Mcnt_b, nborfs, mbins, &Mcnt_b_q); 	 

	 mi_cnt_ab   = CalculateMIbasic(Mcnt_a_q, Mcnt_b_q, nborfs, mbins, mbins);	 
	 pe_ab       = pearson(Mcnt_a, Mcnt_b, nborfs);	 
	 fprintf(fpfo, "%s\t%s\t%5.4f\tnan\tnan\t%5.4f\n", kmers[a], kmers[b], mi_cnt_ab, pe_ab);
	 
       }
       
     }

     fclose(fpfo);
     
     printf("Done.\n");
     
   } // if (outfullmimatrix != 0)

   //						
   //  END: reading RNA sequences
   //


   
   taken = (int*)calloc(nbkmers, sizeof(int));
   for (i=0; i<nbkmers; i++) {
     taken[i] = -1;
   }

   passed = (int**)malloc(nbkmers * sizeof(int*));
   for (i=0; i<nbkmers; i++) {
     passed[i] = (int*)malloc(nbkmers * sizeof(int));
     for (j=0; j<nbkmers; j++) {
       passed[i][j] = -1;  // don't know
     }
   }
   

   cnt = -1;
   a   = 0;
   while (a < nbkmers) {
     
     //
     // get the ath column of the gene_kmer matrix
     //

     if (taken[a] == -1) {
       
       taken[a] = ++cnt;

       if (outfile != 0) {       
	 fprintf(fp, "%s\t%d\n", kmers[a], cnt);
       }
       printf("%s\t%d\n", kmers[a], cnt);
       
       Mcnt_a_f = itof_matrix_column(gene_kmer_cnt,        a, nborfs);
       Mcnt_a   = itof_matrix_column(gene_kmer_cnt_masked, a, nborfs);
       quantize_M_zero_eqpop(Mcnt_a, nborfs, mbins, &Mcnt_a_q); 

       
       b = a + 1;
       while (b < nbkmers) {

	 // continue only if b is not taken
	 if (taken[b] == -1) {

	   //
	   // get the bth column of the gene_kmer matrix
	   //
	   Mcnt_b_f = itof_matrix_column(gene_kmer_cnt,        b, nborfs);
	   Mcnt_b   = itof_matrix_column(gene_kmer_cnt_masked, b, nborfs);
	   quantize_M_zero_eqpop(Mcnt_b, nborfs, mbins, &Mcnt_b_q); 
	   mi_cnt_ab   = CalculateMIbasic(Mcnt_a_q, Mcnt_b_q, nborfs, mbins, mbins);
	   max_rank_test (mi_cnt_ab, Mcnt_a_q, mbins, Mcnt_b_q, mbins, nborfs, shuffle, 1, &rank_cnt_ab, 1, &z_cnt_ab); 
	   pe_ab       = pearson(Mcnt_a, Mcnt_b, nborfs);
	   
	   if (outmatrix != 0)
	     fprintf(fpo, "%s\t%s\t%5.4f\t%d\t%5.4f\t%5.4f", kmers[a], kmers[b], mi_cnt_ab, rank_cnt_ab, z_cnt_ab, pe_ab);

	   printf("%s\t%s\t%5.4f\t%d\t%5.4f\t%5.4f", kmers[a], kmers[b], mi_cnt_ab, rank_cnt_ab, z_cnt_ab, pe_ab);

	   if ((rank_cnt_ab == 0) && (z_cnt_ab > minz) && (pe_ab > PE_T) && (dna_rna[a] == dna_rna[b])) {
	     
	     //cluster_rank = testMotifsInteraction(Mcnt_a, Mcnt_b, gene_kmer_pos, a, b, nborfs, (dna_rna[a]==0?mbins_dist_dna:mbins_dist_rna), shuffle, &rank, &z, &pe);

	     Ei = intersect_binary_vector_f(Edata[a], Edata[b], nborfs);
	     //Eu = union_binary_vector_f(Edata[a], Edata[b], nborfs);
	     
	     //
	     //  we use the original profiles, here, not the ones based on over-rep
	     //  why ? because we need a background. Imagine the case of a cluster where all motifs are
	     //          close to each other; the entropy would be zero
	     // 
	     cluster_rank = testMotifPairDistance(Ei, quantized, Mcnt_a_f, Mcnt_b_f, gene_kmer_pos, a, b, nborfs, shuffle, divbins, &rank, &z, &pe);
	     fprintf(fpo, "\t%d\t%5.4f", rank, z);
	     printf("\t%d\t%5.4f", rank, z);

	     //cluster_rank = testMotifPairDistance(Eu, quantized, Mcnt_a_f, Mcnt_b_f, gene_kmer_pos, a, b, nborfs, shuffle, divbins, &rank, &z, &pe);
	     fprintf(fpo, "\t%d\t%5.4f", rank, z);
	     printf("\t%d\t%5.4f", rank, z);

	     fprintf(fpo, "\n");
	     printf("\n");

	     free(Ei);
	     //free(Eu);
	     
	     

	     
	   } else {

	     fprintf(fpo, "\tnan\tnan\tnan\tnan\n"); //\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\n");
	     printf(      "\tnan\tnan\tnan\tnan\n"); //\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\n");

	   }

	   mi_calculated[a][b] = 1;
	   mi_calculated[b][a] = 1;
	   
	   

	   if ((rank_cnt_ab == 0) && (z_cnt_ab > minz) && ( (poscor == 0) || ((poscor == 1) && (pe_ab > PE_T))))  {

	     //
	     //  make sure that this guy is connected to all the members of the current cluster, ie a < x < b				     //

	    
	     
	     good = 1;
	     for (c=a+1; c<=b-1; c++) {
	       
	       //	       printf("c=%d\n", c);

	       if (taken[c] == cnt) {

		 Mcnt_c   = itof_matrix_column(gene_kmer_cnt_masked, c, nborfs);
		 Mcnt_c_f = itof_matrix_column(gene_kmer_cnt,        c, nborfs);
       
		 quantize_M_zero_eqpop(Mcnt_c, nborfs, mbins, &Mcnt_c_q); 

		 mi_cnt_bc   = CalculateMIbasic(Mcnt_b_q, Mcnt_c_q, nborfs, mbins, mbins);
		 max_rank_test (mi_cnt_bc, Mcnt_b_q, mbins, Mcnt_c_q, mbins, nborfs, shuffle, 1, &rank_cnt_bc, 1, &z_cnt_bc); 
		 pe_bc       = pearson(Mcnt_b, Mcnt_c, nborfs);

		 if (outmatrix != 0)
		   fprintf(fpo, "%s\t%s\t%5.4f\t%d\t%5.4f\t%5.4f", kmers[b], kmers[c], mi_cnt_bc, rank_cnt_bc, z_cnt_bc, pe_bc);

		 printf("%s\t%s\t%5.4f\t%d\t%5.4f\t%5.4f", kmers[b], kmers[c], mi_cnt_bc, rank_cnt_bc, z_cnt_bc, pe_bc);
				 
		 if ((rank_cnt_bc == 0) && (z_cnt_bc > minz) && (pe_bc > PE_T) && (dna_rna[b] == dna_rna[c])) {

		   Ei = intersect_binary_vector_f(Edata[b], Edata[c], nborfs);
		   //Eu = union_binary_vector_f(Edata[b], Edata[c], nborfs);

		   cluster_rank = testMotifPairDistance(Ei, quantized, Mcnt_b_f, Mcnt_c_f, gene_kmer_pos, b, c, nborfs, shuffle, divbins, &rank, &z, &pe);

		   fprintf(fpo, "\t%d\t%5.4f", rank, z);
		   printf("\t%d\t%5.4f", rank, z);

		   //cluster_rank = testMotifPairDistance(Eu, quantized, Mcnt_b_f, Mcnt_c_f, gene_kmer_pos, b, c, nborfs, shuffle, divbins, &rank, &z, &pe);
		   fprintf(fpo, "\t%d\t%5.4f", rank, z);
		   printf("\t%d\t%5.4f", rank, z);

		   fprintf(fpo, "\n");
		   printf("\n");

		   free(Ei);
		   //free(Eu);
		   
	     		   
		 } else {

		   fprintf(fpo, "\tnan\tnan\tnan\tnan\n"); //\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\n");
		   printf(      "\tnan\tnan\tnan\tnan\n"); //\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\n");

		 }
  
		 mi_calculated[b][c] = 1;
		 mi_calculated[c][b] = 1;
		 
		 if ((rank_cnt_bc == 0) && (z_cnt_bc > minz) && ( (poscor == 0) || ((poscor == 1) && (pe_bc > PE_T))))  {
		   // good one, do nothing
		 } else {
		   good = 0;
		   break;	
		 }
		 
	       } // endif

	     }

	     if (outfile != 0) {
	       
	       if (good == 1) {
		 fprintf(fp, "%s\t%d\n", kmers[b], cnt);
		 printf("%s\t%d\n", kmers[b], cnt);		      
		 taken[b] = cnt;	 
	       }
      
	     }  // endif
	       
	   } //endif rank_cnt_ab == 0

	 } // endif taken[b] == -1
	 
	 b++;
       }
       
       
     }

     a++;
   }
   
   if (outfile != 0)
     fclose(fp);

   if ((outmatrix != 0) && (doallstats == 1)) {
     
     printf("Finish calculated rank, z-scores for remaining entries ... ");

     for (a=0; a<nbkmers; a++) {
       
       Mcnt_a   = itof_matrix_column(gene_kmer_cnt_masked, a, nborfs);
       Mcnt_a_f = itof_matrix_column(gene_kmer_cnt,        a, nborfs);
       
       quantize_M_zero_eqpop(Mcnt_a, nborfs, mbins, &Mcnt_a_q); 
       
       for (b=a; b<nbkmers; b++) {

	 if (mi_calculated[a][b] == 0) {
	   
	   Mcnt_b   = itof_matrix_column(gene_kmer_cnt_masked, b, nborfs);
	   Mcnt_b_f = itof_matrix_column(gene_kmer_cnt,        b, nborfs);

	   quantize_M_zero_eqpop(Mcnt_b, nborfs, mbins, &Mcnt_b_q); 
	   
	   mi_cnt_ab   = CalculateMIbasic(Mcnt_a_q, Mcnt_b_q, nborfs, mbins, mbins);
	   max_rank_test (mi_cnt_ab, Mcnt_a_q, mbins, Mcnt_b_q, mbins, nborfs, shuffle, 1, &rank_cnt_ab, 1, &z_cnt_ab); 
	   pe_ab       = pearson(Mcnt_a, Mcnt_b, nborfs);
	   
	   fprintf(fpo, "%s\t%s\t%5.4f\t%d\t%5.4f\t%5.4f", kmers[a], kmers[b], mi_cnt_ab, rank_cnt_ab, z_cnt_ab, pe_ab);

	   printf("%s\t%s\t%5.4f\t%d\t%5.4f\t%5.4f", kmers[a], kmers[b], mi_cnt_ab, rank_cnt_ab, z_cnt_ab, pe_ab);
	   

	   if ((rank_cnt_ab == 0) && (z_cnt_ab > minz) && (pe_ab > PE_T) && (dna_rna[a] == dna_rna[b])) {
	     
	     Ei = intersect_binary_vector_f(Edata[a], Edata[b], nborfs);
	     
	     if (a != b) {
	       testMotifPairDistance(Ei, quantized, Mcnt_a_f, Mcnt_b_f, gene_kmer_pos, a, b, nborfs, shuffle, divbins, &rank, &z, &pe);
	     } else {	       
	       testMotifColocWithSelf(Ei, quantized, Mcnt_a_f, gene_kmer_pos, a, nborfs, shuffle, divbins, &rank, &z, &pe);
	     }

	     fprintf(fpo, "\t%d\t%5.4f", rank, z);
	     printf("\t%d\t%5.4f", rank, z);	     

	     //cluster_rank = testMotifPairDistance(Eu, quantized, Mcnt_a_f, Mcnt_b_f, gene_kmer_pos, a, b, nborfs, shuffle, divbins, &rank, &z, &pe);
	     fprintf(fpo, "\t%d\t%5.4f", rank, z);
	     printf("\t%d\t%5.4f", rank, z);

	     fprintf(fpo, "\n");
	     printf("\n");



	     //free(Eu);		   
	     free(Ei);

	   } else {

	     fprintf(fpo, "\tnan\tnan\tnan\tnan\n"); //\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\n");
	     printf(      "\tnan\tnan\tnan\tnan\n"); //\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\n");

	   }
	 }
       }
     }
     
     printf("Done.\n");
   }

   if (outmatrix != 0)
     fclose(fpo);

   return 0;
}






//
// test MI(dist bet M1,M2;E)
//
int  testMotifPairRelativeOrientation(float* E, int quantized, float* M1_cnt, float* M2_cnt, 
				      int*** dist, int*** ori, int a, int b, int n,
				      int shuffle, int divbins, int* rank, double* z, float* p)
{
  int     i,j,k;
  int     nbgenes;
  double  mi_med;
  int     rank_med;
  double  z_med;
  float*  D;
  int*    O_q;
  float   d;
  int     mbins_dist = 2;
  int     ebins_dist;
  float*  E_dist;
  int*    E_q_dist;
  float*  E_q_bins_dist;
  
  E_dist = (float*)malloc(n * sizeof(float));
  D      = (float*)malloc(n * sizeof(float));
  O_q    = (int*)malloc(n * sizeof(int));
  nbgenes = 0;
  for (i=0; i<n; i++) {
    if ((M1_cnt[i] > 0) && (M2_cnt[i] > 0)) {
      E_dist[nbgenes] = E[i];
      D     [nbgenes] = 100000;
      for (j=0; j<(int)(M1_cnt[i]); j++) {
	for (k=0; k<(int)(M2_cnt[i]); k++) {
	  d = fabs( dist[i][a][j] - dist[i][b][k] );
	  if (d < D[nbgenes]) {
	    D[nbgenes] = d;
	    
	    if (ori[i][a][j] == ori[i][b][k]) {
	      O_q[nbgenes] = 1;	      
	    } else {
	      O_q[nbgenes] = 0;
	    }
	    
	  }
	}
      }
      nbgenes  ++;
    }
    
  }

  if (quantized == 0) {
    ebins_dist = (int)(0.5 + nbgenes / ( mbins_dist * divbins ));
  }

  quantize_E(E_dist, nbgenes, quantized, &ebins_dist,  &E_q_dist, &E_q_bins_dist); 
  
  mi_med = CalculateMIbasic(O_q, E_q_dist, nbgenes, mbins_dist, ebins_dist);       
  max_rank_test (mi_med, O_q, mbins_dist, E_q_dist, ebins_dist, nbgenes, shuffle, 1, &rank_med, 1, &z_med); 

  *z    = z_med;
  *rank = rank_med;
  *p    = -1000.0;
  
  return (rank_med == 0);

}



//
// test MI(dist bet M1,M2;E)
//
int  testMotifPairOrder(float* E, int quantized, float* M1_cnt, float* M2_cnt, 
			int*** dist, int a, int b, int n,
			int shuffle, int divbins, int* rank, double* z, float* p)
{
  int     i,j,k;
  int     nbgenes;
  double  mi_med;
  int     rank_med;
  double  z_med;
  float*  D;
  float   d;
  int     mbins_dist;
  int     ebins_dist;
  float*  E_dist;
  int*    E_q_dist;
  float*  E_q_bins_dist;

  int*    RO;  // relative order on the chr
    
  E_dist = (float*)malloc(n * sizeof(float));
  D      = (float*)malloc(n * sizeof(float));
  RO     = (int*)  malloc(n * sizeof(int  ));
  

  nbgenes = 0;
  for (i=0; i<n; i++) {
    if ((M1_cnt[i] > 0) && (M2_cnt[i] > 0)) {
      E_dist[nbgenes] = E[i];
      D     [nbgenes] = 100000;
      for (j=0; j<(int)(M1_cnt[i]); j++) {
	for (k=0; k<(int)(M2_cnt[i]); k++) {
	  d = dist[i][a][j] - dist[i][b][k];
	  if (fabs(d) < fabs(D[nbgenes])) {
	    D[nbgenes] = d;
	  }
	}
      }
      
      RO[nbgenes] = (D[nbgenes]>0?1:0);

      nbgenes  ++;
    }
    
  }
  
  mbins_dist = 2;

  if (quantized == 0) {
    ebins_dist = (int)(0.5 + nbgenes / ( mbins_dist * divbins ));
  }
  
  quantize_E(E_dist, nbgenes, quantized, &ebins_dist,  &E_q_dist, &E_q_bins_dist); 
  
  mi_med = CalculateMIbasic(RO, E_q_dist, nbgenes, mbins_dist, ebins_dist);       
  max_rank_test (mi_med, RO, mbins_dist, E_q_dist, ebins_dist, nbgenes, shuffle, 1, &rank_med, 1, &z_med); 
  
  *z    = z_med;
  *rank = rank_med;
  *p    = -1000.0;
  
  return (rank_med == 0);

}


//
// test MI(dist bet M1,M2;E)
//
int  testMotifColocWithSelf(float* E, int quantized, float* M_cnt, 
			   unsigned short*** dist, int a, int n,
			   int shuffle, int divbins, int* rank, double* z, float* p)
{
  int     i,j,k;
  int     nbgenes;
  double  mi_med;
  int     rank_med;
  int     best_rank_med;
  double  z_med;
  double  best_z_med;
  float*  D;
  int*    D_q;
  float*  D_q_bins;
  float   d;
  int     mbins_dist;
  int     ebins_dist;
  float*  E_dist;
  int*    E_q_dist;
  float*  E_q_bins_dist;
    
  E_dist = (float*)malloc(n * sizeof(float));
  D      = (float*)malloc(n * sizeof(float));

  nbgenes = 0;
  for (i=0; i<n; i++) {
    if (M_cnt[i] >= 2) { 
      //printf("Gene %d has >= 2 copies\n", i);
      E_dist[nbgenes] = E[i];
      D     [nbgenes] = 100000;
      for (j=0; j<(int)(M_cnt[i])-1; j++) {
	for (k=j+1; k<(int)(M_cnt[i]); k++) {
	  d = fabs( dist[i][a][j] - dist[i][a][k] );
	  if (d < D[nbgenes]) {
	    D[nbgenes] = d;
	  }
	}
      }
      //printf("d=%f\n", D[nbgenes]);
      nbgenes  ++;
    }
    
  }

  
  
  best_rank_med = shuffle + 100;
  best_z_med    = -100;
  
  if (nbgenes >= 5) {  // 5 is arbitrary, but doesn't matter
    
   add_small_values_to_identical_floats(D, nbgenes);	  
   
   
   for (mbins_dist=2; mbins_dist<=5; mbins_dist++) {
     
     D_q   = Quantize(D, nbgenes, mbins_dist, &D_q_bins); 
     
     if (quantized == 0) {
       ebins_dist = (int)(0.5 + nbgenes / ( mbins_dist * DIVBINS_DIST ));
     }
     
     quantize_E(E_dist, nbgenes, quantized, &ebins_dist,  &E_q_dist, &E_q_bins_dist); 
     
     mi_med = CalculateMIbasic(D_q, E_q_dist, nbgenes, mbins_dist, ebins_dist);    
     
     //CalculateAndShowCountMatrix(D_q,  E_q_dist, nbgenes, mbins_dist, ebins_dist, 0, 0, 0);
     
     max_rank_test (mi_med, D_q, mbins_dist, E_q_dist, ebins_dist, nbgenes, shuffle, 1, &rank_med, 1, &z_med); 
     
     if ((rank_med < best_rank_med) || ((rank_med == best_rank_med) && (rank_med == 0) && (z_med > best_z_med))) {
       best_z_med    = z_med;
       best_rank_med = rank_med;      
     }
     
     free(D_q);
     free(E_q_dist);
     //free(E_q_bins_dist);

   }
   
  }

  *z    = best_z_med;
  *rank = best_rank_med;
  *p    = -1000.0;
  
  free(D);
  free(E_dist);

  return (best_rank_med == 0);

}




//
// test MI(dist bet M1,M2;E)
//
int  testMotifPairDistance(float* E, int quantized, float* M1_cnt, float* M2_cnt, 
			   unsigned short*** dist, int a, int b, int n,
			   int shuffle, int divbins, int* rank, double* z, float* p)
{
  int     i,j,k;
  int     nbgenes;
  double  mi_med;
  int     rank_med;
  int     best_rank_med;
  double  z_med;
  double  best_z_med;
  float*  D;
  int*    D_q;
  float*  D_q_bins;
  float   d;
  int     mbins_dist;
  int     ebins_dist;
  float*  E_dist;
  int*    E_q_dist;
  float*  E_q_bins_dist;
    
  E_dist = (float*)malloc(n * sizeof(float));
  D      = (float*)malloc(n * sizeof(float));

  nbgenes = 0;
  for (i=0; i<n; i++) {
    if ((M1_cnt[i] > 0) && (M2_cnt[i] > 0)) {
      E_dist[nbgenes] = E[i];
      D     [nbgenes] = 100000;
      for (j=0; j<(int)(M1_cnt[i]); j++) {
	for (k=0; k<(int)(M2_cnt[i]); k++) {
	  d = fabs( dist[i][a][j] - dist[i][b][k] );
	  if (d < D[nbgenes]) {
	    D[nbgenes] = d;
	  }
	}
      }
      nbgenes  ++;
    }
    
  }
  
  add_small_values_to_identical_floats(D, nbgenes);	  
  
  best_rank_med = shuffle + 100;
  best_z_med    = -100;
  
  for (mbins_dist=2; mbins_dist<=5; mbins_dist++) {
    
    D_q   = Quantize(D, nbgenes, mbins_dist, &D_q_bins); 

    if (quantized == 0) {
      ebins_dist = (int)(0.5 + nbgenes / ( mbins_dist * DIVBINS_DIST ));
    }

    quantize_E(E_dist, nbgenes, quantized, &ebins_dist,  &E_q_dist, &E_q_bins_dist); 

    mi_med = CalculateMIbasic(D_q, E_q_dist, nbgenes, mbins_dist, ebins_dist);       
    max_rank_test (mi_med, D_q, mbins_dist, E_q_dist, ebins_dist, nbgenes, shuffle, 1, &rank_med, 1, &z_med); 
    
    if ((rank_med < best_rank_med) || ((rank_med == best_rank_med) && (rank_med == 0) && (z_med > best_z_med))) {
      best_z_med    = z_med;
      best_rank_med = rank_med;      
    }

  }

  *z    = best_z_med;
  *rank = best_rank_med;
  *p    = -1000.0;
  
  return (best_rank_med == 0);

}





//
// 
//
int  testMotifsInteraction(float* M1_cnt, float* M2_cnt, int*** dist, int a, int b, 
			   int n,
			   int mbins_dist, int shuffle, int* rank, double* z, float* p)
{
  int    i, k;
  float  med1, med2;
  float* MM1;
  int*   MM1_q;
  float* MM2;
  int*   MM2_q;
  float* MM1_q_bins;
  float* MM2_q_bins;

  int    nbgenes;
  double mi_med;
  int    rank_med;
  int    best_rank_med;
  double z_med;
  double best_z_med;
  float  pe;

  //
  // determine the median distance for the two motifs
  //

  MM1   = (float*)malloc(n * sizeof(float));
  MM2   = (float*)malloc(n * sizeof(float));

  nbgenes = 0;
  for (i=0; i<n; i++) {
    if ((M1_cnt[i] > 0) && (M2_cnt[i] > 0)) {
      med1       = average_int(dist[i][a], (int)(M1_cnt[i]));
      med2       = average_int(dist[i][b], (int)(M2_cnt[i]));
      MM1[nbgenes] = med1;
      MM2[nbgenes] = med2;
      nbgenes  ++;
    }
  }
  
  pe       = pearson(MM1, MM2, nbgenes);
  *p = pe;
  
  add_small_values_to_identical_floats(MM1, nbgenes);	  
  add_small_values_to_identical_floats(MM2, nbgenes);	  
  
  //printf("quantizing using %d bins.\n", mbins_dist);

  best_rank_med = shuffle + 100;
  best_z_med    = -100;
  for (mbins_dist=2; mbins_dist<=5; mbins_dist++) {

    MM1_q   = Quantize(MM1, nbgenes, mbins_dist, &MM1_q_bins); 
    MM2_q   = (int*)malloc(nbgenes * sizeof(int));
    
    //for (k=0; k<=mbins_dist; k++) {
    //  printf("O=%4.3f\n", MM1_q_bins[k]);
    //}

    for (i=0; i<nbgenes; i++) {
      for (k=0; k<mbins_dist; k++) {
	if (MM2[i] <= MM1_q_bins[k+1]) {
	  MM2_q[i] = k;
	  break;
	}
      }
      if (MM2[i] > MM1_q_bins[mbins_dist]) {
	MM2_q[i] = mbins_dist-1;
      }
      //printf("%5.4f\t%d\n", MM2[i], MM2_q[i]);
    }
    
    // calc MI
    mi_med  = CalculateMIbasic(MM1_q, MM2_q, nbgenes, mbins_dist, mbins_dist);       
    max_rank_test (mi_med, MM1_q, mbins_dist, MM2_q, mbins_dist, nbgenes, shuffle, 1, &rank_med, 1, &z_med); 
  
    if ((rank_med < best_rank_med) || ((rank_med == best_rank_med) && (rank_med == 0) && (z_med > best_z_med))) {
      best_z_med    = z_med;
      best_rank_med = rank_med;      
    }
    
    free(MM1_q);
    free(MM2_q);

    // 
    // other way
    //

    
    MM1_q   = (int*)malloc(nbgenes * sizeof(int));
    MM2_q   = Quantize(MM2, nbgenes, mbins_dist, &MM2_q_bins); 

    for (i=0; i<nbgenes; i++) {
      for (k=0; k<mbins_dist; k++) {
	if (MM1[i] <= MM2_q_bins[k+1]) {
	  MM1_q[i] = k;
	  break;
	}
      }
      if (MM1[i] > MM2_q_bins[mbins_dist]) {
	MM1_q[i] = mbins_dist-1;
      }
    }
    
    // calc MI
    mi_med  = CalculateMIbasic(MM1_q, MM2_q, nbgenes, mbins_dist, mbins_dist);       
    max_rank_test (mi_med, MM1_q, mbins_dist, MM2_q, mbins_dist, nbgenes, shuffle, 1, &rank_med, 1, &z_med); 
  
    if ((rank_med < best_rank_med) || ((rank_med == best_rank_med) && (rank_med == 0) && (z_med > best_z_med))) {
      best_z_med    = z_med;
      best_rank_med = rank_med;      
    }
    
    free(MM1_q);
    free(MM2_q);    

    
  }



  free(MM1); 
  free(MM2);
 

  *z = best_z_med;
  *rank = best_rank_med;

  return (best_rank_med == 0);

}






//
// 
//
int  testMotifPairPosition(float* M1_cnt, float* M2_cnt, int*** dist, int a, int b, 
			   int n,
			   int mbins_dist, int shuffle, int* rank, double* z, float* p)
{
  int    i;
  float  med1, med2;
  float* MM1;
  int*   MM1_q;
  float* MM2;
  int*   MM2_q;
  float* MM1_q_bins;
  float* MM2_q_bins;

  int    nbgenes;
  double mi_med;
  int    rank_med;
  int    best_rank_med;
  double z_med;
  double best_z_med;
  float  pe;

  //
  // determine the median distance for the two motifs
  //

  MM1   = (float*)malloc(n * sizeof(float));
  MM2   = (float*)malloc(n * sizeof(float));

  nbgenes = 0;
  for (i=0; i<n; i++) {
    if ((M1_cnt[i] > 0) && (M2_cnt[i] > 0)) {
      med1       = average_int(dist[i][a], (int)(M1_cnt[i]));
      med2       = average_int(dist[i][b], (int)(M2_cnt[i]));
      MM1[nbgenes] = med1;
      MM2[nbgenes] = med2;
      nbgenes  ++;
    }
  }
  
  pe       = pearson(MM1, MM2, nbgenes);
  *p = pe;
  
  add_small_values_to_identical_floats(MM1, nbgenes);	  
  add_small_values_to_identical_floats(MM2, nbgenes);	  
  
  //printf("quantizing using %d bins.\n", mbins_dist);

  best_rank_med = shuffle + 100;
  best_z_med    = -100;
  for (mbins_dist=2; mbins_dist<=5; mbins_dist++) {

    MM1_q   = Quantize(MM1, nbgenes, mbins_dist, &MM1_q_bins); 
    MM2_q   = Quantize(MM2, nbgenes, mbins_dist, &MM2_q_bins); 
      
    // calc MI
    mi_med  = CalculateMIbasic(MM1_q, MM2_q, nbgenes, mbins_dist, mbins_dist);       
    max_rank_test (mi_med, MM1_q, mbins_dist, MM2_q, mbins_dist, nbgenes, shuffle, 1, &rank_med, 1, &z_med); 
  
    if ((rank_med < best_rank_med) || ((rank_med == best_rank_med) && (rank_med == 0) && (z_med > best_z_med))) {
      best_z_med    = z_med;
      best_rank_med = rank_med;      
    }

  }


  *z = best_z_med;
  *rank = best_rank_med;

  return (best_rank_med == 0);

}



//
// add LEN !
//
int  testMotifClustering(float* M1_cnt, float* M2_cnt, int*** dist, int a, int b, int* gene_length, int n, float* d_avg)
{
  int g, h, i, j, k, cnt;
  double avg_d, d, my_avg_d;
  int  cnt_smaller;
  int* kmer_pos;
  int* kmer_pos1;
  int* kmer_pos2;

  //
  // determine the median distance between the median position
  //
  avg_d = 0.0;
  cnt   = 0;
  for (i=0; i<n; i++) {
    
    if ((M1_cnt[i] > 0) && (M2_cnt[i] > 0)) {
      // coccur here
      
      /*
      med1   = median_int(dist[i][a], (int)(M1_cnt[i]));
      med2   = median_int(dist[i][b], (int)(M2_cnt[i]));
      d      = (float)(fabs(med1 - med2) );
      avg_d += d;
      cnt   ++;
      */
      
      for (g=0; g<(int)M1_cnt[i]; g++) {
	for (h=0; h<(int)M2_cnt[i]; h++) {
	  d      = fabs(dist[i][a][g] - dist[i][b][h]);
	  avg_d += d;
	  cnt ++;
	}
      }
      
    }
    
  }

  avg_d /= cnt;

  *d_avg = avg_d;

  //
  // now test whether avg_d is smaller than expected
  //

  kmer_pos  = (int*)malloc(10000 * sizeof(int));
  kmer_pos1 = (int*)malloc(10000 * sizeof(int));
  kmer_pos2 = (int*)malloc(10000 * sizeof(int));

  cnt_smaller = 0;
  for (j=0; j<10000; j++) {    
    my_avg_d = 0.0;
    cnt   = 0;

    for (i=0; i<n; i++) {
      if ((M1_cnt[i] > 0) && (M2_cnt[i] > 0)) {
	
	
	for (k=0; k<(int)(M1_cnt[i]); k++) {
	  kmer_pos1[k] = (int) (0.5 + gene_length[i] * default_rand()); 
	}
	for (k=0; k<(int)(M2_cnt[i]); k++) {
	  kmer_pos2[k] = (int) (0.5 + gene_length[i] * default_rand()); 
	}

	// coccur here	
	for (g=0; g<(int)M1_cnt[i]; g++) {	  
	  for (h=0; h<(int)M2_cnt[i]; h++) {
	    d         = fabs(kmer_pos1[g] - kmer_pos2[h]);
	    my_avg_d += d;
	    cnt ++;
	  }
	}
	
	/*
	for (k=0; k<(int)(M1_cnt[i]); k++) {
	  kmer_pos[k] = (int) (0.5 + gene_length[i] * default_rand()); 
	}
	med1      = median_int(kmer_pos, (int)(M1_cnt[i]));
	for (k=0; k<(int)(M2_cnt[i]); k++) {
	  kmer_pos[k] = (int) (0.5 + gene_length[i] * default_rand()); 
	}
	med2      = median_int(kmer_pos, (int)(M2_cnt[i]));
	d         = (float)(fabs(med1 - med2) );
	my_avg_d += d;
	cnt ++;
	*/
      }
    }

    my_avg_d /= cnt;

    //printf("avg_d = %f, my_avg_d = %f\n", avg_d, my_avg_d);
    if (my_avg_d <= avg_d) {
      cnt_smaller ++;
    }
    
  }

  //printf("Res: cnt_smaller = %d\n", cnt_smaller);
  
  return cnt_smaller;
  
}

     

void readKmers_this (char* kmerfile, int max_nb_kmers, int inc_nb_kmers, char*** kmers, int* nbkmers, int* kmersize, int** dna_rna) 
{
  char* buff;
  FILE* fp;
  char* s;
  int   len = 10000;

  buff = (char*)malloc(1000 * sizeof(char));
    
  //
  // START READING THE LIST OF KMERS
  //
  *kmers  = (char**)malloc(max_nb_kmers * sizeof(char*));

  *dna_rna  = (int*)malloc(max_nb_kmers * sizeof(int));
  
  
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
    (*kmers)[*nbkmers] = (char*)calloc(100, sizeof(char));   
    strcat((*kmers)[*nbkmers], s);    
    free(s);

    s = mystrtok(0, '\t');
    (*dna_rna)[ *nbkmers ] = atoi(s);

    (*nbkmers)++;

    // allocate more memory if we need it
    if (*nbkmers == max_nb_kmers) {
      die("Olivier: allocate more space for the kmers\n");
      *kmers = allocate_more_kmers(*kmers, *nbkmers, max_nb_kmers, inc_nb_kmers);
      max_nb_kmers += inc_nb_kmers;
    }
  }
  
  // motif size becomes the length of the first k-mer in the list
  if (*nbkmers > 0) {
    *kmersize = strlen((*kmers)[0]);    
  } 
  fclose(fp);     
  //
  // END READING THE LIST OF KMERS
  //
  
}     
