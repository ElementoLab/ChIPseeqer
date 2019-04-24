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

#define DIVBINS_DIST 10.0               // for position bias 
//#define DIVBINS_ELSE 50.0               // for orientation, copy number


int CmpFunc(const void* _a, const void* _b);
void add_to_dist_report(FILE* fpr, char* kmer, int* M_q, int* E_q, int nborfs, float* M_q_bins, float* E_q_bins, int mbins, int ebins, int maxmbins, int quantized, char* label_cor, char** labels);


int main(int argc, char ** argv)
{
     

  // GENERAL
  int i, j, k;
  FILE  *fp;
  ENTRY e;
  ENTRY *ep;
  
  // general parameters
  int     verbose = 0;
  Params  p;
  int     shuffle;

  // kmer stufff
  char**   kmers;

  int      kmersize;
  int      nbkmers;
  char*    kmerfile;
  int      max_nb_kmers = 10000;
  int      inc_nb_kmers = 10000;
   
  // sequence stuff
  int nborfs = 0;
  int         nborfs_nom;

    
  char* realname;
  int    size;
  char*  seq;
  char* name;
  int   nextSequence_started = 0;
  int   nextSequence_ended   = 0;
  char* nextSequence_currentLine=0; 
  char* fastafile;
  int     true_idx = 0;
  int     idx;
  int         singlestrand = 0;
  int         l;
  
  // motif stuff
  int*        kmer_pos;
  int*        kmer_ori;
  float**     gene_kmer_avgpos;
  float**     gene_kmer_medpos;
  float**     gene_kmer_minpos;
  short**     gene_kmer_sinori;
  short**     gene_kmer_ori3;
  short**     gene_kmer_ori5;
  short**     gene_kmer_cnt;


  int*        Mavg_q;

  int*        Mcnt_q;

  int         no, np;

  float*      Mavg_nom;
  double      mi_avg;
  double      mi_k;
  double      mi_o5, mi_o3, z_o5, z_o3;
  double      z_avg;
  double      z_k;
  float*      Mavg;
  float*      Mcnt;
  int         mindist;
  int         ori1;
  int         ori2;
  float*      Mavg_q_bins;
  int         mbins;
  float*      O5;
  float*      O3;
  int*        O5_q;
  int*        O3_q;
  float*      E;
  float*      E_q_bins;
  int*        E_q;
  int         ebins;
  float*      E_q_bins_dist;
  int*        E_q_dist;
  int         ebins_dist; 


  //
  // expression stuff
  //
  char*     expfile   = 0;
  char*     expfile_b = 0;

  int       m, n;
  char**    rownames;
  char**    colnames;
  float**   data;

  int       m_ori, n_ori;
  char**    rownames_ori;
  char**    colnames_ori;
  float**   data_ori;
  float*    E_ori;
  
 
  int       quantized = 0;
 
  float*    E_nom;


  int     shuffle_rank;

  char*   outfile;
  char*   oldoutfile = 0;
  int     report;
  char*   outrepfile;
  FILE*   fpr = 0, *ofp = 0;

  int     outprm = 0;
  char*   outprmfile;
  int     seed = 12345;  // random
 

  
  int     me;
  extern   int CmpInt();

  int     rank_avg;

  int     rank_o5, rank_o3, rank_k;
  char**  labels;

  int     avg_len            = 0;
  int     mbins_interval     = 100;
  int     mbins_dist;
  int     nbcopies_shuffle   = 10000;
  int     donbcopies         = 1;
  int     best_rank_avg;
  double  best_z_avg         = -10000000.0;
  int*    best_Mavg_q        = 0;
  double  best_mi_avg        = -10000000.0;
  int     best_mbins_dist    = -10000000;
  int*    best_E_q_dist      = 0; 
  int     best_ebins_dist    = -10000000;
  float*  best_Mavg_q_bins   = 0;
  float*  best_E_q_bins_dist = 0;
  
  float** Edata_tmp;
  float** Edata;
  float    divbins           = 50.0;
   
   if (argc == 1) {
     printf("Usage : mi_dist -expfile FILE -expfile_b FILE -kmerfile FILE -fastafile FILE -k INT -rna [01] -newoutfile FILE -mbins_dist \n");
     exit(0);
   }

   
   p.quantized    = 0;
   p.shuffle      = 10000;
   p.shuffle_rank = 0;
   p.singlestrand = 0;
   p.verbose      = 0;
   p.mbins        = 2;
   p.outfile      = 0;
   p.report       = 1;
   p.mbins_dist   = 0;



   get_Params(argc, argv, &p); 

   kmerfile     = p.kmerfile;
   fastafile    = p.fastafile;
   kmersize     = p.kmersize;
   expfile      = p.expfile;
   singlestrand = p.singlestrand;
   quantized    = p.quantized;
   shuffle      = p.shuffle;
   shuffle_rank = p.shuffle_rank;
   verbose      = p.verbose;
   mbins        = p.mbins;
   mbins_dist   = p.mbins_dist;
   report       = p.report;
   outfile      = p.outfile;
      
   if (exist_parameter(argc, argv, "-expfile_b")) 
     expfile_b        = get_parameter(argc, argv, "-expfile_b");
   
   if (exist_parameter(argc, argv, "-divbins")) 
     divbins          = atof(get_parameter(argc, argv, "-divbins"));
   
   if (exist_parameter(argc, argv, "-nbcopies_shuffle")) 
     nbcopies_shuffle = atoi(get_parameter(argc, argv, "-nbcopies_shuffle"));

   if (exist_parameter(argc, argv, "-seed")) 
     seed = atoi(get_parameter(argc, argv, "-seed"));
   default_set_seed(seed);
   
   if (exist_parameter(argc, argv, "-donbcopies")) 
     donbcopies       = atoi(get_parameter(argc, argv, "-donbcopies"));

   if (exist_parameter(argc, argv, "-outprm")) 
     outprm           = atoi(get_parameter(argc, argv, "-outprm"));

   if (exist_parameter(argc, argv, "-mbins_interval")) 
     mbins_interval = atoi(get_parameter(argc, argv, "-mbins_interval"));
   
   
   if (exist_parameter(argc, argv, "-oldoutfile")) 
     oldoutfile    = get_parameter(argc, argv, "-oldoutfile");

   if (shuffle_rank == 0) {
     shuffle_rank = shuffle;  // default
   }

   readKmers_general (kmerfile, max_nb_kmers, inc_nb_kmers, &kmers, &nbkmers, &kmersize); 

   if (nbkmers == 0) {
     printf("no motif in input file, exiting.\n");
     exit(0);
   } else {
     printf("Read %d motifs.\n", nbkmers);
   }

   //
   //  read in expression data
   //
   
   readFloatTable(expfile_b, &m, &n, &data, &rownames, &colnames, 0, 1);
   printf("Read %d (genes) x %d (conditions) matrix.\n", n, m);
   
   
   //
   //  read in original data
   //
   readFloatTable(expfile  , &m_ori, &n_ori, &data_ori, &rownames_ori, &colnames_ori, 0, 1);

   printf("Read %d (genes) x %d (conditions) matrix.\n", n_ori, m_ori);
   
   E_ori            = (float* ) malloc( n * sizeof(float));      
   gene_kmer_avgpos = (float**) malloc( n * sizeof(float*));
   Edata_tmp        = (float**) malloc( n * sizeof(float*));
   gene_kmer_medpos = (float**) malloc( n * sizeof(float*));
   gene_kmer_minpos = (float**) malloc( n * sizeof(float*));
   gene_kmer_sinori = (short**) malloc( n * sizeof(short*));
   gene_kmer_ori3   = (short**) malloc( n * sizeof(short*));
   gene_kmer_ori5   = (short**) malloc( n * sizeof(short*));
   gene_kmer_cnt    = (short**) malloc( n * sizeof(short*));
   
   kmer_pos      = (int*)    malloc( 10000 * sizeof(int));
   kmer_ori      = (int*)    malloc( 10000 * sizeof(int));

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

   printf("Reading sequences and building motif profiles ...");

   nextSequence_currentLine = (char*)malloc(200000 * sizeof(char));
   nborfs  = 0;
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
       continue;
     } 

     // get the index for the gene
     idx = (int)ep->data;       

     // update expression vector
     E_ori[ true_idx ]         = data_ori[idx][0];
     Edata_tmp[ true_idx ]     = data[idx];
     
     
     // allocate mem for position profile     
     gene_kmer_avgpos[ true_idx ] = (float*)calloc(nbkmers, sizeof(float));
     gene_kmer_medpos[ true_idx ] = (float*)calloc(nbkmers, sizeof(float));
     gene_kmer_minpos[ true_idx ] = (float*)calloc(nbkmers, sizeof(float));

     gene_kmer_sinori[ true_idx ] = (short*)calloc(nbkmers, sizeof(short));

     gene_kmer_ori5  [ true_idx ] = (short*)calloc(nbkmers, sizeof(short));
     gene_kmer_ori3  [ true_idx ] = (short*)calloc(nbkmers, sizeof(short));

     gene_kmer_cnt   [ true_idx ] = (short*)calloc(nbkmers, sizeof(short));

     
     if (verbose == 2) printf("%s\n", realname);
     
     l  = strlen(seq);

     avg_len += l;

     // create position profile
     for (i=0; i<nbkmers; i++) {
       
       no = 0; np = 0;
       
       kmer_pos = (int*)malloc(strlen(seq) * sizeof(int));
       kmer_ori = (int*)malloc(strlen(seq) * sizeof(int));

     
       findSites(kmers[i], seq, kmer_pos, &np, kmer_ori, &no, l, singlestrand, 1) ;

       gene_kmer_cnt[ true_idx ][ i ] = np;

       if (np > 0) {

	 
	 //
	 // calculate median position
	 //
	 me = median_int(kmer_pos, np);
	 if (singlestrand == 0) 
	   gene_kmer_medpos[ true_idx ][ i ] = (l - me);
	 else 
	   gene_kmer_medpos[ true_idx ][ i ] = me;
	 
	 //
	 // calculate min position
	 //
	 mindist = INT_MAX;	
	 
	 for (j=0; j<np; j++) {

	   if (singlestrand == 0) {
	     if ((l-kmer_pos[j]) < mindist) {
	       mindist = l - kmer_pos[j];
	     }
	     
	   } else {
	     if (kmer_pos[j] < mindist) {
	       mindist = kmer_pos[j];
	     }

	   }
	 }
	 
	 gene_kmer_minpos[ true_idx ][ i ] = (float)mindist; 
	 
	 //
	 // calculate average position
	 //
	 gene_kmer_avgpos[ true_idx ][ i ] = 0;
	 for (j=0; j<np; j++) {
	   if (singlestrand == 0) 
	     gene_kmer_avgpos[ true_idx ][ i ] += (l - kmer_pos[ j ]);
	   else 
	     gene_kmer_avgpos[ true_idx ][ i ] += kmer_pos[ j ];
	 }
	 gene_kmer_avgpos[ true_idx ][ i ] = gene_kmer_avgpos[ true_idx ][ i ]  / np;

       } else {
	 gene_kmer_minpos[true_idx][i] = NAN;
	 gene_kmer_avgpos[true_idx][i] = NAN;
	 gene_kmer_medpos[true_idx][i] = NAN;
       }

       //
       // orientation bias					
       //
            
       if (verbose == 2) printf("%d:BEF: got %d ori\n", i, no);
       
       if (singlestrand == 1) {
	 // different : must recalculate, but in double strand mode
	 //kmer_pos = (int*)malloc(strlen(seq) * sizeof(int));
	 //kmer_ori = (int*)malloc(strlen(seq) * sizeof(int));

	 //no = 0; np = 0;
	 //printf("relook for %s\n", kmers[i]);
	 findSites(kmers[i], seq, kmer_pos, &np, kmer_ori, &no, l, 0, 1) ;
	 
       } 
	 
       if (verbose == 2) printf("%d:AFT: got %d ori\n\n", i, no);
            
       if (no > 0) {
	 
	 ori1 = 0;
	 ori2 = 0;
	 
	 for (j=0; j<no; j++) {
	   if (kmer_ori[ j ] == 1) {
	     gene_kmer_ori5[ true_idx ][ i ] ++;
	     ori1 ++;
	   } else {
	     gene_kmer_ori3[ true_idx ][ i ] ++;
	     ori2 ++;
	   }
	 }
	 
	 if (ori1 == ori2)
	   gene_kmer_sinori[ true_idx ][ i ] = 1;
	 else if (ori1 > ori2)
	   gene_kmer_sinori[ true_idx ][ i ] = 0;
	 else
	   gene_kmer_sinori[ true_idx ][ i ] = 2;
	 
	 
	 
       } else {
	 gene_kmer_sinori[true_idx][i] = SHRT_MAX;
	 gene_kmer_ori5  [true_idx][i] = 0; //SHRT_MAX;
	 gene_kmer_ori3  [true_idx][i] = 0; //SHRT_MAX;
       }
       
       
       
       free(kmer_pos);
       free(kmer_ori);     
       
     }
     
     true_idx ++;
     
     free(seq);
     free(name);

   }

   if (verbose == 1) 
     printf("read and processed sequences, avg and min pos for each kmer.\n");

   printf("Done.\n");




   nborfs = true_idx;

   //
   transpose_f(Edata_tmp, nborfs, m, &Edata);

   avg_len = (int)( 0.5 + avg_len / (float)nborfs );
   printf("Average sequence length is %d.\n", avg_len);

   fclose(fp);

   //
   //  open outfile and report file
   //
   if (outfile != 0) {
     fp = fopen(outfile, "w");
     if (fp == 0) {
       printf("cannot open outfile: %s\n", outfile);
     }
   }

   if (oldoutfile != 0) {
     ofp = fopen(oldoutfile, "w");
     if (ofp == 0) {
       printf("cannot open oldoutfile: %s\n", oldoutfile);
     }
   }

   if (report == 1) {
     if (outfile != 0) {
       outrepfile = (char*)calloc(1000, sizeof(char));
       strcat(outrepfile, outfile);
       strcat(outrepfile, ".rep");
       fpr = fopen( outrepfile, "w");
     } else {
       fpr = fopen( "report.txt", "w");
     }
   }

   
   if (quantized == 0) {
     // add a little random number to each value in the E vector
     add_small_values_to_identical_floats(E_ori, nborfs);
  }

   //
   //  START LOOKING AT MOTIFS
   //
   
   for (i=0; i<nbkmers; i++) {

     //
     // get the ith column of the gene_kmer matrix
     //
     Mavg       = f_matrix_column   (gene_kmer_avgpos, i, nborfs);
     O5         = stof_matrix_column(gene_kmer_ori5,   i, nborfs);
     O3         = stof_matrix_column(gene_kmer_ori3,   i, nborfs);
     Mcnt       = stof_matrix_column(gene_kmer_cnt,    i, nborfs);

     //
     // E is specific to that motif here
     //
     E          = Edata[i];

     //									
     // remove missing values 1) allocate memory
     //
     Mavg_nom   = (float*)malloc(nborfs * sizeof(float));
     E_nom      = (float*)malloc(nborfs * sizeof(float));
     nborfs_nom = 0;
	    
     
     for (j=0; j<nborfs; j++) {
       if ((Mcnt[j] > 0) && (!isnan(E[j]))) {
	 Mavg_nom[nborfs_nom] = Mavg[j];
	 E_nom   [nborfs_nom] = E[j];
	 nborfs_nom++;
       }
     }

     //
     // variable stuff
     //
     add_small_values_to_identical_floats(Mavg_nom, nborfs_nom);	  

     
     if (quantized == 0) {
       ebins      = (int)(0.5 + nborfs / ( mbins * divbins ));
     }

     best_rank_avg = shuffle + 100;

     for (mbins_dist=2; mbins_dist<=5; mbins_dist++) {

       // 1 means 'always quantized'
       quantize_E(E_nom, nborfs_nom, 1, &ebins_dist,  &E_q_dist,     &E_q_bins_dist); 
       
       Mavg_q   = Quantize(Mavg_nom, nborfs_nom, mbins_dist, &Mavg_q_bins); 
       mi_avg   = CalculateMIbasic(Mavg_q, E_q_dist, nborfs_nom, mbins_dist, ebins_dist);   

       max_rank_test (mi_avg, Mavg_q, mbins_dist, E_q_dist, ebins_dist, nborfs_nom, shuffle, 1, &rank_avg, 1, &z_avg); 

       //CalculateAndShowCountMatrix(Mavg_q,  E_q_dist, nborfs_nom, mbins_dist, ebins_dist, 0, 0, 0);

       if ((rank_avg < best_rank_avg) || ( (rank_avg == best_rank_avg) && (rank_avg == 0) && (z_avg > best_z_avg)) ) {
	 best_rank_avg      = rank_avg;
	 best_z_avg         = z_avg;
	 best_Mavg_q        = Mavg_q;
	 best_mi_avg        = mi_avg;
	 best_mbins_dist    = mbins_dist;
	 best_E_q_dist      = E_q_dist;
	 best_E_q_bins_dist = E_q_bins_dist;
	 best_ebins_dist    = ebins_dist;
	 best_Mavg_q_bins   = Mavg_q_bins;
       }
     }


     //
     // orientation stuff					       
     //

     // here use input quantization instructions


     quantize_E(E_ori, nborfs, quantized, &ebins, &E_q, &E_q_bins); 
     printf("Quantized E vector into %d bins\n", ebins);
     //for (j=0; j<nborfs; j++) {
     //  printf("e=%5.4f, eq=%d\n", E_ori[j], E_q[j]);
     //}
     quantize_M_zero_eqpop(O5, nborfs, mbins, &O5_q); 
     quantize_M_zero_eqpop(O3, nborfs, mbins, &O3_q); 

     mi_o5    = CalculateMIbasic(O5_q, E_q, nborfs, mbins, ebins);   
     mi_o3    = CalculateMIbasic(O3_q, E_q, nborfs, mbins, ebins);     

     max_rank_test (mi_o5,  O5_q,   mbins, E_q, ebins, nborfs, shuffle, 1, &rank_o5, 1, &z_o5); 
     max_rank_test (mi_o3,  O3_q,   mbins, E_q, ebins, nborfs, shuffle, 1, &rank_o3, 1, &z_o3); 

     //if ( (mi_min>c_pmin) || (mi_avg>c_pavg) || (mi_o_sin>c_osin)) {
     
     printf("Motif %s, present in %d genes.\n", kmers[i], nborfs_nom);
     printf("       Average distance: MI=%5.4f, rank=%d/%d, Z=%5.4f, mbins=%d\n", best_mi_avg,   best_rank_avg, shuffle, best_z_avg, best_mbins_dist); 
     printf("  Orientation bias {5'): MI=%5.4f, rank=%d/%d, Z=%5.4f\n", mi_o5,    rank_o5 , shuffle, z_o5 ); 
     printf("  Orientation bias (3'): MI=%5.4f, rank=%d/%d, Z=%5.4f\n", mi_o3,    rank_o3 , shuffle, z_o3 ); 

   
     if (outfile != 0) {
       fprintf(fp, "%s\t%s\t%d\t%5.4f\t%d\t%5.4f\n", kmers[i], "d_avg", nborfs_nom,   best_mi_avg,   best_rank_avg, best_z_avg);
       fprintf(fp, "%s\t%s\t%d\t%5.4f\t%d\t%5.4f\n", kmers[i], "o5",    nborfs,       mi_o5,    rank_o5,  z_o5);
       fprintf(fp, "%s\t%s\t%d\t%5.4f\t%d\t%5.4f\n", kmers[i], "o3",    nborfs,       mi_o3,    rank_o3,  z_o3);
     }


     Mcnt_q     = (int*)  malloc(nborfs * sizeof(int));

     if (donbcopies == 1) {
       
       for (k=1; k<=10; k++) {
	 np = 0;
	 
	 // stop if we don't want to do nbcopies
	 if ((donbcopies == 0) && (k > 1)) {
	   break;
	 }
	 
	 for (j=0; j<nborfs; j++) {
	   if (Mcnt[j] >= k) {
	     Mcnt_q[j] = 1;
	     np ++;
	   } else
	     Mcnt_q[j] = 0;
	 }
	 
	 if (np > 0) {
	   
	   mi_k = CalculateMIbasic(Mcnt_q, E_q, nborfs, mbins, ebins);   
	   max_rank_test (mi_k,  Mcnt_q, mbins, E_q, ebins, nborfs, nbcopies_shuffle, 1, &rank_k, 1, &z_k); 
	   printf("       Copy number (%d+): %d genes, MI=%5.4f, rank=%d/%d, Z=%5.4f\n", k, np, mi_k, rank_k, nbcopies_shuffle, z_k); 
	   
	   if (outfile != 0) {
	     fprintf(fp, "%s\tc%d\t%d\t%5.4f\t%d\t%5.4f\n", kmers[i], k, nborfs, mi_k, rank_k, z_k);
	   }
	   
	 } else {
	   
	   break;
	   
	 }
       }
    
     }

     printf("\n");
     


     //
     //  write report file (count matrix)
     //

     if (report == 1) {

       labels = (char**)malloc( mbins * sizeof(char*));

       add_to_dist_report(fpr, kmers[i], best_Mavg_q, best_E_q_dist, nborfs_nom, best_Mavg_q_bins, best_E_q_bins_dist,
			  best_mbins_dist, best_ebins_dist, best_mbins_dist, 
			  1, "d_avg", 0);  // exp vector is always quantized here
       
       labels[0] = "0";
       labels[1] = "1+";
       add_to_dist_report(fpr, kmers[i], O5_q, E_q, nborfs, 0, E_q_bins,
			  mbins, ebins, mbins, 
			  quantized, "o5", labels);

       add_to_dist_report(fpr, kmers[i], O3_q, E_q, nborfs, 0, E_q_bins,
			  mbins, ebins, mbins, 
			  quantized, "o3", labels);

       
        
     }

     free(Mavg);
     //free(Mmed);
     //free(Mmin_q);
     free(Mavg_q);
     //free(Mmed_q);
     free(E_nom);
     //free(Mmin_nom);
     free(Mavg_nom);
     //free(Mmed_nom);
     free(E_q_dist);
     
     
     free(E_q);
     

   

    
     
   
   }
   
   if (outfile != 0) 
     fclose(fp);
   
   if (oldoutfile != 0) 
     fclose(ofp);
   
   if (outfile != 0) {

     if (outprm == 1) {
       
       outprmfile = (char*)calloc(1000, sizeof(char));
       strcat(outprmfile, outfile);
       strcat(outprmfile, ".prm");
       
       fp = fopen(outprmfile, "w"); 
       
       fprintf(fp, "expfile\t%s\n" , expfile  );
       fprintf(fp, "quantized\t%d\n" , quantized  );
       fprintf(fp, "shuffle\t%d\n" , shuffle  );
       fprintf(fp, "rna\t%d\n" , singlestrand  );
       fprintf(fp, "fastafile\t%s\n" , fastafile  );
       fprintf(fp, "seed\t%d\n" , seed  );
       fprintf(fp, "ebins\t%d\n" , ebins  );
       fprintf(fp, "mbins\t%d\n" , mbins  );
       fprintf(fp, "sample_size\t%d\n", nborfs);
       
       fclose(fp);
       
     }
   }
    

   if (report == 1) {
     fclose(fpr);
   }

   return 0;
}


void add_to_dist_report(FILE* fpr, char* kmer, int* M_q, int* E_q, int nborfs, float* M_q_bins, float* E_q_bins, int mbins, int ebins, int maxmbins, int quantized, char* label_cor, char** labels)
{
  
  int** tmpCounts;
  int   l, k;
  
  tmpCounts   = ConstructCountMatrix (M_q, E_q, nborfs, mbins, ebins);       

  //
  // header line : kmer, correl type
  //
  fprintf(fpr, "%s\t%s\tnan", kmer, label_cor);
  
  for (l=0; l<mbins; l++) {
    if (M_q_bins == 0) 
      fprintf(fpr, "\t%s", labels[l]);
    else
      fprintf(fpr, "\t[%3.2f;%3.2f]", M_q_bins[l], M_q_bins[l+1]);
  }

  for (l=mbins; l<maxmbins; l++) {
    fprintf(fpr, "\tnan");
  }       
  
  if (quantized == 0) {
    fprintf(fpr, "\tnan");
  }
  fprintf(fpr, "\n");	 
  //
  // end header line
  //

  //
  // traverse all bins
  //
  for (k=0; k<ebins; k++) {
    fprintf(fpr, "%s\t%s\t%d", kmer, label_cor, k);
    
    for (l=0; l<mbins; l++) {
      fprintf(fpr, "\t%d", tmpCounts[l][k]);
    }
    for (l=mbins; l<maxmbins; l++) {
      fprintf(fpr, "\tnan");
    }	 
    if (quantized == 0) {
      if (E_q_bins == 0)
	fprintf(fpr, "\tnan");
      else
	fprintf(fpr, "\t[%4.3f;%4.3f]", E_q_bins[k], E_q_bins[k+1]);

    }	 
    fprintf(fpr, "\n");
  }
  free(tmpCounts);
}
     
