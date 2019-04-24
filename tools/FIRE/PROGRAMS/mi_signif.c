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


typedef struct _Kmer {

  int    index;      // where is it found
  double score;
  float  red;
  int    robustness;
  double zscore;
  int    rank;

} Kmer;

int CmpFunc(const void* _a, const void* _b);


int evalSeed(int ii, short** gene_kmer, int nborfs, double score, int mbins, int* E_q, int ebins, int shuffle, int jn, int jn_t); 


int main(int argc, char ** argv)
{
     

  // GENERAL
  int i;
  char** kmers;
  char** seeds;
  int kmersize;
  int nbkmers;
  char* kmerfile = 0;

  int nborfs = 0;

  FILE  *fp;
  int    size;
  char*  seq;
  char* name;

  int   nextSequence_started = 0;
  int   nextSequence_ended   = 0;
  char* nextSequence_currentLine=0; 
  char* fastafile;

  ENTRY e;
  ENTRY *ep;
  char* realname;


  char*   expfile;
  int     m, n;
  char**  rownames;
  char**  colnames;
  float** data;
   
  int     kmer_np;
  int     kmer_no;
  int     true_idx = 0;
  int     idx;
  char*   outprmfile;

  short** gene_kmer;
  short** gene_seed;

  int*    idx_gene_name;
  int*    idx_gene_name_rank;   

  float*  M;
  int*    M_q;

  float*  M_motif;
  int*    M_motif_q;

  float*  M_seed;
  int*    M_seed_q;

  float*  E;
  int*    E_q;
  float*  E_q_bins;


  int     verbose = 1;

  int*     kmer_counts;

  int*     kmer_pos;
  int*     kmer_ori;
   
  int         singlestrand = 0;

  int         gap = 0;
   
  


  //double    mymi;

  int      quantized = 0;
  float    F,G;
  int      ebins   = 0;
  int      mbins   = 2;

  int      shuffle = 10000;

  extern   int CmpInt();
  int      max_nb_kmers = 200000;
  int      inc_nb_kmers = 10000;
  double*   a_mi_shu;
  //double    z;

  int      seed = 12345;  // random
  char*    outfile = 0;


  int      report  = 1;
   
  FILE*    fpr_motifs = 0;
  FILE*    fpr_seeds  = 0;

  char*    outrepfile;

  float    divbins     = 50.0;
   
  int      count_seq;
  char*    count_seq_index;

  int      jn   = 5;
  int      jn_t = 3;
  int      jn_f = 3;

  char     stmp[400];
  int      pass1;
  int      pass2;
  //int      rank;
  //int      ro;
  //int      sortbymi = 0;
  int      optimout = 0;

  int      outprm   = 0;
  
  int      the_signif_cat = -1;
  char**   signif_cats;
  
  int      doreportonly = 0;
  
  int      motif_ro;
  double   motif_z;
  int      motif_rank;
  float    motif_red;
  double   motif_mi;

  int      seed_ro;
  double   seed_z;
  int      seed_rank;
  float    seed_red = 0.0;
  double   seed_mi = 0.0;

  signif_cats = (char**)malloc(3 * sizeof(char*));
  signif_cats[0] = "OK";
  signif_cats[1] = "OK-NO-SEED";
  signif_cats[2] = "NOT-SIGNIF";
  
  if (argc == 1) {
    printf("Usage : mi_signif -expfile FILE -motiffile FILE -fastafile FILE -rna INT\n");
    exit(0);
  }

  expfile         = get_parameter(argc, argv, "-expfile");

  if (exist_parameter(argc, argv, "-motiffile")) 
    kmerfile        = get_parameter(argc, argv, "-motiffile");

  if (exist_parameter(argc, argv, "-report")) 
    report = atoi(get_parameter(argc, argv, "-report"));

  if (exist_parameter(argc, argv, "-optimout")) 
    optimout = atoi(get_parameter(argc, argv, "-optimout"));

  if (exist_parameter(argc, argv, "-divbins")) 
    divbins = atof(get_parameter(argc, argv, "-divbins"));
  
  if (exist_parameter(argc, argv, "-doreportonly")) 
    doreportonly = atoi(get_parameter(argc, argv, "-doreportonly"));
  
  if (exist_parameter(argc, argv, "-jn")) 
    jn = atoi(get_parameter(argc, argv, "-jn"));

  if (exist_parameter(argc, argv, "-jn_f")) 
    jn_f = atoi(get_parameter(argc, argv, "-jn_f"));

  if (exist_parameter(argc, argv, "-jn_t")) 
    jn_t = atoi(get_parameter(argc, argv, "-jn_t"));
 
  if (exist_parameter(argc, argv, "-outfile")) 
    outfile        = get_parameter(argc, argv, "-outfile");

  fastafile       = get_parameter(argc, argv, "-fastafile");

  if (exist_parameter(argc, argv, "-seed")) 
    seed = atoi(get_parameter(argc, argv, "-seed"));
  default_set_seed(seed);
   
  if (exist_parameter(argc, argv, "-quantized")) 
    quantized = atoi(get_parameter(argc, argv, "-quantized"));

  if (exist_parameter(argc, argv, "-shuffle")) 
    shuffle = atoi(get_parameter(argc, argv, "-shuffle"));
  if (exist_parameter(argc, argv, "-ebins")) 
    ebins            = atoi(get_parameter(argc, argv, "-ebins"));

  if (exist_parameter(argc, argv, "-mbins")) {
    mbins            = atoi(get_parameter(argc, argv, "-mbins"));
  }

  if (exist_parameter(argc, argv, "-singlestrand")) 
    singlestrand    = atoi(get_parameter(argc, argv, "-singlestrand"));   
   
  if (exist_parameter(argc, argv, "-rna")) 
    singlestrand    = atoi(get_parameter(argc, argv, "-rna"));   
   
  if (exist_parameter(argc, argv, "-verbose")) 
    verbose         = atoi(get_parameter(argc, argv, "-verbose"));   

  

  if (optimout == 1) {
    readKmers_general_special_optim (kmerfile, max_nb_kmers, inc_nb_kmers, &kmers, &seeds, &nbkmers, &kmersize); 
  } else {
    readKmers_general               (kmerfile, max_nb_kmers, inc_nb_kmers, &kmers, &nbkmers, &kmersize); 
  }
  
  if (nbkmers == 0) {
    die("Please a non-empty motif file.\n");
  }

  kmer_counts = (int*)calloc(nbkmers, sizeof(int));
   

  //
  //  read in expression data
  //
  readFloatTable(expfile, &m, &n, &data, &rownames, &colnames, 0, 1);

  //m = 1;

  E             = (float* ) malloc( n * sizeof(float));
  gene_kmer     = (short**) malloc( n * sizeof(short*));
  gene_seed     = (short**) malloc( n * sizeof(short*));

  kmer_pos      = (int*)    malloc( 10000 * sizeof(int));
  kmer_ori      = (int*)    malloc( 10000 * sizeof(int));

  idx_gene_name      = (int*)malloc( n * sizeof(int));
  idx_gene_name_rank = (int*)malloc( n * sizeof(int));

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
    gene_kmer    [ true_idx ] = (short*)calloc(nbkmers, sizeof(short));
       
    for (i=0; i<nbkmers; i++) {
      
      kmer_np = 0;
      kmer_no = 0;
      if (mbins > 2) 
	findSites(kmers[i], seq, kmer_pos, &kmer_np, kmer_ori, &kmer_no, 10000, singlestrand, 1);
      else
	kmer_np =  re_matches(kmers[i], seq, singlestrand); 
	  
      gene_kmer[ true_idx ][ i ] = (short)kmer_np;    

    }

    
    if (optimout == 1) {
      
      gene_seed    [ true_idx ] = (short*)calloc(nbkmers, sizeof(short));

      for (i=0; i<nbkmers; i++) {
	
	kmer_np = 0;
	kmer_no = 0;
	
	if (mbins > 2) 
	  findSites(seeds[i], seq, kmer_pos, &kmer_np, kmer_ori, &kmer_no, 10000, singlestrand, 1);
	else
	  kmer_np =  re_matches(seeds[i], seq, singlestrand); 

	gene_seed[ true_idx ][ i ] = (short)kmer_np;    

      }
      
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
   
  
  
  if (report == 1) {
    if (outfile != 0) {
      
      outrepfile = (char*)calloc(1000, sizeof(char));
      strcat(outrepfile, outfile);
      strcat(outrepfile, ".motifs.rep");
      fpr_motifs = fopen( outrepfile, "w");
      free(outrepfile);
      
      if (optimout == 1) {
	outrepfile = (char*)calloc(1000, sizeof(char));
	strcat(outrepfile, outfile);
	strcat(outrepfile, ".seeds.rep");
	fpr_seeds = fopen( outrepfile, "w");
	free(outrepfile);
      }

    } else {
      fpr_motifs = fopen( "report.motifs.txt", "w");
      if (optimout == 1) {
	fpr_seeds = fopen( "report.seeds.txt", "w");
      }      
    }
  }



  //
  //  PROCESS EXPRESSION DATA
  //
  quantize_E(E, nborfs, quantized, &ebins, &E_q, &E_q_bins); 
   

  if (shuffle > 0) {
    a_mi_shu = (double*)malloc(sizeof(double) * shuffle );  
  }

  //						
  //  PROCESS ALL MOTIFS[/SEED PAIRS]
  //
  for (i=0; i<nbkmers; i++) {
    
    //
    // START: MOTIF
    //
    M = stof_matrix_column(gene_kmer, i, nborfs);     
    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

    motif_mi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins);

    if (quantized == 0) {
      fitSimpleLinearModel(M, E, nborfs, &F, &G);
      motif_red   = F;
    } else {
      motif_red   = -1;
    }
    
    motif_rank = -1;
    motif_ro   = -1;
    motif_z    = -1.0;
    
    if (doreportonly == 0) {
      pass1 = max_rank_test   (motif_mi, M_q, mbins, E_q, ebins, nborfs, shuffle, 1, &motif_rank, 1, &motif_z); 
      pass2 = jacknife_max_rank_test(M_q, mbins, E_q, ebins, nborfs, shuffle, jn, jn_f, jn_t, 1, &motif_ro);       
    }
    
    free(M);
    free(M_q);
    
    //
    // END: MOTIF
    //

    //
    // START: SEED
    //
    
    if (optimout == 1) {
      
      
      M = stof_matrix_column(gene_seed, i, nborfs);     
      quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

      seed_mi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins);
    
      if (quantized == 0) {
	fitSimpleLinearModel(M, E, nborfs, &F, &G);
	seed_red   = F;
      } else {
	seed_red   = -1;
      }
      
      seed_rank = -1;
      seed_z    = -1.0;
      seed_ro   = -1;

      if (doreportonly == 0) {
	pass1 = max_rank_test   (seed_mi, M_q, mbins, E_q, ebins, nborfs, shuffle, 1, &seed_rank, 1, &seed_z); 
	pass2 = jacknife_max_rank_test(M_q, mbins, E_q, ebins, nborfs, shuffle, jn, jn_f, jn_t, 1, &seed_ro); 
      }
      
      free(M);
      free(M_q);
      
    }
    
    //
    // END: SEED
    //

 
    if ((motif_rank == 0) && (motif_ro >= jn_t)) {
      the_signif_cat  = 0;  // ok
    } else {
      the_signif_cat  = 2;  // non-signif
    }

    //
    // may update down the category
    //
    if (optimout == 1) {
      if ((seed_rank > 0) || (seed_ro < jn_t)) {
	if (the_signif_cat == 0) 
	  the_signif_cat = 1;    // ok-no-seed
      }

    }

    
    //
    //  STANDARD OUTPUT					       
    //

    printf("Processing motif %s.\n", kmers[i]);
    printf("MI=%4.3f\n", motif_mi);
    if (quantized == 0) 
      printf("Lin reg coef=%4.3f\n", motif_red);	  
    printf("Shuffle rank=%d\n", motif_rank);
    printf("Z-score=%4.3f\n", motif_z);
    printf("Sigificance category= %s\n", signif_cats[ the_signif_cat ]);
    printf("Robustness (jn_f=%d)=%d/%d\n", jn_f, motif_ro, jn);    

           
    sprintf(stmp, "%5.4f\t%d\t%4.3f\t%d\t%s",  motif_mi, motif_rank, motif_z, motif_ro, signif_cats[ the_signif_cat ]); 
    if (outfile != 0) {
      fprintf(fp, "%s\t%s", kmers[i], stmp);      
      if (quantized == 0) {
      	fprintf(fp, "\t%4.3f", motif_red);
      }
    }
    
    if (optimout == 1) {

      sprintf(stmp, "%s\t%5.4f\t%d\t%4.3f\t%d", seeds[i], seed_mi, seed_rank, seed_z, seed_ro); 

      if (outfile != 0) {
	fprintf(fp, "\t%s", stmp);
	if (quantized == 0)
	  fprintf(fp, "\t%4.3f", seed_red);
      }
      
    }
    
    if (outfile != 0) 
      fprintf(fp, "\n");
    
    printf("\n");
    
    //
    // REPORT OUTPUT
    //

    if (report == 1) {
      
      // MOTIFS
      M_motif = stof_matrix_column(gene_kmer, i, nborfs);     
      quantize_M_zero_eqpop(M_motif, nborfs, mbins, &M_motif_q); 

      add_to_report(fpr_motifs, kmers[i], kmersize, gap, M_motif_q, mbins, E_q, ebins, (quantized==1?0:E_q_bins), nborfs);
      
      // SEEDS
      if (optimout == 1) {
	
	M_seed = stof_matrix_column(gene_seed, i, nborfs);     
	quantize_M_zero_eqpop(M_seed, nborfs, mbins, &M_seed_q); 
	
	add_to_report(fpr_seeds, seeds[i], kmersize, gap, M_seed_q, mbins, E_q, ebins, (quantized==1?0:E_q_bins), nborfs);
	
	free(M_seed);
	free(M_seed_q);
	
      }
      
      free(M_motif);
      free(M_motif_q);
      
    }
  
  }

  fclose(fp);
  fclose(fpr_motifs);

  if (optimout == 1) {
    fclose(fpr_seeds);
  }

  if (outfile != 0) {
    
    if (outprm == 1) {
      outprmfile = (char*)calloc(1000, sizeof(char));
      strcat(outprmfile, outfile);
      strcat(outprmfile, ".prm");
      fp = fopen(outprmfile, "w"); 
      fprintf(fp, "expfile\t%s\n" , expfile  );
      fprintf(fp, "kmerfile\t%s\n" , kmerfile  );
      fprintf(fp, "quantized\t%d\n" , quantized  );
      fprintf(fp, "shuffle\t%d\n" , shuffle  );
      fprintf(fp, "gap\t%d\n" , gap  );
      fprintf(fp, "rna\t%d\n" , singlestrand  );
      fprintf(fp, "fastafile\t%s\n" , fastafile  );
      fprintf(fp, "ebins\t%d\n" , ebins  );
      fprintf(fp, "mbins\t%d\n" , mbins  );
      fprintf(fp, "sample_size\t%d\n", nborfs);
      fclose(fp);
    }
  }
  

  free(E_q);
  free(E);
  

  return 0;
}


int CmpFunc(const void* _a, const void* _b)
{
  const Kmer* a = (const Kmer*) _a;
  const Kmer* b = (const Kmer*) _b;
  
  if (a->score < b->score) 
    return 1; 
  else if(a->score == b->score) 
    return  0;
  else         
    return -1;
}



