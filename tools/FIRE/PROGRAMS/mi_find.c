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
#include "sequences.h"

#define MIN_KMER_COUNT 20


typedef struct _Kmer {
  int   index;      // where is it found
  double score;
  float red;
} Kmer;

int CmpFunc(const void* _a, const void* _b);
void getPresPosOriVectors(char* seq1, int motifsize, int gap, PT* pt, int singlestrand, short** kmer_pre, float** kmer_pos, short** kmer_ori, int nbkmers);


int evalSeed(int ii, short** gene_kmer, int nborfs, double score, int mbins, int* E_q, int ebins, int shuffle, int jn, int jn_f, int jn_t); 


int main(int argc, char ** argv)
{
     

  // GENERAL
  int i, j;
  char** kmers;
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
   
  //int     kmer_np;
  //int     kmer_no;
  int     true_idx = 0;
  int     idx;

  short** gene_kmer;

  int*    idx_gene_name;
  int*    idx_gene_name_rank;   

  float*   M;
  float*   E;

  int      verbose = 1;
  char*    kmer;
   

  int*     kmer_counts;
  int      pos;

  int*     kmer_pos;
  int*     kmer_ori;
   
  float**  gene_kmer_pos;
  short**  gene_kmer_ori;
  //int      fast = 1;

   
  int      singlestrand = 0;

  int      gap = 0;
   
  int*     E_q;
  int*     E_q_shu;
  float*   E_q_bins;

  int*     M_q;
  double   mi_shu_avg;
  double   mi_shu_std;
  double   mymi;

  int      quantized = 0;
  float    F,G;
  int      ebins   = 0;
  int      mbins   = 2;
  int      ii;
  Kmer*    kmers_array;
  int      shuffle = -1;
  int      save_kmer_profile = 1;
  extern   int CmpInt();
  int      max_nb_kmers = 200000;
  int      inc_nb_kmers = 10000;
  double*   a_mi_shu = 0;
  double    z = 0.0;
  int      max_selected = -1;
  float    max_nbsdev   = 0.0;

  int      seed = 12345;  // random
  char*    outfile = 0;
  char*    outprmfile;
  int      outprm  = 0;
  int      report  = 1;
   
  FILE*    fpr = 0;
  char*    outrepfile;


  float*     seq_lens = 0;
  float*     seq_gc   = 0;
  float*     seq_cpg  = 0;
  int        docor   = 0;
  int*     V_q;
  int      l;

  int      gcbins      = 1;
  int      cpgbins     = 1;
  int      lenbins     = 1;
  int      correctcpg  = 0;
  int      correctgc   = 0;
  int      correctlen  = 0;
  double      c;

  float    divbins     = 50.0;

  // CpG island stuff
  char*    cpgfile = 0;
  int      cpgfile_bins = 1;

   
  int      count_seq;
  char*    count_seq_index;
  int      nb_prev_bad;

  //int      l;
  int      pass = -1;
  int      nborfs_retained;
  int*     M_q_cross;
  int*     E_q_cross;
  float    r;
  double    z_cross;
  int      jn   = 10;
  int      jn_t = 3;
  int      jn_f = 3;

  double    c_cross;
   
  int      fastthreshold      = 1;
  int      idx_eval           = -1;
  int      idx_eval_up;
  int      idx_eval_do;
  int      fastthreshold_jump = 100;

  int      shuffle_rank       = -1;
  int      last_idx_eval;

  int*     seed_pass;
  PT       pt;
  
  int      outputallseeds     = 1;

  char*    firedir            = 0;
  
  // debug
  int      debug              = 1;
  float    max_exp_val        = -1.0;
  float    min_exp_val        = 100000.0;

  int      nbkmers_act        = -1;

  if (argc == 1) {
    printf("Usage : mi_find -expfile FILE -kmerfile FILE -fastafile FILE -k INT -singlestrand [01] -fast [01]\n");
    exit(0);
  }

  //
  //  get FIREDIR environment variable
  //
  firedir = getenv("FIREDIR");
  if (firedir == 0) {
    die("please set the FIREDIR environment variable\n");
  }


  expfile         = get_parameter(argc, argv, "-expfile");

  if (exist_parameter(argc, argv, "-kmerfile")) 
    kmerfile        = get_parameter(argc, argv, "-kmerfile");

  if (exist_parameter(argc, argv, "-report")) 
    report = atoi(get_parameter(argc, argv, "-report"));

  if (exist_parameter(argc, argv, "-outputallseeds")) 
    outputallseeds = atoi(get_parameter(argc, argv, "-outputallseeds"));

  if (exist_parameter(argc, argv, "-fastthreshold")) 
    fastthreshold = atoi(get_parameter(argc, argv, "-fastthreshold"));

  if (exist_parameter(argc, argv, "-fastthreshold_jump")) 
    fastthreshold_jump = atoi(get_parameter(argc, argv, "-fastthreshold_jump"));

  if (exist_parameter(argc, argv, "-divbins")) 
    divbins = atof(get_parameter(argc, argv, "-divbins"));
	
  if (exist_parameter(argc, argv, "-jn")) 
    jn = atoi(get_parameter(argc, argv, "-jn"));

  if (exist_parameter(argc, argv, "-jn_t")) 
    jn_t = atoi(get_parameter(argc, argv, "-jn_t"));

  if (exist_parameter(argc, argv, "-jn_f")) 
    jn_f = atoi(get_parameter(argc, argv, "-jn_f"));
 
  //
  //  conditional MI options
  //
  if (exist_parameter(argc, argv, "-correctcpg")) {
    correctcpg = atoi(get_parameter(argc, argv, "-correctcpg"));
    cpgbins    = 3;
  }

  if (exist_parameter(argc, argv, "-cpgfile")) {
    cpgfile = get_parameter(argc, argv, "-cpgfile");
    cpgfile_bins = 2;
  }


  if (exist_parameter(argc, argv, "-correctgc")) {
    correctgc = atoi(get_parameter(argc, argv, "-correctgc"));
    gcbins    = 3;
  }

  if (exist_parameter(argc, argv, "-correctlen")) {
    correctlen = atoi(get_parameter(argc, argv, "-correctlen"));
    lenbins    = 3;
  }

  if (exist_parameter(argc, argv, "-outfile")) 
    outfile        = get_parameter(argc, argv, "-outfile");


  if (exist_parameter(argc, argv, "-outprm")) 
    outprm        = atoi(get_parameter(argc, argv, "-outprm"));

   
  if (exist_parameter(argc, argv, "-docor")) 
    docor        = atoi(get_parameter(argc, argv, "-docor"));


  fastafile       = get_parameter(argc, argv, "-fastafile");
  kmersize        = atoi(get_parameter(argc, argv, "-k"));



  kmer            = get_parameter(argc, argv, "-kmer");
  pos             = atoi(get_parameter(argc, argv, "-pos"));   
   
  if (exist_parameter(argc, argv, "-seed")) 
    seed = atoi(get_parameter(argc, argv, "-seed"));
  
  
  default_set_seed(seed);
  
      
  if (exist_parameter(argc, argv, "-quantized")) 
    quantized = atoi(get_parameter(argc, argv, "-quantized"));


  if (exist_parameter(argc, argv, "-max_selected")) 
    max_selected = atoi(get_parameter(argc, argv, "-max_selected"));
   
  if (exist_parameter(argc, argv, "-max_nbsdev")) 
    max_nbsdev = atof(get_parameter(argc, argv, "-max_nbsdev"));
   
  
  //if (exist_parameter(argc, argv, "-fast")) 
  //  fast             = atoi(get_parameter(argc, argv, "-fast"));

  if (exist_parameter(argc, argv, "-ebins")) 
    ebins            = atoi(get_parameter(argc, argv, "-ebins"));

  if (exist_parameter(argc, argv, "-mbins")) {
    mbins            = atoi(get_parameter(argc, argv, "-mbins"));
  }

  if (exist_parameter(argc, argv, "-gap")) 
    gap            = atoi(get_parameter(argc, argv, "-gap"));
   
  if (exist_parameter(argc, argv, "-singlestrand")) 
    singlestrand    = atoi(get_parameter(argc, argv, "-singlestrand"));   
   
  if (exist_parameter(argc, argv, "-rna")) 
    singlestrand    = atoi(get_parameter(argc, argv, "-rna"));   
   
  if (exist_parameter(argc, argv, "-verbose")) 
    verbose         = atoi(get_parameter(argc, argv, "-verbose"));   

  if (exist_parameter(argc, argv, "-save_kmer_profile")) 
    save_kmer_profile         = atoi(get_parameter(argc, argv, "-save_kmer_profile"));   
   
  if (kmerfile == 0) {

    if (kmersize == 0) {
      kmersize = 7;
    }
     
    kmerfile = (char*)calloc(1000 , sizeof(char));

    if (singlestrand == 0) {
      sprintf(kmerfile, "%s/FIRE_DATA/KMERFILES/%dmers_dna.txt", firedir, kmersize);
    } else {
      sprintf(kmerfile, "%s/FIRE_DATA/KMERFILES/%dmers_rna.txt", firedir, kmersize);
    }
     
  }
   
 
  readKmers_general (kmerfile, max_nb_kmers, inc_nb_kmers, &kmers, &nbkmers, &kmersize); 


  PT_createPrefixTree(&pt, nbkmers, kmersize);
  PT_populate(&pt, kmers, nbkmers);
  



  if (exist_parameter(argc, argv, "-shuffle")) {
    shuffle = atoi(get_parameter(argc, argv, "-shuffle"));
    
    if (exist_parameter(argc, argv, "-shuffle_rank")) 
      shuffle_rank = atoi(get_parameter(argc, argv, "-shuffle_rank"));
    else
      shuffle_rank = shuffle;
  }

  if (shuffle == -1) {
    
    shuffle      = nbkmers;
    shuffle_rank = shuffle;
    
    if (verbose == 1) {
      printf("Setting -shuffle to %d.\n", shuffle);
    }
    
  } else {
    
    printf("-shuffle set to %d by user.\n", shuffle);

  }

  if ((shuffle == 0) && (verbose == 1)) {
    printf("Warning: -shuffle set to 0.\n");
  }

  if ((jn > 0) && (verbose == 1)) {
    printf("Jack-knife enabled.\n");
  }

  if (verbose == 2)
    printf("Read %d candidate seeds\n", nbkmers);

  kmer_counts = (int*)calloc(nbkmers, sizeof(int));
   
  kmers_array = (Kmer*)malloc(nbkmers * sizeof(Kmer));

  //
  //  read in expression data
  //
  readFloatTable(expfile, &m, &n, &data, &rownames, &colnames, 0, 1);

  if (debug == 1) {
    
    printf("Expression file is a %d by %d matrix.\n", n, m);
    
    for (i=0; i<n; i++) {
      if (data[i][0] > max_exp_val) 
	max_exp_val = data[i][0];
      if (data[i][0] < min_exp_val) 
	min_exp_val = data[i][0];
    }
    
    if (quantized == 0) {
      printf("Min expression value read is %3.2f, max is %3.2f\n", min_exp_val, max_exp_val);
    } else {
      printf("Min expression value read is %d, max is %d\n", (int)min_exp_val, (int)max_exp_val);
    }	
    
  }

  E             = (float* ) malloc( n * sizeof(float));
  gene_kmer_pos = (float**) malloc( n * sizeof(float*));
  gene_kmer     = (short**) malloc( n * sizeof(short*));
  gene_kmer_ori = (short**) malloc( n * sizeof(short*));

  kmer_pos      = (int*)    malloc( 10000 * sizeof(int));
  kmer_ori      = (int*)    malloc( 10000 * sizeof(int));

  idx_gene_name      = (int*)malloc( n * sizeof(int));
  idx_gene_name_rank = (int*)malloc( n * sizeof(int));

  // calculate correlations (does not influence MI values, this is just for information)
  if (docor == 1) {
    seq_gc        = (float*)malloc( n * sizeof( float ));
    seq_cpg       = (float*)malloc( n * sizeof( float ));
    seq_lens      = (float*)malloc( n * sizeof( float ));
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

  if (verbose == 1) {
    printf("Reading sequences ... ");
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
    E[ true_idx ]             = data[idx][0];
    idx_gene_name[ true_idx ] = idx;  
     
    //						
    // build a kmer profile for this sequence .. use regexp or not ? YES, because k-mers are preselected
    //     
    gene_kmer    [ true_idx ] = (short*)calloc(nbkmers, sizeof(short));
    if (gene_kmer[true_idx] == 0) {
      printf("Cannot allocate memory for gene_kmer[%d]\n", true_idx);
      exit(0);
    }

    //
    //  gc content / seq len
    //
    if (docor == 1) { 
      seq_cpg [ true_idx ] = calc_CpG_content(seq, strlen(seq));
      seq_gc  [ true_idx ] = calc_gc_content(seq, strlen(seq));
      seq_lens[ true_idx ] = (float)(strlen(seq));
    }


   
       
    for (i=0; i<nbkmers; i++) {
      gene_kmer    [true_idx][i] = 0;
    }
    
    getPresPosOriVectors(seq, kmersize, gap, &pt, singlestrand, &(gene_kmer[true_idx]), &(gene_kmer_pos[true_idx]), &(gene_kmer_ori[true_idx]), nbkmers); 

    
    
    // increase the kmer count when possible
    for (i=0; i<nbkmers; i++) {
      if ( gene_kmer[true_idx][i] > 0 )
	kmer_counts[i] ++;
    }
    
    

     
    true_idx ++;
     
    free(seq);
    free(name);

  }

  fclose(fp);

  if (verbose == 1) {
    printf("Done\n");
  }

  nborfs = true_idx;
  printf("Number of ORFs = %d\n", nborfs);

  if ((quantized == 0) && (ebins == 0)) {
    ebins = (int)(0.5 + (float)nborfs / ( divbins * mbins ));
    if (verbose == 1)
      printf("Number of bins for expression data will be %d.\n", ebins);
  }

  
  if (quantized == 0) {
    // add a little random number to each value in the E vector
    add_small_values_to_identical_floats(E, nborfs);
    

    //  printf("%d\t%f\n", i, E[i]);
    //}
    //exit(0);
    
  }

   
  //
  //  PROCESS EXPRESSION DATA
  //
  quantize_E(E, nborfs, quantized, &ebins, &E_q, &E_q_bins); 
   

  if (shuffle > 0) {
    a_mi_shu = (double*)malloc(sizeof(double) * shuffle );  
  }

   
  


   
  //						
  //  PROCESS ALL KMERS
  //
  if (verbose == 1) 
    printf("Calculating Mutual Information for all seeds ... ");

  for (i=0; i<nbkmers; i++) {
     
    //
    // get the ith column of the gene_kmer matrix
    //
    M = stof_matrix_column(gene_kmer, i, nborfs);     
    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 
    mymi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins);

    if (quantized == 0) {
      fitSimpleLinearModel(M, E, nborfs, &F, &G);
      kmers_array[i].red   = F;
    } else {
      kmers_array[i].red   = -1;
    }
     
    kmers_array[i].index = i;
    kmers_array[i].score = mymi;
          
    
          
    free(M);
    free(M_q);
  }

  

  qsort((void*)kmers_array, nbkmers, sizeof(Kmer), CmpFunc);

  if (verbose == 1) 
    printf("Done\n");
   
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
      strcat(outrepfile, ".rep");
      fpr = fopen( outrepfile, "w");
    } else {
      fpr = fopen( "report.txt", "w");
    }
  }


  //
  //  two tactics: 1) as before 2) sample
  //

  if (fastthreshold == 1) {

    if (max_selected >= 0)
      nbkmers = max_selected+1;

    if (verbose == 1) 
      printf("Determining threshold ...\n");

    seed_pass = (int*)malloc(nbkmers * sizeof(int));
    for (i=0; i<nbkmers; i++) {
      seed_pass[i] = -1; // don't know
    }	

    idx_eval      =  0;
    last_idx_eval = -1;

    // phase 1, find the lower limit     
    for (idx_eval=0; idx_eval<nbkmers; idx_eval+=fastthreshold_jump) {

      // evaluate idx eval


      pass = evalSeed(kmers_array[idx_eval].index, gene_kmer, nborfs, kmers_array[idx_eval].score, mbins, E_q, ebins, shuffle, jn, jn_f, jn_t); 

      if (verbose == 1) 
	printf("Evaluating seed #%d, %s ... %s.\n", idx_eval, getGappedKmer(kmers[kmers_array[idx_eval].index], gap), (pass==1?"passed":"not passed"));

      //printf("eval %d, pass=%d\n", idx_eval, pass);
      if (pass == 0) {
	seed_pass[ idx_eval] = 0;
	break;
      } else {
	seed_pass[ idx_eval] = 1;
	last_idx_eval = idx_eval;
      }
     
    }
    

    printf("Decreasing intervals phase.\n");
	
    if ( last_idx_eval >= 0 ) {

      // do line-in-desert
      idx_eval_up  = max(0, idx_eval - fastthreshold_jump);
      idx_eval_do  = min(nbkmers - 1, idx_eval);
      
       
      while ( (idx_eval_do - idx_eval_up) > 10) {
       
	idx_eval     = idx_eval_up + (idx_eval_do - idx_eval_up) / 2;

	pass = evalSeed(kmers_array[idx_eval].index, gene_kmer, nborfs, kmers_array[idx_eval].score, mbins, E_q, ebins, shuffle, jn, jn_f, jn_t); 
	
	if (verbose == 1) 
	  printf("Evaluating seed #%d, %s ... %s.\n", idx_eval, getGappedKmer(kmers[kmers_array[idx_eval].index], gap), (pass==1?"passed":"not passed"));

	if (pass == 1) {
	 
	  // go down half interval
	  idx_eval_up   = idx_eval;
	  last_idx_eval = idx_eval; 

	  // record that it passed
	  seed_pass[ idx_eval] = 1;
	 
	} else {
	   
	  // go up half interval
	  idx_eval_do  = idx_eval;

	  // record that it did not pass
	  seed_pass[ idx_eval] = 0;
	}
	 
      }
    
    } 

    //
    //  eval next 10 to be sure
    //
    
    //  last_idx_eval = idx_eval;
    
    printf("Searches for 10 consecutive 'not passed'.\n");

    if (last_idx_eval >= 0) {
      nb_prev_bad   = 0;
    } else {
      nb_prev_bad   = 1;
    }
    
    i             = min(nbkmers-1, last_idx_eval + 10);
	  
    while ((i < nbkmers) && (nb_prev_bad < 10)) {
      
      if (seed_pass[i] != -1) 
	pass = seed_pass[i];
      else
	pass = evalSeed(kmers_array[i].index, gene_kmer, nborfs, kmers_array[i].score, mbins, E_q, ebins, shuffle, jn, jn_f, jn_t); 
      if (verbose == 1) 
	printf("Evaluating seed #%d, %s ... %s.\n", i, getGappedKmer(kmers[kmers_array[i].index], gap), (pass==1?"passed":"not passed"));
      

      if (pass == 1) {
	
	// position of the last good one
	last_idx_eval = i;

	// reset nb bad
	nb_prev_bad = 0;

	// store
	seed_pass[i] = 1;

	// jump 10 down (+1 because the current one is good)
	i += 10 + 1; // add 1 because 1 is going to be removed immediately below
	
      } else {
	nb_prev_bad ++;
	seed_pass[i] = 0;
      }

      
      i --;
    }
    

    if (verbose == 1) {
      if (last_idx_eval >= 0)
	printf("Threshold set to seed #%d (included).\n", last_idx_eval );
      else
	printf("No seeds passed the tests.\n");
    }
    //printf("threshold would be %d\n", last_idx_eval);
    
    for (i=0; i<=last_idx_eval; i++) {
      
      if ((outputallseeds == 0) && (seed_pass[i] == 0)) {
	continue;
      }

      ii = kmers_array[i].index;
      
      //printf("%s", kmers[ii]);
      if (outfile != 0) 
	fprintf(fp, "%s", kmers[ii]);

      printf("%s", getGappedKmer(kmers[ii], gap));
      
      if (outfile != 0)
	fprintf(fp, "\t%5.4f", kmers_array[i].score);
      printf("\t%5.4f", kmers_array[i].score);
      
      if (quantized == 0) {
	if (outfile != 0)
	  fprintf(fp, "\t%4.3f", kmers_array[i].red);
	printf("\t%4.3f", kmers_array[i].red);
      }
      
      if (outfile != 0) 
	fprintf(fp, "\n"); 
      printf("\n");
      
      if (report == 1) {
	
	M = stof_matrix_column(gene_kmer, ii, nborfs);     
	quantize_M_zero_eqpop(M, nborfs, mbins, &M_q);
	
	add_to_report(fpr, kmers[ii], kmersize, gap, M_q, mbins, E_q, ebins, (quantized==1?0:E_q_bins), nborfs);
	
	free(M);
	free(M_q);

      }
      
    }
    
     
  } else {


    nb_prev_bad = 0;

    for (i=0; i<min(nbkmers, max_selected); i++) {

      if (nb_prev_bad == 10) {
	break;
      }

      ii = kmers_array[i].index;
     
      M = stof_matrix_column(gene_kmer, ii, nborfs);     
      quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

      //
      // reshuffle N times the kmer profile						
      //
      if (shuffle > 0) {
     
	for (j=0; j<shuffle; j++) {

	  E_q_shu = shuffleInt(E_q, nborfs);
	  
	  mymi = CalculateMIbasic(M_q, E_q_shu, nborfs, mbins, ebins);

	  a_mi_shu[j] = mymi;       
	  free(E_q_shu);

	}
       
	mi_shu_avg = average_dbl(a_mi_shu, shuffle);
	mi_shu_std = stddev_dbl (a_mi_shu, shuffle);
       
	z = ( kmers_array[i].score - mi_shu_avg ) / mi_shu_std;

	qsort((void*)a_mi_shu, shuffle, sizeof(double), CmpDblRegular);
       
	c  = a_mi_shu[ shuffle_rank - 1 ];
       
	if (kmers_array[i].score < c) {
	  nb_prev_bad ++;  
	  free(M); free(M_q);
	  continue;
	}
       

      }
     
     
      if (jn > 0) {
       
	M_q_cross = (int*)malloc(nborfs * sizeof(int));
	E_q_cross = (int*)malloc(nborfs * sizeof(int));

	pass = 0;
	for (l=0; l<jn; l++) {
	 
	  //printf("jn %d\n", l);

	  // create a list of objects to keep
	  nborfs_retained = 0;
	  for (j=0; j<nborfs; j++) {
	    r =  default_rand(); //rand()/(RAND_MAX+1.0); //printf("%f\n", r);
	    //printf("%f\n", r);
	    if (r > 1.0/jn_f) {
	      M_q_cross[ nborfs_retained ] = M_q[j];
	      E_q_cross[ nborfs_retained ] = E_q[j];
	      nborfs_retained ++;
	    }
	   
	  }
	 
	  //printf("retained %d\n", nborfs_retained);

	  // eval mi
	  mymi = CalculateMIbasic(M_q_cross, E_q_cross, nborfs_retained, mbins, ebins);
	 
	  // shufflings
	  z_cross    = get_zscore_and_rank_value(mymi, M_q_cross, mbins, E_q_cross, ebins, nborfs_retained, shuffle, shuffle_rank, &c_cross); 
	 
	  // check if passes
	  //if ((z_cross > max_nbsdev) && (mymi > c_cross))

	  if (mymi > c_cross) 
	    pass ++;
	 
	}
     
	free(M_q_cross);
	free(E_q_cross);

	if ( pass < jn_t ) {
	  nb_prev_bad ++;  // increase the number of bads
	  free(M); free(M_q);
	  continue;
	}
      }
   
     
      //
      // if we are here, the k-mer passed all tests
      //
     
      if (outfile != 0) 
	fprintf(fp, "%s", kmers[ii]);

      printf("%s", getGappedKmer(kmers[ii], gap));
     
      if (outfile != 0)
	fprintf(fp, "\t%5.4f", kmers_array[i].score);
      printf("\t%5.4f", kmers_array[i].score);
     
      if (quantized == 0) {
	if (outfile != 0)
	  fprintf(fp, "\t%4.3f", kmers_array[i].red);
	printf("\t%4.3f", kmers_array[i].red);
      }
     
      if (shuffle > 0) {
	if (outfile != 0)
	  fprintf(fp, "\t%4.3f", z);
	printf("\t%4.3f", z);
      }
     
      if (jn > 0) {
	if (outfile != 0)
	  fprintf(fp, "\t%d", pass);
	printf("\t%d", pass);
      }
     
      if (outfile != 0)
	fprintf(fp, "\n");
      printf("\n");
     
      if (report == 1) {
	add_to_report(fpr, kmers[ii], kmersize, gap, M_q, mbins, E_q, ebins, (quantized==1?0:E_q_bins), nborfs);
       
      }
     
      nb_prev_bad = 0;  // reset the number of bads
     
      if (nb_prev_bad == 10) {
	break;
      }
     
      free(M); free(M_q);
     
     
    }
   

  }


  free(E_q);

  if (outfile != 0) {
    fclose(fp);
     
    if (outprm == 1) {
      outprmfile = (char*)calloc(1000, sizeof(char));
      strcat(outprmfile, outfile);
      strcat(outprmfile, ".prm");
      fp = fopen(outprmfile, "w"); 

      fprintf(fp, "expfile\t%s\n" , expfile  );
      fprintf(fp, "quantized\t%d\n" , quantized  );
      fprintf(fp, "shuffle\t%d\n" , shuffle  );
      fprintf(fp, "gap\t%d\n" , gap  );
      fprintf(fp, "rna\t%d\n" , singlestrand  );
      fprintf(fp, "fastafile\t%s\n" , fastafile  );
      fprintf(fp, "max_nbsdev\t%f\n" , max_nbsdev  );
      fprintf(fp, "max_selected\t%d\n" , max_selected  );
      //fprintf(fp, "fast\t%d\n" , fast  );
      fprintf(fp, "seed\t%d\n" , seed  );
      fprintf(fp, "ebins\t%d\n" , ebins  );
      fprintf(fp, "mbins\t%d\n" , mbins  );
      fprintf(fp, "sample_size\t%d\n", nborfs);
       
      //
      //  PROCESS BASIC CORRELATIONS
      //
       
      if (docor == 1) { 
	 
	if (quantized == 0) {
	  ebins = (int)( 0.5 + sqrt( nborfs / divbins ) );
	  mbins = ebins;
	} else {
	  mbins = (int)( 0.5 + nborfs / ( divbins * ebins ) );
	}
	 
	quantize_E(E, nborfs, quantized, &ebins, &E_q, NULL); 


	V_q   = Quantize(seq_lens, nborfs, mbins, NULL);
	mymi  = CalculateMIbasic(V_q, E_q, nborfs, mbins, ebins);

	z     = get_zscore(mymi, V_q, mbins, E_q, ebins, nborfs, shuffle); 

	fprintf(fp, "cor_seq_len_mi\t%5.4f\ncor_seq_len_z\t%4.3f\n", mymi, z);
	free(V_q);
	 
	V_q   = Quantize(seq_gc, nborfs, mbins, NULL);
	mymi  = CalculateMIbasic(V_q, E_q, nborfs, mbins, ebins);
	z     = get_zscore(mymi, V_q, mbins, E_q, ebins, nborfs, shuffle); 
     
	fprintf(fp, "cor_gc_mi\t%5.4f\ncor_gc_z\t%4.3f\n", mymi, z);
	free(V_q);


	//
	// CpG rate
	//
	V_q   = Quantize(seq_cpg, nborfs, mbins, NULL);
	mymi  = CalculateMIbasic(V_q, E_q, nborfs, mbins, ebins);
	z     = get_zscore(mymi, V_q, mbins, E_q, ebins, nborfs, shuffle); 
	fprintf(fp, "cor_cpg_mi\t%5.4f\ncor_cpg_z\t%4.3f\n", mymi, z);
	free(V_q);

      }
       
      fclose(fp);
    }
     
  }

  //free(E);
  

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



//
//  get three vectors of presence, position (avg), orientation (avg)
//
void getPresPosOriVectors(char* seq1, int motifsize, int gap, PT* pt, int singlestrand, short** kmer_pre, float** kmer_pos, short** kmer_ori, int nbkmers) {
  
  int   l, i, j;
  int   stop;
  char* mce;
  char* c_mce;
  int   nb1;
  int   nb2;

  int   sid = 0;
  
  mce = (char*)calloc(motifsize+1, sizeof(char));  
  l   = strlen(seq1);

  stop = strlen(seq1) - motifsize - gap + 1;
  if (stop < 0)
    stop = 0;
  
  for (j=sid; j<stop-sid; j++) {
    
    if ((sid == 1) && ((seq1[j-1] == 'N') || (seq1[j+motifsize+1] == 'N')))
      continue;

    if (gap > 0) {
      strncpy(mce, seq1+j,  (motifsize / 2));
      strncpy(mce+(motifsize / 2), seq1+j + (motifsize / 2) + gap, motifsize - (motifsize / 2));
      
    } else {
      strncpy(mce, seq1+j, motifsize);
    }
    
    mce[motifsize] = '\0';
    c_mce          = complement(mce);
    
    nb1 = PT_existWord(pt, mce);

    nb2 = -1;
    if (singlestrand == 0) {
      nb2 = PT_existWord(pt, c_mce);    
    }    

    //printf("%s\n", mce);
    
    // if any of them has been found, no need to go on
    if ((nb1 > 0) || (nb2 > 0)) {
      
      if (nb1 > 0)
	i = (pt->node2index)[nb1];
      else
	i = (pt->node2index)[nb2];
      
      (*kmer_pre)[i] ++;
    }

    free(c_mce);
  }
  
  /*
    for (i=0; i<nbkmers; i++) {
    if (!isnan((*kmer_pos)[i]))
    (*kmer_pos)[i] = (*kmer_pos)[i] / (*kmer_pre)[i];
    }
  */

}



//
//  
//
int evalSeed(int ii, short** gene_kmer, int nborfs, double score, int mbins, int* E_q, int ebins, int shuffle, int jn, int jn_f, int jn_t) 
{
  float* M;
  int    pass1, pass2;
  int*   M_q;
  
 
  M = stof_matrix_column(gene_kmer, ii, nborfs);     
  quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

  pass1 = max_rank_test(score, M_q, mbins, E_q, ebins, nborfs, shuffle, 0, 0, 0, 0); 
  //printf("pass1=%d\n", pass1);
  if (pass1 == 0) {
    free(M); free(M_q);
    return 0;
  }

  // at this stage, if jn == 0, seeds passes.
  if (jn == 0) {
    free(M); free(M_q);
    return 1;
  }

  //printf("jn=%d\tjn_t=%d\n", jn, jn_t);
  pass2 = jacknife_max_rank_test(M_q, mbins, E_q, ebins, nborfs, shuffle, jn, jn_f, jn_t, 0, 0); 
  //printf("pass2=%d\n", pass2);
  if (pass2 == 0) {
    free(M); free(M_q);
    return 0;
  }
  
  free(M); free(M_q);
  
  return 1;
}
