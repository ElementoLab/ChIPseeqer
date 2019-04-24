#define _GNU_SOURCE
#include <search.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>

#ifdef USEMPATROL
#include <mpatrol.h>
#endif


#include "dataio.h"
#include "regexp.h"
#include "statistics.h"
#include "prefix.h"
#include "information.h"
#include "mi_library.h"
#include "sequences.h"


#define MIN_KMER_COUNT 20

//double   minCondInfoNormalized(char **A_seq, int** A_q, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, int shuffle, double minr, int v, int* midx); 

int minCondInfoNormalized(char **A_seq, int** A_q, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, int shuffle, double minr, int v, float** A_raw, float* M, int* midx, double* theminratio); 


double   evaluate_motif_mi(char* ch_str, float** M, int* E_q, int mbins, int ebins, int nborfs, char** mykmers, char** mykmers_c, int* masked, int true_nbkmers, unsigned char** kmer_seq, int singlestrand, int* kmer_count, int* overcount, int* goodkmers, int nbgoodkmers, int** newgoodkmers, int* nbnewgoodkmers, int correctcpg, int* V_q, int gcbins ); 
void     get_motif_profile(char* ch_str, float** M, int nborfs, int** matched_kmers, int *nb_matched_kmers, char** mykmers, char** mykmers_c, int true_nbkmers, unsigned char** kmer_seq, int singlestrand ); 
void     encode_kmer(char* kmer, int kmersize, int add5, int add3, char** d_ch, char** ch, int nbchars, int** encoded_kmer, char** encoded_re);

typedef struct _RE_struct {
  char* str;
  int*  arr;
  int   len;
  int   cluster;
  int*  M_q;
} RE_struct;


char**   ch;
char**   d_ch;
int      nbchars = 15;
float**  p_ch;


int main(int argc, char ** argv)
{
     

  // GENERAL
  int i, j;
  char** kmers;
  int kmersize;
  int nbkmers = 0;
  char* kmerfile;

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
  int     m;
  char**  rownames;
  char**  colnames;
  float** data;
   

  int     true_idx = 0;


  int*    idx_gene_name;
   
  float*  M;
  float*  E;

   
  int     verbose = 0;
  int     k;
  char*   kmer;
   
  int     kmer_count;

  int         singlestrand = 0;

  int         gap = 0;
   
  int*     E_q;
  int*     E_q_shu;
  float*   E_q_bins = 0;
  int*     M_q;
  double    mymi;

  int      l;
  int      kk;
   
  char*   ch_str;
  int*    ch_array;
  char*   init_ch_str;
  int*    init_ch_array;
   
  double   init_best_mymi;


   
  
  int best_l;
  double best_mymi;
  char* best_ch_str;
  int*  best_ch_array;
  int*  best_best_ch_array;
  int   change;

  char* best_best_ch_str = 0;
  double best_best_mymi;
	 
  FILE*    fpr = 0;
  

  PT       pt;
  char*    mce;
  char*    c_mce;
  int      stop;
  int      nb1;
  int      nb2;

  unsigned char** kmer_seq;
  //unsigned char** kmer_seq2;

  int      max_seq   = 20000;
  int      max_kmers = 500000;
  char**   mykmers;
  char**   mykmers_c;

  int      true_nbkmers;
  int         idx;

  int*         k_inc;
  int*         k_shu;
  int          myk;

  int          quantized = 0;
   
  int          nbrepeats = 10;
  int          nr;


  int          best_count;

  double       z = -100000;
  int          overcount;
  int*         init_goodkmers;
   
  int*         newgoodkmers;
  int          nbnewgoodkmers;
  int          ebins = 0;
  int          mbins   = 2;

  int*         matched_kmers; 
  int          nb_matched_kmers;
  int          nk;
  double*      a_mi_shu;
  double       mi_shu_avg;
  double       mi_shu_std;
  int          shuffle = 10000;
    
  int*         masked;
  int          max_nb_kmers = 10000, inc_nb_kmers=10000;
  int          add3 = 1;
  int          add5 = 1;
  int          add  = 2;
  int          kmersize_add;

  int          mask = 0;
  float*       oldmis;
  float        minmi;
  float        max_nbsdev = 5;
  int          t = 1;

  char*        outfile = 0;
  int          seed = 12345;
  char*        outprmfile;
  int          outprm = 0;
  int          report = 1;
  char*        outrepfile;

  int          cnt_change;
  
  int*         V_q = 0;
  int          gcbins    = 1;

  int          max_nb_changes = 0;
  
  
  float        divbins     = 50.0;
     
  int          correct = 0;
  int          hasit   = 0;
  int          lch     = -1;
  double        c = -1;

  int          shuffle_rank = shuffle;
  char**       store_optimized_motifs;
  int          optim_robustness;
  char*        best_best_ch_str_c;

  float        maxfreq = 0.5;
  float        myfreq;
  float        lastmyfreq;
  float        best_lastmyfreq;
  
  int**        opt_mot_profiles;   // store the optimized motif profiles
  float**      opt_mot_profiles_nonq; // non-quantized version
  char**       opt_mot_seq;
  int          nb_opt_mot = 0;

  float        minr    = 5.0; // min imi for accepting a motif
  double       minratio;

  char*        logfile;
  int          dolog = 1;
  FILE*        fpl = 0;
  int          didchange = 1;
  int          midx;
  int          cnt_ones;

  int          pass_cond_test = 0;
  int          maxdegeneracy  = -1;
  int          no_neg_cor     = 0;
  int          cntreadseq = 0;
  kmer_seq     = (unsigned char**)malloc( max_kmers *  sizeof( unsigned char* ));
  mykmers      = (char**)malloc( max_kmers *  sizeof( char* ));
  mykmers_c    = (char**)malloc( max_kmers *  sizeof( char* ));
  d_ch         = (char**)malloc(nbchars * sizeof(char*));
  ch           = (char**)malloc(nbchars * sizeof(char*));

  d_ch[0] = "N";  ch[0]   = ".";
  d_ch[1] = "A";  ch[1]   = "A";
  d_ch[2] = "C";  ch[2]   = "C";
  d_ch[3] = "G";  ch[3]   = "G";
  d_ch[4] = "T";  ch[4]   = "T";
  d_ch[5] = "M";  ch[5]   = "[AC]";
  d_ch[6] = "R";  ch[6]   = "[AG]";
  d_ch[7] = "W";  ch[7]   = "[AT]";
  d_ch[8] = "S";  ch[8]   = "[CG]";
  d_ch[9] = "Y";  ch[9]   = "[CT]";
  d_ch[10]= "K";  ch[10]  = "[GT]";
  d_ch[11]= "V";  ch[11]  = "[ACG]";
  d_ch[12]= "H";  ch[12]  = "[ACT]";
  d_ch[13]= "B";  ch[13]  = "[CGT]";
  d_ch[14]= "D";  ch[14]  = "[AGT]";
   
  if (argc == 1) {
    printf("Usage : mi_optimize -kmerfile FILE -expfile FILE -quantized 0/1 -fastafile FILE -rna 0/1\n");
    exit(0);
  }

  expfile         = get_parameter(argc, argv, "-expfile");
  kmerfile        = get_parameter(argc, argv, "-kmerfile");
  fastafile       = get_parameter(argc, argv, "-fastafile");
  kmer            = get_parameter(argc, argv, "-kmer");

   
  if (exist_parameter(argc, argv, "-divbins")) 
    divbins = atof(get_parameter(argc, argv, "-divbins"));
  
  if (exist_parameter(argc, argv, "-max_nb_changes")) 
    max_nb_changes = atoi(get_parameter(argc, argv, "-max_nb_changes"));
  
  if (exist_parameter(argc, argv, "-maxdegeneracy")) {
    maxdegeneracy = atoi(get_parameter(argc, argv, "-maxdegeneracy"));
    printf("Set maxdegeneracy to %d.\n", maxdegeneracy);
  }

  if (exist_parameter(argc, argv, "-no_neg_cor")) {
    no_neg_cor = atoi(get_parameter(argc, argv, "-no_neg_cor"));
    if (no_neg_cor == 1) 
      printf("Ignore negatively correlated motifs.\n");
  }


  if (exist_parameter(argc, argv, "-nbrepeats")) 
    nbrepeats = atoi(get_parameter(argc, argv, "-nbrepeats"));

  if (exist_parameter(argc, argv, "-report")) 
    report = atoi(get_parameter(argc, argv, "-report"));
    

  if (exist_parameter(argc, argv, "-max_nbsdev")) 
    max_nbsdev = atoi(get_parameter(argc, argv, "-max_nbsdev"));
    
  if (exist_parameter(argc, argv, "-outfile")) 
    outfile        = get_parameter(argc, argv, "-outfile");

  if (exist_parameter(argc, argv, "-outprm")) 
    outprm        = atoi(get_parameter(argc, argv, "-outprm"));

  if (exist_parameter(argc, argv, "-log")) 
    dolog = atoi(get_parameter(argc, argv, "-log"));

  if (exist_parameter(argc, argv, "-add")) {
    add3 = atoi(get_parameter(argc, argv, "-add"));
    add5 = add3;
    add  = add3 + add5;
  }

  
  if (exist_parameter(argc, argv, "-add5")) 
    add5 = atoi(get_parameter(argc, argv, "-add5"));
  
  if (exist_parameter(argc, argv, "-add3")) 
    add3 = atoi(get_parameter(argc, argv, "-add3"));


  if (exist_parameter(argc, argv, "-maxfreq")) 
    maxfreq = atof(get_parameter(argc, argv, "-maxfreq"));


  if (exist_parameter(argc, argv, "-seed")) 
    seed = atoi(get_parameter(argc, argv, "-seed"));
  default_set_seed(seed);


  if (exist_parameter(argc, argv, "-mask")) 
    mask = atoi(get_parameter(argc, argv, "-mask"));

  if (exist_parameter(argc, argv, "-minr")) 
    minr = atof(get_parameter(argc, argv, "-minr"));

  if (exist_parameter(argc, argv, "-t")) 
    t = atoi(get_parameter(argc, argv, "-t"));

  if (exist_parameter(argc, argv, "-quantized")) 
    quantized = atoi(get_parameter(argc, argv, "-quantized"));
 
  if (exist_parameter(argc, argv, "-shuffle")) 
    shuffle = atoi(get_parameter(argc, argv, "-shuffle"));
 
  if (exist_parameter(argc, argv, "-shuffle_rank")) 
    shuffle_rank = atoi(get_parameter(argc, argv, "-shuffle_rank"));


  if (exist_parameter(argc, argv, "-gap")) 
    gap            = atoi(get_parameter(argc, argv, "-gap"));
   
  if (exist_parameter(argc, argv, "-rna")) 
    singlestrand    = atoi(get_parameter(argc, argv, "-rna"));   

  if (exist_parameter(argc, argv, "-singlestrand")) 
    singlestrand    = atoi(get_parameter(argc, argv, "-singlestrand"));   


  if (exist_parameter(argc, argv, "-verbose")) 
    verbose         = atoi(get_parameter(argc, argv, "-verbose"));   
 
   
  if (exist_parameter(argc, argv, "-ebins")) 
    ebins            = atoi(get_parameter(argc, argv, "-ebins"));
   
  if (exist_parameter(argc, argv, "-mbins")) {
    mbins            = atoi(get_parameter(argc, argv, "-mbins"));
  }
  

  printf("Using minr = %5.4f.\n", minr);
      
  readKmers (kmer, kmerfile, max_nb_kmers, inc_nb_kmers, &kmers, &oldmis, &nbkmers, &kmersize); 

  if (nbkmers == 0) {
    die("Please a non-empty seed file as input\n");
  }
  
  minmi = oldmis[ nbkmers - 1];

  kmersize_add = kmersize + add3 + add5;

  printf("Motif length will be %d.\n", kmersize_add);


  printf("Creating prefix tree for storing patterns ... ");
  PT_createPrefixTree(&pt, 200000, kmersize_add);
  pt.n_allowed = 0;

  printf("Done.\n");

  //
  //  read in expression data
  //

  printf("Reading expression data ... ");
  readFloatTable(expfile, &m, &max_seq, &data, &rownames, &colnames, 0, 1);
  printf("Done.\n");
  
  //printf("Found %d genes in expression file.\n", max_seq);


  idx_gene_name = (int*)malloc( max_seq * sizeof(int));
  if (idx_gene_name == 0) {
    die("Cannot allocate idx_gene_name ... \n");
  }
  E             = (float*)malloc( max_seq * sizeof(int));
  if (E == 0) {
    die("Cannot allocate E ... \n");
  }

  //
  //  create a hash out of the gene names
  //
  printf("Creating hash table for gene names ... ");
  hcreate(100000);
  for (i=0; i<max_seq; i++) {
    e.key = strdup(rownames[i]); 
    e.data = (char*)i;
    hsearch(e, ENTER);
  }
  printf("Done.\n");


  //
  //  read in the sequences that are also in the microarray
  //

  printf("Reading sequences and building the index.\n");
  fp = fopen(fastafile, "r");
  if (!fp) {
    printf("cannot open %s ..\n", fastafile);
    exit(0);
  }

  nextSequence_currentLine = (char*)malloc(200000 * sizeof(char));
  nborfs  = 0;
  realname = (char*)calloc(100, sizeof(char));

  true_nbkmers = 0;

  cntreadseq = 0;

  while ( (seq = nextSequence(fp, &name, &size, &nextSequence_started, &nextSequence_ended, nextSequence_currentLine)) ) {   
    
    //if (cntreadseq % 1000 == 0)
    //  printf("%d           \r", cntreadseq);
    //printf("name = %s, seq = %s\n", name, seq);
    
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

    idx = (int)ep->data;       

    if (true_idx > max_seq) {
      printf("non-sense: true_idx=%d > max_seq=%d\n", true_idx, max_seq);
    }

    E[ true_idx ] = data[idx][0];
    free( data[ idx ] );

    l   = strlen(seq);

    mce = (char*)calloc(kmersize_add, sizeof(char));  
    if (mce == 0) {
      die("sorry cannot allocate mce\n");
    }

    stop = l - kmersize_add - gap + 1;

    if (stop < 0)
      stop = 0;
    for (j=0; j<stop; j++) {

      if (gap > 0) {
	strncpy(mce, seq+j,  (kmersize_add / 2));
	strncpy(mce   + (kmersize_add / 2), 
		seq+j + (kmersize_add / 2) + gap, kmersize_add - (kmersize_add / 2)); 
      } else {
	strncpy(mce, seq + j, kmersize_add);
      }

      mce[kmersize_add] = '\0';

      c_mce  = complement(mce);
      nb1    = PT_existWord(&pt, mce);
      nb2    = -1;

      if (singlestrand == 0) {
	nb2 = PT_existWord(&pt, c_mce);    
      }    
      if ((nb1 != -2 ) && (nb2 != -2)) {
	if (nb1 > 0) 
	  i = nb1;
	else if (nb2 > 0)
	  i = nb2;
	else {
	  i = PT_addWord(&pt, mce);
	  
	  mykmers [ true_nbkmers ] = (char*)calloc( kmersize_add+1, sizeof( char ));
	  memcpy(mykmers[true_nbkmers], mce, kmersize_add); 
	  mykmers[true_nbkmers][kmersize_add] = '\0';

	  if (singlestrand == 0) {
	    mykmers_c [ true_nbkmers ] = (char*)calloc( kmersize_add+1, sizeof( char ));
	    memcpy(mykmers_c[true_nbkmers], c_mce, kmersize_add); 
	    mykmers_c[true_nbkmers][kmersize_add] = '\0';
	  }
	  
	  kmer_seq[ true_nbkmers ] = create_binarized_array(max_seq);

	  if (kmer_seq[ true_nbkmers ] == 0) {
	    die("sorry cannot allocate memory for kmer_seq[ true_nbkmers ]\n");
	  }
	  

	  (pt.node2index)[ i ]     = true_nbkmers;
	  true_nbkmers ++;
	}
	 
	idx = (pt.node2index)[ i ];
	 
	if (idx >= max_kmers) {
	  die("max_kmers reached, die ..\n");
	}
      
	//what follows replaced 
	//kmer_seq2[ idx ][ true_idx ] = 1;

	set_entry_in_binarized_array(kmer_seq[ idx ], true_idx);

	
      }
       
      free(c_mce);
    }
    true_idx ++;
    free(seq);
    
    free(name);
    
  }
  fclose(fp);

  
  printf("%d k-mers in index\n", true_nbkmers);
  
  
  cnt_ones = 0;
  for (i=0; i<true_nbkmers; i++) {
    for (j=0; j<true_idx; j++) {
      if (get_entry_in_binarized_array(kmer_seq[i], j) > 0)
	cnt_ones ++;
    }
  }
  
  printf("Matrix occupancy = %d / %f\n", cnt_ones, (double)true_nbkmers * true_idx);
  

  printf("Done.                                \n");
  
  nborfs = true_idx;
  
  //					       
  //  determine number of bins
  //
  if ((quantized == 0) && (ebins == 0)) {
    ebins = (int)(0.5 + (float)nborfs / ( divbins * mbins * gcbins));
  }

   
  if (quantized == 0) {
    // add a little random number to each value in the E vector
    add_small_values_to_identical_floats(E, nborfs);
  }

  


  //					      
  //  create an index of masked kmers						
  //
  masked = (int*)calloc(true_nbkmers, sizeof(int));


  //						
  //  quantize expression profile
  //

  quantize_E(E, nborfs, quantized, &ebins, &E_q, &E_q_bins); 
  
  //
  //  memory alloc for init kmer, good kmer store, etc
  //
  init_ch_array      = (int*) malloc(kmersize_add * sizeof(int));   
  init_ch_str        = (char*)calloc(100, sizeof(char));
  ch_array           = (int*) malloc(kmersize_add * sizeof(int));
  best_ch_array      = (int*) malloc(kmersize_add * sizeof(int));
  best_best_ch_array = (int*) malloc(kmersize_add * sizeof(int));
  k_inc              = (int*) malloc(kmersize_add * sizeof(int));
  init_goodkmers     = (int*) malloc(true_nbkmers * sizeof(int));
  newgoodkmers       = (int*) malloc(true_nbkmers * sizeof(int));

  // 
  //  loop over the k-mers
  //


  if (outfile != 0) {
      fp = fopen(outfile, "w");
      if (fp == 0) {
	printf("cannot open %s\n", outfile);
	exit(0);
      }
  } else {
    fp = stdout;
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
  
  
  if (dolog == 1) {
    if (outfile != 0) {
      logfile = (char*)calloc(1000, sizeof(char));
      strcat(logfile, outfile);
      strcat(logfile, ".log");
      fpl = fopen( logfile, "w");
    } else {
      fpl = fopen( "log.txt", "w");
    }
  } 

  
  // 
  opt_mot_profiles      = (int**)   malloc(nbkmers * sizeof(int*)  );
  opt_mot_profiles_nonq = (float**) malloc(nbkmers * sizeof(float*));
  opt_mot_seq           = (char**)  malloc(nbkmers * sizeof(char*) );
    
  for (nk=0; nk<nbkmers; nk++) {

    if (verbose == 1) {
      printf("Studying %s.\n", kmers[nk]);
    }

    // encode the kmer
    init_ch_str[0] = '\0';
    encode_kmer(kmers[nk], kmersize, add5, add3, d_ch, ch, nbchars, &init_ch_array, &init_ch_str);
     
    // 
    // COMPUTE initial profile
    //
    M = (float*)calloc(nborfs,  sizeof(float));
    get_motif_profile(init_ch_str, &M, nborfs, &matched_kmers, &nb_matched_kmers, mykmers, mykmers_c, true_nbkmers, kmer_seq, singlestrand ); 
    free(matched_kmers);

    //
    // calculate the start frequency
    //
    kmer_count = 0;
    for (i=0; i<nborfs; i++) {
      if (M[i] > 0) 
	kmer_count ++;
    }
    lastmyfreq = kmer_count / (float)nborfs;
    best_lastmyfreq = lastmyfreq;
    
     
    //
    //  quantize motif profile
    //
    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 


    //
    //  the initial number of good kmer is all kmers
    //
    for (i=0; i<true_nbkmers; i++) {
      init_goodkmers[ i ] = i;
    }
    for (k=0; k<kmersize_add; k++) {
      k_inc[k] = k;
    }
     
    //
    //  check how much information it adds to the previous guys
    //
    printf("\nEvaluating seed %s ... \n", getGappedKmer(kmers[nk], gap));
    
    
    if (no_neg_cor == 1) {

      // pearson correlation
      float pe      = pearson(M, E, nborfs);

      printf("pe = %3.2f\n", pe);

      if (pe <  1e-3) {
	printf("Negatively correlated motif, skipping.\n");
	continue;
      }

    }
    
    pass_cond_test = minCondInfoNormalized(opt_mot_seq, opt_mot_profiles, mbins, nb_opt_mot, 
					   M_q, mbins, E_q, ebins, nborfs, shuffle, minr, 1, 
					   opt_mot_profiles_nonq, M, &midx, &minratio); 

    
    if ((pass_cond_test == 0) && (nb_opt_mot > 0)) {
      printf("not optimized, minratio = %4.3f for %s (< minr=%4.3f).\n", minratio, opt_mot_seq[midx], minr); 
      fprintf(fpl, "%s not optimized, minratio = %4.3f for %s (< minr=%4.3f).\n", kmers[nk], minratio, opt_mot_seq[midx], minr); 
      fflush(fpl);

      continue;
    } else {
      printf("optimizing.\n"); // because r_min > minr.\n");
    }

    //									
    //  initial mi value
    //
    init_best_mymi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins); 

     if (dolog == 1) {
       fprintf(fpl, "%s\t%d\t%d\t%s\t%4.3f\t%4.3f\n", kmers[nk], -1, -1, init_ch_str, init_best_mymi, lastmyfreq);
     }
	 
    
    if (nk > 0) {
      
      
      //
      //  this is obsolete and 
      //
      if (shuffle > 0) {
	 
	//
	// obtain the profile for the optimized motif
	//
	 
	a_mi_shu = (double*)malloc(sizeof(double) * shuffle );  
	 
	for (j=0; j<shuffle; j++) {
	  E_q_shu = shuffleInt(E_q, nborfs); 
	  mymi = CalculateMIbasic(M_q, E_q_shu, nborfs, mbins, ebins);
	  a_mi_shu[j] = mymi;       
	  free(E_q_shu);
	}
	 
	 
	mi_shu_avg = average_dbl(a_mi_shu, shuffle);
	mi_shu_std = stddev_dbl (a_mi_shu, shuffle);
	z = ( init_best_mymi - mi_shu_avg ) / mi_shu_std;

	qsort((void*)a_mi_shu, shuffle, sizeof(double), CmpDblRegular);
	c  = a_mi_shu[ shuffle_rank - 1 ];


	free(a_mi_shu);
      }
      
    }

    free(M_q);

    //					       
    //  loop thru nb of repeats
    //
     
    memcpy( best_best_ch_array, init_ch_array, kmersize_add * sizeof(int)); 
    best_best_ch_str = strdup(init_ch_str);
    best_best_mymi   = init_best_mymi;

  
    //
    //  store 
    //
    store_optimized_motifs = (char**)malloc( nbrepeats * sizeof(char*));
     
    for (nr=0; nr<nbrepeats; nr++) {
      
      if (verbose == 1)
	printf("Repeat %d.\n", nr);
       
      best_mymi   = init_best_mymi;
      best_ch_str = strdup(init_ch_str);
      best_count  = -1;
      memcpy(ch_array, init_ch_array, kmersize_add * sizeof(int));


      cnt_change = 0;
	
      while (1) {
	 
	change = 0;
	 

	
	//
	// create a random index
	//
	k_shu = shuffleInt(k_inc, kmersize_add); 
         
	//
	// change each column
	//
	for (k=0; k<kmersize_add; k++) {
	   
	  myk       = k_shu[ k ];
	  best_l    = ch_array[ myk ];
	   
	  //
	  //  change to everything possible, if it is compatible with the initial seed
	  //

	  didchange = 0;

	  for (l=0; l<nbchars; l++) {
	    
	    // 
	    // the following code is to make sure that the RE stays compatible with the original motif
	    //
	    if ((l != 0) && (myk > add5-1) && (myk<kmersize_add-add3)) {
	      lch = strlen(ch[l]);
	      hasit = 0;
	      for (kk=0; kk<lch; kk++) {
		if ( kmers[nk][ myk-add5] == ch[l][kk] ) {
		  hasit = 1;
		  break;
		}
	      }
	      if (hasit == 0) {
		continue;
	      }

	    }

	    //
	    // this code make sure that the degeneracy level is below threshold
	    //
	    if (maxdegeneracy > 0) {
	      lch = strlen(ch[l]);
	      
	      // [XX] and [XXX] should have two less characters
	      if (lch >= 4) {
		lch = lch - 2;
	      } 
	      if (lch > maxdegeneracy) {
		continue;
	      }
	    }


	    ch_array[ myk ] = l;   // corresponding letter is ch[ ch_array[kk] ]
	    ch_str          = (char*)calloc(100, sizeof(char));

	    // 
	    // reconstitute a new RE
	    //
	    for (kk=0; kk<kmersize_add; kk++) {
	      strcat(ch_str, ch[ ch_array[kk] ]);
	    }



	    //
	    //  evaluate its MI
	    //
	    kmer_count = 0;
	    mymi =  evaluate_motif_mi(ch_str, &M, E_q, mbins, ebins, nborfs, mykmers, mykmers_c, masked, true_nbkmers, kmer_seq, singlestrand, &kmer_count, &overcount, (l==0?init_goodkmers:newgoodkmers), (l==0?true_nbkmers:nbnewgoodkmers), (l==0?(&newgoodkmers):NULL), (l==0?&nbnewgoodkmers:0), (correct==1?1:0), V_q, gcbins);
	    	     
	    
	    
	    kmer_count = 0;
	    for (i=0; i<nborfs; i++) {
	      if (M[i] > 0) 
		kmer_count ++;
	    }
	    
	    

	    myfreq = kmer_count / (float)nborfs;

	     
	     //
	    //  evaluate its MI wrt other found motifs
	    //
	    
	    pass_cond_test = 1;  // pass by default
	    if (nb_opt_mot > 0) {

	      quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 
	      pass_cond_test = minCondInfoNormalized(opt_mot_seq, opt_mot_profiles, mbins, nb_opt_mot, M_q, mbins, E_q, ebins, nborfs, shuffle, minr, 0, opt_mot_profiles_nonq, M, &midx, &minratio); 
	      free(M_q);

	    } 

	    //
	    // accept ONLY if mi increases AND (freq < T OR freq decreases) AND pass constraints
	    //
	    if ((mymi > best_mymi) && ((myfreq < maxfreq) || (myfreq < lastmyfreq)) && (pass_cond_test == 1)) {

	      best_mymi   = mymi;
	      best_l      =  l;
	      best_ch_str = strdup(ch_str);
	      best_count  = kmer_count;
	      memcpy( best_ch_array, ch_array, kmersize_add * sizeof(int)); 
	      change      = 1;
	      cnt_change ++;
	      
	      if (verbose == 1)
		printf("Improved '%s', mi = %f             \n", ch_str, mymi);
	  
	      lastmyfreq  = myfreq;	   
	      didchange   = 1;
   
	      if ((max_nb_changes > 0) && (cnt_change == max_nb_changes)) {
		break;
	      }
	      
	    }
	  
	    free(ch_str);
	  }
	  
	  
	  if (dolog == 1) {
	    if (didchange == 1) 
	      fprintf(fpl, "%s\t%d\t%d\t%s\t%4.3f\t%4.3f\n", kmers[nk], nr, k, best_ch_str, best_mymi, lastmyfreq);
	    else 
	      fprintf(fpl, "%s\t%d\t%d\n", kmers[nk], nr, k);
	  }
	 
	  ch_array[ myk ] = best_l;
	   
	  if ((max_nb_changes > 0) && (cnt_change == max_nb_changes)) {
	    change = 0;
	    break;
	  }
	  
	   
	}
	 
	free(k_shu);
	 
	if (change == 0) {
	   
	  //
	  //  store the motif
	  //
	  store_optimized_motifs[nr] = strdup( best_ch_str );


	  //
	  //  
	  //
	  if (best_mymi > best_best_mymi) {
	    best_lastmyfreq  = lastmyfreq;
	    best_best_ch_str = strdup( best_ch_str );
	    best_best_mymi   = best_mymi;
	    memcpy( best_best_ch_array, best_ch_array, kmersize_add * sizeof(int)); 
	  }
	   
	  break;
	   
	}	
	
	
      }  // end while(1)
       
      free(best_ch_str);
    }  // end loop over nbrepeats
     
     
    
   

     
    //
    //  get motif profile
    //
     
    nb_matched_kmers = 0;
    for (i=0; i<nborfs; i++) {
      M[i] = 0.0;
    }


    //printf("Frequency of the optimal motif is %3.2f\n", best_lastmyfreq);

    optim_robustness = 0;
    best_best_ch_str_c = complement(best_best_ch_str);

    for (nr=0; nr<nbrepeats; nr++) {

      if ((strcmp(store_optimized_motifs[nr], best_best_ch_str) == 0) || (strcmp(store_optimized_motifs[nr], best_best_ch_str_c) == 0)) {
	optim_robustness ++;
      }
    
    }

    free(best_best_ch_str_c);

    get_motif_profile(best_best_ch_str, &M, nborfs, &matched_kmers, &nb_matched_kmers, mykmers, mykmers_c, true_nbkmers, kmer_seq, singlestrand ); 

    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 

    //					       
    //  shuffle
    //
    if (shuffle > 0) {
       
      //
      // obtain the profile for the optimized motif
      //
       
      a_mi_shu = (double*)malloc(sizeof(double) * shuffle );  
       
      for (j=0; j<shuffle; j++) {
	E_q_shu = shuffleInt(E_q, nborfs); 
	
	mymi = CalculateMIbasic(M_q, E_q_shu, nborfs, mbins, ebins);
	a_mi_shu[j] = mymi;       
	free(E_q_shu);
      }
      
       
      mi_shu_avg = average_dbl(a_mi_shu, shuffle);
      mi_shu_std = stddev_dbl (a_mi_shu, shuffle);
      z = ( best_best_mymi - mi_shu_avg ) / mi_shu_std;
      free(a_mi_shu);
    }     
    
    //						
    // store this motif 
    //
    opt_mot_seq     [ nb_opt_mot ] = strdup(best_best_ch_str);

    //					       
    // store the motif profile
    //
    opt_mot_profiles[ nb_opt_mot ]      = M_q;
    opt_mot_profiles_nonq[ nb_opt_mot ] = M;
    nb_opt_mot ++;
  
     
    fprintf(fp, "%s\t%5.4f\t%d",  getGappedMotif(best_best_ch_str, kmersize_add, gap), best_best_mymi, optim_robustness);
    
   
    
    if (outfile != 0)
      printf("Best optimized motif: %s\nMutual information: %f\nOptimization stability: %d",  getGappedMotif(best_best_ch_str, kmersize_add, gap), best_best_mymi, optim_robustness);
      
    if (shuffle > 0) {
      fprintf(fp, "\t%4.3f", z);
      if (outfile != 0)
	printf("\nZ-score: %4.3f", z);
    }

    fprintf(fp, "\t%s\t%5.4f", getGappedKmer(kmers[nk], gap), init_best_mymi);
    if (outfile != 0)
      printf("\nStarting seed: %s\nSeed MI: %f", getGappedKmer(kmers[nk], gap), init_best_mymi);

    fprintf(fp, "\n");
    if (outfile != 0)
      printf("\n");
    
    
    if (report == 1) {
      add_to_report(fpr, best_best_ch_str, strlen(best_best_ch_str), gap, M_q, mbins, E_q, ebins, E_q_bins, nborfs);
    }

    //free(M);


  } // end loop over seeds


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
      fprintf(fp, "seed\t%d\n" , seed  );
      fprintf(fp, "ebins\t%d\n" , ebins  );
      fprintf(fp, "mbins\t%d\n" , mbins  );
      fprintf(fp, "sample_size\t%d\n", nborfs);
      fprintf(fp, "minr\t%f\n", minr);

      fclose(fp);
      
    }
    
  }
  
  if (report == 1) {
    fclose( fpr );
  }

  if (dolog == 1)
    fclose(fpl);

  return 0;
}



 
void get_motif_profile(char* ch_str, float** M, int nborfs, int** matched_kmers, int *nb_matched_kmers, char** mykmers, char** mykmers_c, int true_nbkmers, unsigned char** kmer_seq, int singlestrand ) 
{

  int         i, j;
  pcre*       re;
  const char* error;
  int         erroffset;

  *matched_kmers    = (int*)malloc( true_nbkmers * sizeof(int));
  *nb_matched_kmers = 0;

  for (i=0; i<nborfs; i++) {
    (*M)[i] = 0.0;
  }
  
  re = pcre_compile(ch_str, 0, &error, &erroffset, NULL);

  for (j=0; j<true_nbkmers; j++) {

    if (raw_re_matches(re, mykmers[j], mykmers_c[j], singlestrand) == 1) {

      for (i=0; i<nborfs; i++) {
	// replaces (*M)[i] += (float)(kmer_seq[j][i]);
	(*M)[i] += get_entry_in_binarized_array(kmer_seq[j], i);
      }

      (*matched_kmers)[ *nb_matched_kmers ] = j;
      (*nb_matched_kmers) ++;
    }

  }
  
  pcre_free(re);
  
}



double evaluate_motif_mi(char* ch_str, float** M, int* E_q, int mbins, int ebins, int nborfs, char** mykmers, char** mykmers_c, int* masked, int true_nbkmers, unsigned char** kmer_seq, int singlestrand, int* kmer_count, int* overcount, int* goodkmers, int nbgoodkmers, int** newgoodkmers, int* nbnewgoodkmers, int correct, int* V_q, int gcbins  ) 
{

  int         i, j, k;
  int         my_kmer_count = 0;
  pcre*       re;
  const char* error;
  int         erroffset;
  double      mymi;
  int*        M_q;

  int         m = -1;

  for (i=0; i<nborfs; i++) {
    (*M)[i] = 0.0;
  }
  
  re = pcre_compile(ch_str, 0, &error, &erroffset, NULL);

  if (newgoodkmers != 0) {
    *nbnewgoodkmers = 0;
  }
  
  // loop over the index
  for (k=0; k<nbgoodkmers; k++) {

    j = goodkmers[k];  // get the index for the kmer
    
    //
    // if that k-mer is not masked, and matches the RE, add it to the profile
    //
    if ((masked[j] == 0) && (raw_re_matches(re, mykmers[j], mykmers_c[j], singlestrand) == 1)) {
      
      if (newgoodkmers != 0) {
	(*newgoodkmers)[ *nbnewgoodkmers ] = j;
	(*nbnewgoodkmers) ++;
      }

      for (i=0; i<nborfs; i++) {
	// what follows replaces (*M)[i] += (float)(kmer_seq[j][i]);
	m = get_entry_in_binarized_array(kmer_seq[j], i);
	(*M)[i] += m;

	if (m > 0) {
	  my_kmer_count ++;
	}
      }
      
    }
  }
  
  pcre_free(re);
	   
  if (my_kmer_count < MIN_KMER_COUNT)
    *overcount = 1;
  else 
    *overcount = 0;
    
  quantize_M_zero_eqpop(*M, nborfs, mbins, &M_q); 

  //
  //  correlate with expression
  //
  mymi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins);
  
  
	   
  free(M_q);


  *kmer_count = my_kmer_count;
  return mymi;
}





//
//  returns whether new motif B_q (B) correlates with the A motifs 
//
int minCondInfoNormalized(char **A_seq, int** A_q, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, int shuffle, double minr, int v, float** A_raw, float* M, int* midx, double* theminratio) 
{
  int    i;
  int*   AB_q;
  int    DAE, DAB;
  double mi_ab_e, mi_a_e, mi_a_b;
  double cmi;
  double minratio = 1e6;
  double ratio;
  double pe_t = 1e-3;
  double pe_a_b;

  DAE    = DA * DE;
  DAB    = DA * DB;
  *midx  = -1;

  for (i=0; i<nA; i++) {        

    AB_q        = combineQuantizedVectors(A_q[i], B_q, n, DA, DB);    
    mi_ab_e     = CalculateMIbasic       (AB_q,   E_q, n, DAB, DE);
    mi_a_e      = CalculateMIbasic       (A_q[i], E_q, n, DA, DE);
    cmi         = mi_ab_e - mi_a_e;
    mi_a_b      = CalculateMIbasic       (A_q[i], B_q, n, DA, DB);
  
    if (cmi < 1e-16)
      cmi = 0.0; 
    else if (mi_a_b < 1e-16)
      mi_a_b = 1e-16;
    
    ratio       = cmi / mi_a_b;

    if (i == 0)
      minratio = ratio;
    else {
      
      if (ratio < minratio) {
	minratio = ratio;
      }
    }

    free(AB_q);

    if (v == 1) {
      
      printf("Motif %s : I(B;E|A)/I(A;B)=%5.4f\n", A_seq[i], cmi/mi_a_b);    
      //printf("i=%d I(A,B;E)=%1.4f, I(A;E)=%1.4f, I(B;E|A)=%1.4f, I(A;B)=%1.4f, I(B;E|A)/I(A;B)=%5.4f, s:%s\n", i, mi_ab_e, mi_a_e, cmi, mi_a_b, cmi/mi_a_b, A_seq[i]);    
    }
    

    //					       
    //  break early feature 
    //
    if (minratio < minr) {

      // pearson correlation
      pe_a_b      = pearson(A_raw[i], M, n);
      
      // essentially, if minr low, and if correlation positive, we need to break (motif not good) ..
      if (pe_a_b > pe_t) {
	
	*midx = i;
	
	if (v == 1)
	  printf("Break early (motif is too close to existing optimized motif).\n");
	
	*theminratio = minratio;

	// 0 means the motif does not pass
	return 0;

      }
    }

    
    
  }
  
  // 1 means the currently examined motif is ok
  return 1;
  
}






				      
