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


#define MIN_KMER_COUNT 20


void  get_motif_profile(char* ch_str, char** seqs, int nborfs, int singlestrand, float** M ); 
void  mask_and_reindex (char* re, int kmersize, float* M, char** seqs, int nborfs, int singlestrand);
double evaluate_motif_mi(char* ch_str, char** seqs, int mbins, int ebins, int* E_q, int nborfs, int singlestrand, float** M ); 

void  array_to_RE      (int* arr, int kmersize, char** encoded_re);


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
  int     pos;

  int         singlestrand = 0;

  int         gap = 0;
   
  int*     E_q;
  int*     E_q_shu;
  float*   E_q_bins = 0;
  int*     M_q;
  double    mymi;

  char**   seqs;
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

  unsigned char** kmer_seq;
  int      max_seq   = 20000;
  int      max_kmers = 500000;
  char**   mykmers;
  char**   mykmers_c;

  int      true_nbkmers;
  int      idx;

  int*     k_inc;
  int*     k_shu;
  int      myk;

  int      quantized = 0;
   
  int      nbrepeats = 10;
  int      nr;
  int      best_count;
  double   z = -100000;
  int*     init_goodkmers;
  int*     newgoodkmers;
  int      ebins = 0;
  int      mbins   = 2;
  int      nb_matched_kmers;
  int      nk;
  double*  a_mi_shu;
  double   mi_shu_avg;
  double   mi_shu_std;
  int      shuffle      = 10000;
  int      max_nb_kmers = 10000, inc_nb_kmers=10000;
  int      add  =  2;
  int      add5 =  1;
  int      add3 =  1;
  int      kmersize_add;
  int      mask = 1;
  float*   oldmis;
  float    minmi;
  float    max_nbsdev = 5;
  char*    outfile = 0;

  int      seed = 12345;
  char*    outprmfile;
  int      outprm = 0;
  int      report = 1;
  char*    outrepfile;

  int      cnt_change;
  
  float*   seq_cond = 0;
  int      max_nb_changes = 0;
  

  int      count_seq;
  char*    count_seq_index;
  
  float    divbins     = 50.0;
     
  int      correct = 0;
  int      hasit   = 0;
  int      lch     = -1;
  double   c;

  int      shuffle_rank = shuffle;
  char**   store_optimized_motifs;
  int      optim_robustness;
  char*    best_best_ch_str_c;

  


  kmer_seq     = (unsigned char**)malloc( max_kmers *  sizeof( unsigned char* ));
  mykmers      = (char**)malloc( max_kmers *  sizeof( char* ));
  mykmers_c    = (char**)malloc( max_kmers *  sizeof( char* ));

  d_ch = (char**)malloc(nbchars * sizeof(char*));

  d_ch[0] = ".";
  d_ch[1] = "A";
  d_ch[2] = "C";
  d_ch[3] = "G";
  d_ch[4] = "T";
  d_ch[5] = "M";
  d_ch[6] = "R";
  d_ch[7] = "W";
  d_ch[8] = "S";
  d_ch[9] = "Y";
  d_ch[10]= "K";
  d_ch[11]= "V";
  d_ch[12]= "H";
  d_ch[13]= "B";
  d_ch[14]= "D";
   
   
  ch = (char**)malloc(nbchars * sizeof(char*));
   
   
  ch[0] = ".";
  ch[1] = "A";
  ch[2] = "C";
  ch[3] = "G";
  ch[4] = "T";
  ch[5] = "[AC]";
  ch[6] = "[AG]";
  ch[7] = "[AT]";
  ch[8] = "[CG]";
  ch[9] = "[CT]";
  ch[10]= "[GT]";
  ch[11]= "[ACG]";
  ch[12]= "[ACT]";
  ch[13]= "[CGT]";
  ch[14]= "[AGT]";
        
 
  if (argc == 1) {
    printf("Usage : mi_optimize_slow -kmerfile FILE -expfile FILE -quantized 0/1 -fastafile FILE -rna 0/1\n");
    exit(0);
  }

  expfile         = get_parameter(argc, argv, "-expfile");
  kmerfile        = get_parameter(argc, argv, "-kmerfile");
  fastafile       = get_parameter(argc, argv, "-fastafile");
 
  kmer            = get_parameter(argc, argv, "-kmer");
  pos             = atoi(get_parameter(argc, argv, "-pos"));   
   
  if (exist_parameter(argc, argv, "-divbins")) 
    divbins = atof(get_parameter(argc, argv, "-divbins"));
  
    
  if (exist_parameter(argc, argv, "-max_nb_changes")) 
    max_nb_changes = atoi(get_parameter(argc, argv, "-max_nb_changes"));
   
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


  if (exist_parameter(argc, argv, "-add")) {
    add5 = atoi(get_parameter(argc, argv, "-add"));
    add3 = add5;
    add  = add3 + add5;
  }

  if (exist_parameter(argc, argv, "-add5")) 
    add5 = atoi(get_parameter(argc, argv, "-add5"));
  
  if (exist_parameter(argc, argv, "-add3")) 
    add3 = atoi(get_parameter(argc, argv, "-add3"));


  if (exist_parameter(argc, argv, "-seed")) 
    seed = atoi(get_parameter(argc, argv, "-seed"));
  default_set_seed(seed);
    

  if (exist_parameter(argc, argv, "-mask")) 
    mask = atoi(get_parameter(argc, argv, "-mask"));


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
  

 
      
  readKmers (kmer, kmerfile, max_nb_kmers, inc_nb_kmers, &kmers, &oldmis, &nbkmers, &kmersize); 

  if (nbkmers == 0) {
    die("Please a non-empty seed file as input\n");
  }
  
  minmi = oldmis[ nbkmers - 1];

  kmersize_add = kmersize + add5 + add3;

  PT_createPrefixTree(&pt, 200000, kmersize_add);
   

  //
  //  read in expression data
  //
  readFloatTable(expfile, &m, &max_seq, &data, &rownames, &colnames, 0, 1);

  



  idx_gene_name = (int*)malloc( max_seq * sizeof(int));
  E             = (float*)malloc( max_seq * sizeof(int));

  if (correct == 1) {
    seq_cond      = (float*)malloc( max_seq * sizeof( float ));
  }

  //
  //  create a hash out of the gene names
  //
  hcreate(100000);
  for (i=0; i<max_seq; i++) {
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
  nborfs  = 0;
  realname = (char*)calloc(100, sizeof(char));

  seqs = (char**) malloc (20000 * sizeof(char*));
  true_nbkmers = 0;

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
    
    if (verbose == 2) {
      printf("%s\n", realname);
    }
     
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
    E[ true_idx ] = data[idx][0];
    free( data[ idx ] );

    seqs[ true_idx ] = strdup( seq );
     

    true_idx ++;
    free(seq);
    
    free(name);
    
  }
  fclose(fp);

  if (verbose == 1)
    printf("finished seq loading\n");


  nborfs = true_idx;

  //					       
  //  determine number of bins
  //
  if ((quantized == 0) && (ebins == 0)) {
    ebins = (int)(0.5 + (float)nborfs / ( divbins * mbins));
    if (verbose == 1)
      printf("number of bins for expression data %d\n", ebins);
  }

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
  
    
  for (nk=0; nk<nbkmers; nk++) {

    if (verbose == 1) {
      printf("optimizing seed %s\n", kmers[nk]);
    }


    // encode the kmer
    init_ch_str[0] = '\0';
    encode_kmer(kmers[nk], kmersize, add5, add3, d_ch, ch, nbchars, &init_ch_array, &init_ch_str);
     
    // 
    // COMPUTE initial profile
    //
    M = (float*)calloc(nborfs,  sizeof(float));
    get_motif_profile(init_ch_str, seqs, nborfs, singlestrand, &M); 

     
    //
    //  quantize motif profile
    //
    quantize_M_zero_eqpop(M, nborfs, mbins, &M_q); 


    for (k=0; k<kmersize_add; k++) {
      k_inc[k] = k;
    }
     

    //									
    //  initial mi value
    //
    init_best_mymi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins); 


    //
    //  to do after some guys have been masked
    //
    if (nk > 0) {

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
       
      if ((init_best_mymi < minmi) || ((shuffle > 0) && (init_best_mymi < c))) {
	if (verbose == 1) {
	  printf(" no optimization ran, because new mi too low (%f < %f) or z too low (%f)\n", init_best_mymi, minmi, z);
	}
	continue;
      }
    }

    free(M_q);
     
    memcpy( best_best_ch_array, init_ch_array, kmersize_add * sizeof(int)); 
    best_best_ch_str = strdup(init_ch_str);
    best_best_mymi   = init_best_mymi;

  
    //
    //  store 
    //
    store_optimized_motifs = (char**)malloc( nbrepeats * sizeof(char*));
     
    for (nr=0; nr<nbrepeats; nr++) {
       
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

	  for (l=0; l<nbchars; l++) {

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
		if (verbose)
		  printf("do not study %s at position %d, because it does not match seed\n", ch[l], myk-add5);
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
	    kmer_count = 0;
	     
	    //
	    //  evaluate its MI
	    //
	    mymi =  evaluate_motif_mi(ch_str, seqs, mbins, ebins, E_q, nborfs, singlestrand, &M);
	 
	    if (verbose == 1) 
	      printf("try %s\n", ch_str);

	    if (verbose == 3)
	      printf("test '%s', mi = %f             \n", ch_str, mymi);
	     
	    if (mymi > best_mymi) {
	      best_mymi   = mymi;
	      best_l      =  l;
	      best_ch_str = strdup(ch_str);
	      best_count  = kmer_count;
	      memcpy( best_ch_array, ch_array, kmersize_add * sizeof(int)); 
	      change = 1;
	       
	      cnt_change ++;
	      
	      if (verbose == 1)
		printf(" improved with '%s', mi = %f             \n", ch_str, mymi);
	  
	      
	      if ((max_nb_changes > 0) && (cnt_change == max_nb_changes)) {
		break;
	      }
	      
	    }
	  
	    free(ch_str);
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


    optim_robustness = 0;
    best_best_ch_str_c = complement(best_best_ch_str);

    for (nr=0; nr<nbrepeats; nr++) {

      if ((strcmp(store_optimized_motifs[nr], best_best_ch_str) == 0) || (strcmp(store_optimized_motifs[nr], best_best_ch_str_c) == 0)) {
	optim_robustness ++;
      }
    }

    free(best_best_ch_str_c);

    //printf("get motif profile for %s\n", best_best_ch_str);

    get_motif_profile(best_best_ch_str, seqs, nborfs, singlestrand, &M); 
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
	mymi    = CalculateMIbasic(M_q, E_q_shu, nborfs, mbins, ebins);
	a_mi_shu[j] = mymi;       
	free(E_q_shu);
      }
      
       
      mi_shu_avg = average_dbl(a_mi_shu, shuffle);
      mi_shu_std = stddev_dbl (a_mi_shu, shuffle);
      z = ( best_best_mymi - mi_shu_avg ) / mi_shu_std;
      free(a_mi_shu);
    }
     
    
    
  
     
    fprintf(fp, "%s\t%5.4f\t%d",  best_best_ch_str, best_best_mymi, optim_robustness);
    if (outfile != 0)
      printf("%s\t%5.4f\t%d",  best_best_ch_str, best_best_mymi, optim_robustness);
    

    if (shuffle > 0) {
      fprintf(fp, "\t%4.3f", z);
      if (outfile != 0)
	printf("\t%4.3f", z);
    }

    fprintf(fp, "\t%s\t%5.4f", kmers[nk], init_best_mymi);
    if (outfile != 0)
      printf("\t%s\t%5.4f", kmers[nk], init_best_mymi);

    fprintf(fp, "\n");
    if (outfile != 0)
      printf("\n");
    
    
    if (report == 1) {
      add_to_report(fpr, best_best_ch_str, strlen(best_best_ch_str), 0, M_q, mbins, E_q, ebins, E_q_bins, nborfs);
    }
    
    //
    //   MASKING
    //
    if (mask == 1) {
      mask_and_reindex(best_best_ch_str, kmersize_add, M, seqs, nborfs, singlestrand);
    }

    free(M_q);
    free(M);

  }


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
      
      fclose(fp);
      
    }
    
    
    
  }
  
  if (report == 1) {
    fclose( fpr );
  }


  return 0;
}


 
void get_motif_profile(char* ch_str, char** seqs, int nborfs, int singlestrand, float** M ) 
{
  int         i;
  for (i=0; i<nborfs; i++) {
    if (re_matches(ch_str, seqs[i], singlestrand)) {
      (*M)[i] = 1.0;
    } else {
      (*M)[i] = 0.0;
    }
  }
}



double evaluate_motif_mi(char* ch_str, char** seqs, int mbins, int ebins, int* E_q, int nborfs, int singlestrand, float** M ) 
{
  double       mymi;
  int*        M_q;

  get_motif_profile(ch_str, seqs, nborfs, singlestrand, M); 
    
  quantize_M_zero_eqpop(*M, nborfs, mbins, &M_q); 

  //
  //  correlate with expression
  //
  mymi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins);
  	   
  free(M_q);

  return mymi;
}



void array_to_RE(int* arr, int kmersize, char** encoded_re) {
  int k;
  *encoded_re = (char*)calloc( kmersize * 4, sizeof( char ));
  for (k=0; k<kmersize; k++) {
    strcat( *encoded_re,  ch[ arr[k] ] ); 
  }
}




//
//  change to only mask the entry in the genome
//
void mask_and_reindex(char* re, int kmersize, float* M, char** seqs, int nborfs, int singlestrand)
{
  int i,j,l,k;
  int* a_pos;
  int* a_ori;
  int  np;
  int  no;


  for (i=0; i<nborfs; i++) {
    
    if (M[i] > 0) {
      
      l = strlen(seqs[i]);

      a_pos = (int*)malloc( l * sizeof(int));
      a_ori = (int*)malloc( l * sizeof(int));
      
      findSites(re, seqs[i], a_pos, &np, a_ori, &no, l, singlestrand, 1);
      
      // traverse all sites, masking them as we go
      if (np > 0) {

	for (j=0; j<np; j++) {

	  for (k=a_pos[j]; k<a_pos[j]+kmersize; k++) {
	    seqs[i][k] = 'N';
	  }

	}
	
      }
      
      free(a_pos);
      free(a_ori);
      
      
    }

  }


}
