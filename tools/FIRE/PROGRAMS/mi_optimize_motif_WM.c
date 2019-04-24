// 
// input a regular expression (even kmer), an expression file, a fasta file
//   somthing very similar to mi_optimize in fact
// 


#include <search.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>


#include "dataio.h"
#include "regexp.h"
#include "statistics.h"
#include "prefix.h"
#include "information.h"
#include "mi_library.h"
#include "sequences.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

void findMostInformativeWMThreshold(float* scores, int* E_q, int nborfs, int ebins, float* t, float* mi);




int main(int argc, char ** argv)
{
     

  // GENERAL
  int i, j;
  char** motifs;
  int    motifsize;
  int    nbmotifs    = 0;
  char*  motiffile = 0;

  int    nborfs = 0;

  FILE   *fp_outfile;
  int    size;
  char*  seq;
  char*   name;

  char*   fastafile;

  ENTRY   e;
  ENTRY   *ep;
  char*   realname;

  char*   expfile;
  int     m;
  char**  rownames;
  char**  colnames;
  float** data;
   

  int     true_idx = 0;


  int*    idx_gene_name;
   
  float*  E;

   
  int     verbose = 0;
  int     k;
  char*   kmer;
   



  int         singlestrand = 0;

  int         gap = 0;
   
  int*     E_q;
  float*   E_q_bins = 0;
  

  char**   seqs;
  int      l;
   
	 
  FILE*    fpr = 0;
  


  int      max_seq   = 25000;
  int      true_nbmotifs;
  int      idx;

  int*     k_inc;
  int*     k_shu;

  int      correct;

  int      quantized = 0;
   
  int      nbrepeats = 10;



  int      ebins = 0;
  int      mbins   = 2;

  int      nk;


  int      max_nb_motifs = 10000, inc_nb_motifs=10000;
  int      add  =  2;
  int      add5 =  1;
  int      add3 =  1;


  float    max_nbsdev = 5;
  char*    outfile = 0;

  int      seed = 12345;

  int      outprm = 0;
  int      report = 1;
  char*    outrepfile;

  
  float*   seq_cond = 0;
  

  int      count_seq;
  char*    count_seq_index; 
  float    divbins     = 50.0;     
  float**  logwm;
  int      w;
  seqI     si;
  float    bkg[4] = {0.0, 0.0, 0.0, 0.0};
  int*     stars;
  char*    expanded_motif;
  float*   scores;
  int      c_rand, l_rand_src, l_rand_tgt;
  int      maxinc_max = 30;
  int      maxinc;
  int      inc;
  float    newmi;
  float    curmi;
  int      cnt_nochange;
  float    t, t_best;
  int      numshu = 0;
  int      shu_haschanged = 0;
  char**   mACE;
  int      l_m;
  int**    diwm;
  int      numrows = -1;
  int      numrows_t = -1;
  int**    intwm_0m;
  int**    intwm_1m;  
  float**  logwm_1m;
  float**  logwm_0m;

  int**    intwm;
  int      num_nochange = 10;
  int      markov = 0;
  char*    motif;
  
  int      maxnumshu = 10;
  int      corefirst = 5;

  initialize_nt(); 

  if (argc == 1) {
    printf("Usage : mi_optimize_motif_WM -motiffile FILE -expfile FILE -quantized 0/1 -fastafile FILE -rna 0/1 -outfile FILE -markov 0/1\n");
    exit(0);
  }

  expfile   = get_parameter(argc, argv, "-expfile");

  if (exist_parameter(argc, argv, "-motiffile"))     
    motiffile = get_parameter(argc, argv, "-motiffile");

  if (exist_parameter(argc, argv, "-motif"))     
    motif = get_parameter(argc, argv, "-motif");

  fastafile = get_parameter(argc, argv, "-fastafile");
  kmer      = get_parameter(argc, argv, "-kmer");
  
  if (exist_parameter(argc, argv, "-divbins")) 
    divbins = atof(get_parameter(argc, argv, "-divbins"));
   
  if (exist_parameter(argc, argv, "-nbrepeats")) 
    nbrepeats = atoi(get_parameter(argc, argv, "-nbrepeats"));

  if (exist_parameter(argc, argv, "-num_nochange")) 
    num_nochange = atoi(get_parameter(argc, argv, "-num_nochange"));

  if (exist_parameter(argc, argv, "-max_nbsdev")) 
    max_nbsdev = atoi(get_parameter(argc, argv, "-max_nbsdev"));
    
  if (exist_parameter(argc, argv, "-outfile")) 
    outfile        = get_parameter(argc, argv, "-outfile");

  if (exist_parameter(argc, argv, "-outprm")) 
    outprm        = atoi(get_parameter(argc, argv, "-outprm"));


  //
  //  controls how much to add on each side
  //
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
    
  if (exist_parameter(argc, argv, "-quantized")) 
    quantized = atoi(get_parameter(argc, argv, "-quantized"));
 
  if (exist_parameter(argc, argv, "-gap")) 
    gap            = atoi(get_parameter(argc, argv, "-gap"));
 
  if (exist_parameter(argc, argv, "-markov")) 
    markov         = atoi(get_parameter(argc, argv, "-markov"));

  if (markov == 0)
    numrows_t = 4;
  else 
    numrows_t = 16;
    
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
 
      
  //readMotifs (kmer, motiffile, max_nb_motifs, inc_nb_motifs, &motifs, &oldmis, &nbmotifs, &motifsize); 

  if (motiffile != 0)
    readKmers_general(motiffile, max_nb_motifs, inc_nb_motifs, &motifs, &nbmotifs, &motifsize); 
  else {
    motifs    = (char**) malloc (1 * sizeof(char*)); 
    motifs[0] = motif;
    nbmotifs = 1;
  }

  if (nbmotifs == 0) {
    die("Please a non-empty motif file as input\n");
  }
  
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
  if (seqI_open(&si, fastafile) == 0) {
    die("Error opening file\n");
  }

  nborfs  = 0;
  realname = (char*)calloc(100, sizeof(char));

  // used to store sequences
  seqs = (char**) malloc (max_seq * sizeof(char*));
  true_nbmotifs = 0;

  count_seq       = -1;
  count_seq_index = (char*)calloc(100000,  sizeof(char));

  while (seq = seqI_nextSequence(&si, &name, &size)) {
    
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

    // ?
    count_seq_index[ count_seq ] = 1;

    idx = (int)ep->data;       
    E[ true_idx ] = data[idx][0];
    free( data[ idx ] );

    seqs[ true_idx ] = strdup( seq );
     

    true_idx ++;
    free(seq);
    
    free(name);
    
  }

  seqI_close(&si);

  if (verbose == 1)
    printf("finished seq loading\n");


  nborfs = true_idx;

  //					       
  //  determine number of mbins and ebins
  //
  if ((quantized == 0) && (ebins == 0)) {

    // equation:  nborfs = ebins * mbins * divbins;
    // ebins = mbins = sqrt( nborfs / divbin )

    ebins = (int)(0.5 + (float)nborfs / ( divbins * mbins));
    //ebins = (int)(0.5 + sqrt( nborfs / divbins) ); //
    //mbins = ebins;

    //						
    //  quantize expression profile
    //
    //quantize_E(E, nborfs, quantized, &ebins, &E_q, &E_q_bins); 


  } else {

    //						
    //  quantize expression profile
    //

    
    //mbins = (int)(0.5 + (float)nborfs / ( divbins * ebins));
    //if (verbose == 1)
    //  printf("number of bins for motif data %d\n", mbins);
    
  }


  quantize_E(E, nborfs, quantized, &ebins, &E_q, &E_q_bins); 
  if (verbose == 1)    
    printf("number of bins for expression data %d\n", ebins);
  if (verbose == 1)
    printf("number of bins for motif data %d\n", mbins);

  
  //
  // report writing stuff
  //
  if (outfile != 0) {
    fp_outfile = fopen(outfile, "w");
    if (fp_outfile == 0) {
      printf("cannot open %s\n", outfile);
      exit(0);
    }
  } else {
    fp_outfile = stdout;
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
  //  loop over the motifs
  //
    
  for (nk=0; nk<nbmotifs; nk++) {

    l_m = strlen(motifs[nk]);
    l   = strlen(motifs[nk]) + add5 + add3;

    expanded_motif = (char*)calloc( l+1, sizeof(char));

    for (i=0; i<add5; i++)
      expanded_motif[i] = '.';
    for (i=add5, j=0; j<l_m; i++,j++)
      expanded_motif[i] = motifs[nk][j];
    for (i=add5+l_m, j=0; j<add5; i++,j++)
      expanded_motif[i] = '.';

    expanded_motif[l] = '\0';
    
    
    if (verbose == 1) {
      printf("optimizing motif %s\n", expanded_motif);
    }

    //
    // convert RE to WM
    //
    getIntegerWMfromRegexp(expanded_motif, gap, &intwm_0m, &w);
    printIntegerWM(intwm_0m, w);

    //
    //getFirstOrderMarkovIntegerWMfromRegexp(expanded_motif, gap, &intwm_1m, &w);
    //    
    

    if (markov == 1) {

      if (verbose == 1)
	printf("Convert to 1M.\n");
      
      convert_wm_0m_to_1m(intwm_0m, w, &intwm_1m);
      dint_to_scoringWM  (intwm_1m, w, &diwm); 
      integerWMtoLog(diwm, &logwm_1m, 16, w);       

      logwm = logwm_1m;
      intwm = intwm_1m;

      // printIntegerWM_1M(intwm_1m, w);    
      // printIntegerWM_1M(diwm, w);
      

    } else {
      
      integerWMtoLog(intwm_0m, &logwm_0m, 4, w);
      logwm = logwm_0m;
      intwm = intwm_0m;

    }

    // prelim
    k_inc  = (int*) malloc(w * sizeof(int));
    if (k_inc == 0) {
      printf("Cannot allocate memory for k_inc (w=%d)\n", w);
      exit(0);
    }
    
    stars = (int*)malloc( w * sizeof(int));
    for (i=0; i<w; i++)
      stars[i] = 1;
    
    // 
    // COMPUTE initial profile
    //

    if (markov == 0)
      findAllSeqsMaxWMScores(seqs, nborfs, logwm, 0, w, stars, bkg, singlestrand, &scores, 0);
    else 
      findAllSeqsMaxWMScores(seqs, nborfs, logwm, 1, w, stars, bkg, singlestrand, &scores, 0);

    if (verbose == 1)
      printf("Max score profile calculated\n");

    //									
    //  calculate initial mi value
    //
    t      = 0.0;
    t_best = 0.0;
    findMostInformativeWMThreshold(scores, E_q, nborfs, ebins, &t, &curmi);

    //float curmi_0m; 
    //float t_0m;
    //findMostInformativeWMThreshold(scores_0m, E_q, nborfs, ebins, &t_0m, &curmi_0m);

    //curmi = CalculateMIbasic(M_q, E_q, nborfs, mbins, ebins); 
    //free(M_q);
    //if (verbose == 1)
    printf("Init MI = %f\n", curmi);
    //printf("Init MI 0m = %f\n", curmi_0m);
    //
    // start optimization
    //	

    for (k=0; k<w; k++) {
      k_inc[k] = k;
    }
  
    numshu = 0;
    
    while (1) {

      // flag that determines whether change has occurred
      shu_haschanged = 0;

      //
      // create a random index
      //

      

      k_shu = shuffleInt(k_inc, w); 

      printf("Exploring columns in following order:\n"); fflush(stdout);
      for (k=0; k<w; k++) {
	printf("%3d", k_shu[k]);
      }
      printf("\n");

      for (k=0; k<w; k++) {
      
	// non-randomly pick a column
	c_rand = k_shu[k];
	printf("Column %d", c_rand);

	if (corefirst > 0) {

	  // core only
	  if ((numshu < corefirst) && ((c_rand < add5) || (c_rand >= w - add3))) {
	    printf(" (skipped)\n");
	    continue;
	  }
	  
	  // outer only
	  if ((numshu >= corefirst) && ((c_rand >= add5) || (c_rand < w - add3))) {
	    printf(" (skipped)\n");
	    continue;
	  }


	} else {
	  printf("\n");
	}

	//if (verbose == 1)
	
	
	// treat first column P(Xi) differently of course
	if ((markov == 0) || (c_rand == 0))
	  numrows = 4;
	else 
	  numrows = 16;

	//
	//  make a WM change
	//
	cnt_nochange = 0;
	while (cnt_nochange < num_nochange) {
	  
	  // randomly draw a line that is > 0.0 (substract)
	  l_rand_src = (int) ( numrows * default_rand() );
	  
	  // repeat until the entry has something to take out
	  while (intwm[c_rand][l_rand_src] <= 1) {
	    l_rand_src = (int) (numrows * default_rand() );
	  }
	  	  
	  // randomly draw a different line that is < 1.0
	  l_rand_tgt = (int) ( numrows * default_rand() );
	  while ((intwm[c_rand][l_rand_tgt] == 120) || (l_rand_tgt == l_rand_src)) {
	    l_rand_tgt = (int) ( numrows * default_rand() );
	  }
	  
	  // how much do we transfer from one line to the other one
	  //    maximum is min( l_rand_src, l_rand_tgt )
	  
	  maxinc = min(intwm[c_rand][l_rand_src], 120 - intwm[c_rand][l_rand_tgt]);	  
	  maxinc = min(maxinc, maxinc_max);	  
	  inc    = 1 + (int)(maxinc * default_rand()); 
	  
	  if (verbose == 1)
	    printf("shu = %d, col = %d, row src=%d, row tgt = %d, transfer = %d\n", numshu, c_rand, l_rand_src, l_rand_tgt, inc);

	  // modify
	  intwm[c_rand][l_rand_tgt] += inc;
	  intwm[c_rand][l_rand_src] -= inc;
	  
	  
	  // making log matrix
	  if (markov == 0) {
	    
	    integerWMtoLog(intwm, &logwm, numrows_t, w); 
	    
	  } else {
	    
	    dint_to_scoringWM(intwm, w, &diwm); 
	    integerWMtoLog(diwm, &logwm, numrows_t, w); 

	  } 
	  
	  //
	  // make score profile
	  //
	  if (markov == 0)
	    findAllSeqsMaxWMScores(seqs, nborfs, logwm, 0, w, stars, bkg, singlestrand, &scores, 0);
	  else 
	    findAllSeqsMaxWMScores(seqs, nborfs, logwm, 1, w, stars, bkg, singlestrand, &scores, 0);
	    

	  //
	  // free up logwm
	  //
	  for (i=0; i<w; i++)
	    free(logwm[i]);
	  free(logwm);
	  
	  //									
	  //  calculate new mi value
	  //
	  findMostInformativeWMThreshold(scores, E_q, nborfs, ebins, &t, &newmi);

	  free(scores);
	  
	  if (newmi > curmi) {
	    
	    printf("Increased MI to %f\n", newmi);

	    // printf("Change accepted.\n\n");
	    curmi          = newmi;
	    t_best         = t;
	    cnt_nochange   = 0;  //reset
	    shu_haschanged = 1;  // set changed flag to 1

	  } else {
	    
	    // reverse move
	    intwm[c_rand][l_rand_tgt] -= inc;
	    intwm[c_rand][l_rand_src] += inc;
	    cnt_nochange ++;

	  }

	  
	}

      }

      numshu ++;

      // if no changes, or num shu iterations is 5, exit loop
      if ((shu_haschanged == 0) ||(numshu == maxnumshu))
	break;
    }


    if (markov == 1) {
      printf("Final matrix:\n");
      printIntegerWM_1M(intwm, w);
    } else {

      integerWMtoACE(intwm, w, &mACE); 
    
      fprintf(fp_outfile, "Motif %d\n", nk);
    
      fprintf(fp_outfile, "# init = %s, MI = %f, treshold = %f\n", motifs[nk], curmi, t_best);

      for (i=0; i<120; i++) {
	for (j=0; j<w; j++) {
	  fprintf(fp_outfile, "%c", mACE[i][j]);
	}      
	fprintf(fp_outfile, "\n");
      }
      for (i=0; i<w; i++) {
	fprintf(fp_outfile, "*");
      }

      fprintf(fp_outfile, "\n");
      
      fprintf(fp_outfile, "\n");
    
      if (outfile != 0)
	fflush(fp_outfile);
      
      // free up memory for mACE
      for (i=0; i<120; i++)
	free(mACE[i]);
      free(mACE);
      
      
      free(stars);
      free(k_inc);
      
    }
    
    if (outfile != 0)
      fclose(fp_outfile);

  }

  return 0;
}


void findMostInformativeWMThreshold(float* scores, int* E_q, int nborfs, int ebins, float* t, float* mi)
{
  float* scores_copy;
  float* a_th;
  int    i;
  int    j;
  float  newmi;
  float  bestmi = -1000.0;
  int*   M_q;

  // init
  *t = 0.0;


  // copy score vector
  scores_copy = (float*)malloc( nborfs * sizeof(float));
  memcpy(scores_copy, scores, nborfs * sizeof(float));
  


  // sort scores (ascending)
  qsort((void*)scores_copy, nborfs, sizeof(float), CmpFloat);

  // alloc threshold array
  a_th = (float*)malloc(nborfs * sizeof(float));
  
  // find all thresholds
  a_th[0] = scores_copy[0];

  //printf("Scores_copy[0] = %f\n", scores_copy[0]); getchar();
  
  j = 0;
  for (i=1; i<nborfs; i++) {
    //printf("Comparing scores copy = %f and a_th[j] = %f\n", scores_copy[i], a_th[j]);
    if (scores_copy[i] > a_th[j]) {
      a_th[++j] = scores_copy[i];
    }
  }
  
  /*
  printf("Thresholds\n");
  for (i=0; i<=j; i++)
    printf("%f\n", a_th[i]);
    */

  // iterate thru the tresholds
  for (i=0; i<=j; i++) {
    // quantize accordingly

    M_q   = ThresholdQuantize(scores, nborfs, a_th[i], 1); 

    // calculate MI
    newmi = CalculateMIbasic(M_q, E_q, nborfs, 2, ebins); 

    //printf("Quantize with t = %f, mi=%f\n", a_th[i], newmi);    
    free(M_q);

    // check whether MI improves
    if (newmi > bestmi) {
      bestmi = newmi;
      *t     = a_th[i];
      *mi    = bestmi;
    }

  }


  free(scores_copy);
  free(a_th);

}
