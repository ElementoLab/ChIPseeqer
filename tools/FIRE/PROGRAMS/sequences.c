#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)

#define _GNU_SOURCE
#ifndef CPROTO
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif

#include "sequences.h"
#include "dataio.h"
#include "statistics.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif


//  P(A1,A2,A3) = P(A1)P(A2)P(A3)  (all col indep)
//  P(A1,A2,A3) = P(A1,A2)P(A3|A1,A2)   = P(A1)P(A2,A3|A1) = P(A1)
//                        P(A3|A1,A2) = 
//              = P(A1)P(A2|A1)P(A3|A2)

//  comes down to P(A3|A1,A2) = P(A3|A2) if A1 and A3 are indep, ie P(A1,A3)=0 
 
void initialize_nt() 
{
  int i, j;

  N = (int*)calloc(256, sizeof(int));
  
  N['A']  = 0;
  N['C']  = 1;
  N['T']  = 2;
  N['G']  = 3;
  N['N']  = 4;
  N['X']  = 5;

  C = (int*)calloc(256, sizeof(int));
  
  C['A']          = 2;
  C['C']          = 3;
  C['T']          = 0;
  C['G']          = 1;
  C['N']          = 4;
  C['X']          = 4;
  
  ntoC[0]         = 'A';
  ntoC[1]         = 'C';
  ntoC[2]         = 'T';
  ntoC[3]         = 'G';
  ntoC[4]         = 'N';
  ntoC[5]         = 'X';

  char_compl      = (char*)calloc(256, sizeof(char));
  char_compl['A'] = 'T';
  char_compl['T'] = 'A';
  char_compl['C'] = 'G';
  char_compl['G'] = 'C';

  N1 = (int**)calloc(256, sizeof(int*));
  for (i=0; i<256; i++) {  
    N1[i] = (int*)calloc(256, sizeof(int));
  }

  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      N1[ ntoC[i] ][ ntoC[j] ] = i*4+j;
      ntoC1[ i*4+j ][0]          = ntoC[i];
      ntoC1[ i*4+j ][1]          = ntoC[j]; 
      ntoC1[ i*4+j ][2]          = '\0';
    }
  }
  
  
}


void printWM(float** wm, int w)
{
  int i,j;
  for (j=0; j<4; j++) {
    
    for (i=0; i<w; i++) {
      printf("\t%3.2f", wm[i][j]);
    }
    printf("\n");
  }


}

//
//  read a single WM from a file (in the form of a count matrix)
//
//  transfo = 0 -> integer, 1 -> float, 2 -> log
//
void readACEintWM(char* wmfile, int*** wm, int* w, char*** sites, int* n, int** stars)
{

  FILE*   fp3;
  char*   buff;
  int     len = 100000;
  int**   mywm;
  int     cnt_sites;
  int     l;
  int     col;
  int     myw;
  int     i;
  int     star = 0;
  char**  mysites;
  int*    mystars;

  buff = (char*)malloc(1000*sizeof(char));

  //
  //  read in WM
  //
  fp3 = fopen(wmfile, "r");


  // Motif 1 line
  fgets(buff, len, fp3);

  // line 2, to get the motif width
  fgets(buff, len, fp3);
  myw = strlen(buff)-1; 
  cnt_sites = 1;

  // loop, get the number of sites
  while (!feof(fp3)) {
    //getline(&buff, &len, fp3);
    fgets(buff, len, fp3);  
    if (feof(fp3)) {
      break;
    }
    cnt_sites ++;
  }

  // start again
  rewind(fp3);
  
   
  //
  // input wm, log weights
  //
  mywm = (int**)malloc(myw * sizeof(int*));
  for (col= 0; col < myw; col++) {
    mywm[col] = (int*)malloc(5 * sizeof(int));
    for (i=0; i<5; i++) {
      mywm[col][i] = 0.0;
    }
  }

  mystars = (int*)malloc(myw * sizeof(int));

  //
  //  read in each site
  //
  mysites = (char**)malloc(cnt_sites * sizeof(char*));
  
  cnt_sites = 0;
  fgets(buff, len, fp3);

  while (!feof(fp3)) {
    fgets(buff, len, fp3);
    if (feof(fp3)) {
      break;
    }
    chomp(buff);
    l = strlen(buff);

    star = 0;
    for (i=0; i<l; i++) {
      if ((buff[i] == ' ') || (buff[i] == '*')) {
	star = 1;  break; 
      }
    }

    if (star == 1) {
      for (i=0; i<l; i++) {
	if (buff[i] == '*') 
	  mystars[i] = 1;
	else
	  mystars[i] = 0;
	  
      }
    
      break;
    }

    for (i=0; i<l; i++) {

      mywm[i][ N[ (int)(buff[i]) ] ] += 1;
    }
    
    mysites[ cnt_sites ] = strdup(buff);
    
    cnt_sites++;

  }
  
  //
  //  SKIP: normalize by the number of sites
  //

  /*
  for (col= 0; col < myw; col++) {
      for (i=0; i<4; i++) {	
	//+= bkg[i]; //epsilon;
	int ini = mywm[col][i];
	mywm[col][i]  = log2( ((mywm[col][i] + bkg[i]) / (float)(cnt_sites + 1)) / bkg[i] ); 

      }
  }
  */
  
  //
  //  output 
  //
  *w     = myw;
  *wm    = mywm;
  *sites = mysites;
  *n     = cnt_sites;
  *stars = mystars;

}


//
// 
//
void ACEintWMtologWM(int** intwm, int w, int nsites, float* rawbkg, float*** logwm) {

  int col, i;
  float** mywm;

  mywm = (float**)malloc(w * sizeof(float*));
  for (col= 0; col<w; col++) {
    mywm[col] = (float*)malloc(5 * sizeof(float));
    for (i=0; i<5; i++) {
      mywm[col][i] = 0.0;
    }
  }

  for (col= 0; col < w; col++) {
    for (i=0; i<4; i++) {	
      mywm[col][i]  = log2( ((intwm[col][i] + rawbkg[i]) / (float)(nsites + 1)) / rawbkg[i] ); 
    }
  }
  
  *logwm = mywm;
  
}



//
//  log score
//
float getScoreOnSeq(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn) 
{
  float score1 = 0.0;
  float score2 = 0.0;
  
  int  j, k;

  *hasn = 0;

  for (j=0; j<w; j++) {
    if (stars[j] == 1) {
      if ((seq[i + j] != 'A') && (seq[i + j] != 'T') && (seq[i + j] != 'C') && (seq[i + j] != 'G')) {
	*hasn = 1;
	return 0.0;
      }

      score1 += wml[j][ N[ (int)(seq[i + j]) ] ]; // - bkg[ N[ (int)(seq[i + j]) ] ];
    }
  }

  if (rna == 0) {  // if not RNA, also score reverse complement

    for (j=w-1,k=0; j>=0; j--,k++) {
      if (stars[k] == 1) {
	//printf(" += %f\n", wml[k][ C[ (int)(seq[i + j]) ] ]);
	score2 += wml[k][ C[ (int)(seq[i + j]) ] ]; // - bkg[ C[ (int)(seq[i + j]) ] ];
      }
    }

  }

  //printf("score1 0M = %f, score2 0M = %f\n", score1, score2);

  
  if ((rna == 1) || (score1 >= score2)) {
    *strand = 1;
    return score1;
  } else {
    *strand = -1;
    return score2;
  }
  

}


void findAllWeightMatrixMatches(float** wml, int w, int* stars, float* bkgl, float t, char* seq, int rna, int** matches_pos, int* num_matches, char** matches_ori, int max_num_matches)
{

  int i;
  int l = strlen(seq);
  float s;
  int strand;
  int my_num_matches = 0;
  int my_max_num_matches = max_num_matches;
  void* ptr;
  int  hasn;

  *matches_pos = (int*) malloc(my_max_num_matches * sizeof(int));
  if (!*matches_pos) {
    die("findAllWeightMatrixMatches: not enough memory for matches_pos\n");    
  }
  //*matches_ori = (char*)malloc(my_max_num_matches * sizeof(char));
  //if (!*matches_ori) {
  //  die("findAllWeightMatrixMatches: not enough memory for matches_ori\n");    
  //}



  for (i=0; i<l; i++) {
    s = getScoreOnSeq(wml, w, stars, bkgl, seq, i, rna, &strand, &hasn);
    
    if (!hasn && (s >= t)) {

      /*
      if (strand > 0) 
	printf("%s => %4.3f (%d)\n", substr(seq, i, w), s, strand); 
      else 
	printf("%s => %4.3f (%d)\n", complement(substr(seq, i, w)), s, strand); 
	//printf("%s => %4.3f (%d)\n", , s, strand); 
	*/

      (*matches_pos)[ my_num_matches ] = i;
      //(*matches_ori)[ my_num_matches ] = (char)strand;
      my_num_matches ++;

      if (my_num_matches == my_max_num_matches) {
	my_max_num_matches += max_num_matches;
	ptr = realloc( *matches_pos, my_max_num_matches * sizeof(int));	
	if (ptr == 0)
	  die("findAllWeightMatrixMatches: not enough memory for matches_pos (realloc)\n");
	//ptr = realloc( *matches_ori, my_max_num_matches * sizeof(char));		
	//if (ptr == 0)
	//  die("findAllWeightMatrixMatches: not enough memory for matches_ori (realloc)\n");
	  
      }
      
    }    
  } 

  *num_matches = my_num_matches;
  
}


//
// this function scans all sequences and returns the max WM scores for each
//
void findAllSeqsMaxWMScores(char** seqs, int numseqs, float** logwm, int mt, int w, int* stars, float* bkg, int rna, float** scores, int** a_idx)
{
  int   i;
  float* myscores;
  int   idx;

  int*  my_a_idx;
  
  myscores = (float*)malloc( numseqs* sizeof(float));
  if (myscores == 0) {
    die("findAllSeqsMaxWMScores: cannot allocate memory for myscores\n");
  }
  
  if (a_idx != 0)
    my_a_idx = (int*)malloc( numseqs * sizeof( int ));


  for (i=0; i<numseqs; i++) {    
    myscores[i] = findMaxWMScore(seqs[i], logwm, mt, w, stars, bkg, rna, &idx);
    if (a_idx != 0)
      my_a_idx[i] = idx;
    //printf("%f %s\n", myscores[i], substr(seqs[i], idx, w));
  }
  
  *scores  = myscores;
  
  if (a_idx != 0)
    *a_idx   = my_a_idx;
}


float findMaxWMScore(char* seq, float** logwm, int mt, int w, int* stars, float* bkg, int rna, int* idx)
{

  int i, l;
  float maxs = -1000000.0;
  float s;
  int   hasn, strand;

  l = strlen(seq);
  *idx = -1;

  for (i=0; i<l-w; i++) {    

    if (mt == 0) 
      s = getScoreOnSeq(logwm, w, stars, bkg, seq, i, rna, &strand, &hasn);	    
    else if (mt == 1)
      s = getScoreOnSeq_1M(logwm, w, stars, bkg, seq, i, rna, &strand, &hasn);	 

    if ((hasn == 0) && (s > maxs)) {      
      maxs = s;
      *idx = i;
    }      
  }

  return maxs;
}



//
//  get the score threshold out of a WM (log or not)
//
float getWMThreshold(float** wm, int w, int* stars, float* bkg, int transfo, char** sites, int n) {
  
  float score; 
  int   i, j;
  //char* T    = "ACTG";
  float SX   = 0.0;
  float SX2  = 0.0;
  float AVG;
  float STD  = 0.0;
  float t;
  float* scores;
  int   cnt_above;

  scores = (float*)malloc(n * sizeof(float));
  
  for (i=0; i<n; i++) {

    //printf("site = %s\n", sites[i]);
    
    if (transfo == 0) {
      score = 1.0;
    } else {
      score = 0.0;
    }

    for (j=0; j<w; j++) {
      if (stars[j] == 0)
      	continue;
      if (transfo == 0) {
	//printf(" %c mult pby %4.3f / %4.3f\n",sites[i][j], wm[j][ N[sites[i][j]] ],
	//     bkg[ N[ sites[i][j] ] ]);

	score *= wm[j][ N[  (int)(sites[i][j]) ] ] / bkg[ N[ (int)(sites[i][j]) ] ];
	//printf(" score = %4.3f\n", score);
      } else {
	score += wm[j][ N[ (int)(sites[i][j]) ] ]; // - bkg[ N[ (int)(sites[i][j]) ] ];
	//printf("j=%d %c mult pby %4.3f (%c) - %4.3f\n", j, sites[i][j], wm[j][ N[sites[i][j]] ], sites[i][j], bkg[ N[ sites[i][j]] ]);
	
      }
    }

    ///printf("%s\t%4.3f\n", sites[i], score);
    
    SX  += score;
    SX2 += score * score;

    scores[i] = score;

  }
  
  AVG =  SX/n;
  STD =  (float)sqrt(( SX2 - SX*SX / n )  / (n - 1));
  //printf("AVG = %3.2f, S = %3.2f\n", AVG, STD);

  t = AVG - 2.0 * STD;

  cnt_above = 0;
  for (i=0; i<n; i++) {
    if (scores[i] >= t) {
      cnt_above ++;
    }
  }

  //printf("%d/%d sites score better than t\n", cnt_above, n);

  return t;

}




char* getGappedKmer(char* kmer, int gap) 
{
  
  char*    stmp;
  int      kmersize;
  int      j;
  int      half1, half2;

  
  kmersize = strlen(kmer);
  half1 = kmersize/2;
  half2 = kmersize-half1;
  
  stmp = (char*)calloc((kmersize+gap+1), sizeof(char));
 
  strncpy(stmp, kmer, half1);  
  for (j=0; j<gap; j++) {
    stmp[half1+j] = '.';
  }
  strncpy(stmp+half1+gap, kmer + half1, half2);
  stmp[kmersize+gap] = '\0';
  

  return stmp;
}


char* getGappedMotif(char* motif, int motifsize, int gap) 
{
  
  char*    newmotif;
  int      half1;
  int      l = strlen(motif);
  int      i = 0;  // cnt in ungapped motif
  int      j = 0;  // cnt in   gapped motif
  int      k = 0;  // cnt characters
  int      h = 0;
  int      p_b;

  newmotif = (char*)calloc((l+gap+1), sizeof(char)); 

  half1 = motifsize/2;

  while (i < l) {

    if (motif[i] == '[') {

      p_b = i;
      while (motif[p_b] != ']') {	
	p_b++; 
      } 

      while (i <= p_b) {
	newmotif[j] = motif[i];
	i++; j++;
      }

      k++;

    } else {

      newmotif[j] = motif[i];
      j++; i++; k++;

    }
    
    if (k == half1) {
      for (h=0; h<gap; h++) {
	newmotif[j] = '.';
	j++;
      }
    }

  }

  newmotif[l+gap] = '\0';
  
  return newmotif;
}



void getIntegerWMfromRegexp(char* motif, int gap, int*** wm, int* w)
{
  int   myw;
  int** mywm;
  
  int      l = strlen(motif);
  int      i = 0;  // cnt in ungapped motif
  int      j = 0;  // cnt in   gapped motif
  int      k = 0;  // cnt characters
  int      p_b;
  int      nl;
  int      h;

  myw = getRegexpMotifLength(motif, gap);

  // allocate memory
  mywm = (int**)calloc( myw, sizeof(int*) );
  for (i=0; i<myw; i++) {
    mywm[i] = (int*)calloc(4, sizeof(int));
  }

  // scan the regexp
  i = 0;

  while (i < l) {

    if (motif[i] == '[') {

      p_b = i;
      while (motif[p_b] != ']') {	
	p_b++; 
      } 

      nl = p_b - i - 1;

      if ((nl != 2) && (nl != 3)) {
	printf("Problem: nl should be 2 or 3, it is %d\n", nl);
	exit(0);
      }
      
      nl = 120 / nl;

      for (h=i+1; h<p_b; h++) {
	
	mywm[k][ N[ (int)(motif[h]) ] ] = nl;
	
      }

      while (i <= p_b) {	
	i++; j++;
      }

      k++;

    } else {
      
      if (motif[i] == '.') {
	for (h=0; h<4; h++) {
	  mywm[k][h] = 120 / 4;
	}		
      } else if ((motif[i] == 'A') ||
		 (motif[i] == 'T') ||
		 (motif[i] == 'C') ||
		 (motif[i] == 'G')) {
	mywm[k][ N[ (int)(motif[i]) ] ] = 120;
      }

      j++; i++; k++;
    }    
  }  

  *w  = myw;
  *wm = mywm;
}



void getLinearFloatWMfromRegexp(char* motif, float** lwm, int* lw)
{
  int   w;
  int** wm;
  int   sum;
  int   i, j;

  getIntegerWMfromRegexp(motif, 0, &wm, &w);
  
  *lwm = (float*)calloc(w * 4, sizeof(float));
  
  for (i=0; i<w; i++) {
    sum = 0;
    for (j=0; j<4; j++) {
      sum += wm[i][j];
    }
    for (j=0; j<4; j++) {
      (*lwm)[i*4+j] = wm[i][j] / sum;
    }
  }

  *lw = w;
  
  free(wm);
}

  


//
//  get a P(A)P(B|A)P(C|B) ... matrix
//
void getFirstOrderMarkovIntegerWMfromRegexp(char* motif, int gap, int*** wm, int* w)
{
  int      myw;
  int**    mywm;  
  int      l = strlen(motif);
  int      i = 0;  // cnt in ungapped motif
  int      j = 0;  // cnt in   gapped motif
  int      k = 0;  // cnt characters
  int      m = 0;
  int      p_b;
  int      nl;
  int      h;
  char**   a_c;
  int      nl0, nl1;

  myw = getRegexpMotifLength(motif, gap);

  // allocate memory
  mywm = (int**)calloc( myw, sizeof(int*) );
  for (i=0; i<myw; i++) {
    mywm[i] = (int*)calloc(16, sizeof(int));
  }

  // array of characters
  a_c = (char**)calloc(myw, sizeof(char*));

  // scan the regexp
  i = 0;

  while (i < l) {

    if (motif[i] == '[') {
      
      p_b = i;
      while (motif[p_b] != ']') {	
	p_b++; 
      } 

      nl = p_b - i - 1;

      if ((nl != 2) && (nl != 3)) {
	printf("Problem: nl should be 2 or 3, it is %d\n", nl);
	exit(0);
      }
      
      a_c[k] = (char*)calloc(nl+1, sizeof(char));


      // fill up with all allowed characters at that position
      for (h=i+1,m=0; h<p_b; h++,m++) {	
	a_c[k][m] = motif[h];	
      }
      a_c[k][nl] = '\0';

      while (i <= p_b) {	
	i++; j++;
      }

      k++;

    } else {
      
      if (motif[i] == '.') {

	// . characters
	a_c[k] = (char*)calloc(5, sizeof(char));	      
	for (h=0; h<4; h++) {
	  a_c[k][h] = ntoC[h];
	}		
	a_c[k][4] = '\0';


      } else if ((motif[i] == 'A') ||
		 (motif[i] == 'T') ||
		 (motif[i] == 'C') ||
		 (motif[i] == 'G')) {
	
	a_c[k] = (char*)calloc(2, sizeof(char));
	a_c[k][0] = motif[i];
	a_c[k][1] = '\0';

      }

      j++; i++; k++;
    }    
  }  

  //for (i=0; i<myw; i++) {
  //  printf("a_c[%d] = %s\n", i, a_c[i]);
  //}

  
  // now use the a_c to build the weight matrix

  // first column, single nt   
  nl = strlen(a_c[0]);

  for (j=0; j<16; j++)
    mywm[0][j] = 0;

  for (i=0; i<nl; i++)
    mywm[0][ N[ (int)(a_c[0][i]) ] ] = 120 / nl;
  
  for (i=1; i<myw; i++) {
    nl0 = strlen(a_c[i-1]);
    nl1 = strlen(a_c[i  ]);

    //printf("nl0 = %s, nl1 = %s\n", a_c[i-1], a_c[i]);

    // init column
    for (j=0; j<16; j++)
      mywm[i][j] = 0;

    // 
    for (j=0; j<nl0; j++) {  // the given .. P(.|J)
      for (k=0; k<nl1; k++) {  // P(K|.)
	//printf("%c %c\n", a_c[i-1][j], a_c[i][k]);
	//printf("N1 = %d\n", N1[ (int)(a_c[i-1][j]) ][ (int)(a_c[i][k]) ]);
	// mywm[i][ N1[ (int)(a_c[i-1][j]) ][ (int)(a_c[i][k]) ] ] = (int)(120.0 / (nl0*nl1));

	// records P(Xi,Xi-1)
	mywm[i][ N1[ (int)(a_c[i][k]) ][ (int)(a_c[i-1][j]) ] ] = (int)(120.0 / (nl0*nl1));
	
      }
    }

  }


  //printIntegerWM(mywm, myw);

  *w  = myw;
  *wm = mywm;
}


//
// convert a 0th-order Markov WM into a first-order Markov model
//  this is done by making all previous nt equally probable (dividing by 4 basically)
//
void convert_wm_0m_to_1m(int** intwm_0m, int w, int*** intwm_1m)
{
  
  int lig;
  int** mywm2;
  int i, j;
  int dint_idx;
  

  // allocate memory
  mywm2 = (int**)calloc( w, sizeof(int*) );
  for (i=0; i<w; i++) {
    mywm2[i] = (int*)calloc(16, sizeof(int));
  }

  // 0 is unchanged
  for(lig=0; lig < 4; lig++) {
    mywm2[0][lig] = intwm_0m[0][lig];
  }  

  // 1+
  for (i=1; i<w; i++) {
    for(lig=0; lig < 4; lig++) {      
      for (j=0; j<4; j++) {	
	dint_idx = N1[ ntoC[lig] ][ ntoC[j] ];
	
	//printf("i=%d, lig=%d, j=%d, dint_idx=%d, intwm_0m[i][lig]=%d\n", i, lig, j, dint_idx, intwm_0m[i][lig]);
	
	mywm2[ i ][ dint_idx ] = intwm_0m[i][lig] / 4;
      }
    }
  }
  
  *intwm_1m = mywm2;
}


// transform a dinucleotide count vector into a P(Xi|Xi-1) matrix to be used for scoring.
void dint_to_scoringWM(int** diwm, int w, int*** scwm) 
{
  
  int lig, col;
  int** mywm2;
  int i, j;
  int P[4];
  int dint_idx;
  

  // allocate memory
  mywm2 = (int**)calloc( w, sizeof(int*) );
  for (i=0; i<w; i++) {
    mywm2[i] = (int*)calloc(16, sizeof(int));
  }

  // 0 is unchanged
  for(lig=0; lig < 4; lig++) {
    mywm2[0][lig] = diwm[0][lig];
  }  


  for (col= 1; col < w; col++) {
    
    // calculate P(Xi-1) by summing P(.,Xi-1)
    for (i=0; i<4; i++) {
      P[i] = 0;
      for (j=0; j<4; j++) {
	P[i] += diwm[ col ][ N1[ ntoC[j] ][ ntoC[i] ] ]; 
      }      

      //printf("col=%d, P[i=%d] = %d\n", col, i, P[i]);

    }
    
    // divide P(Xi,Xi-1)/P(Xi-1)
    for (i=0; i<4; i++) {
      for (j=0; j<4; j++) {		
	dint_idx = N1[ ntoC[j] ][ ntoC[i] ];

	if (diwm[ col ][ dint_idx ] == 0) 
	  mywm2[col][ dint_idx ] = 0;
	else
	  mywm2[col][ dint_idx ] = (int)( 120 * diwm[ col ][ dint_idx ] /  (float)P[i] ); 
      }      
    }
    
  }
  
  *scwm = mywm2;
  

}



//
//  log score
//
float getScoreOnSeq_1M(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn) 
{
  float score1 = 0.0;
  float score2 = 0.0;
  
  int   j, k;
  int   nt_i, nt_im1, dint_idx;

  *hasn = 0;

  // exit if sequence to score has Ns
  for (j=0; j<w; j++) {
    //if (stars[j] == 1) {
    if ((seq[i + j] != 'A') && (seq[i + j] != 'T') && (seq[i + j] != 'C') && (seq[i + j] != 'G')) {
      *hasn = 1;
      return 0.0;
    }      
  }
  //}

  
  // score first position
  score1 = wml[0][ N[ (int)(seq[i]) ] ];
  //printf(" += %f\n", wml[0][ N[ (int)(seq[i]) ] ]);
  if (rna == 0) {
    
    //printf(" += %f\n", wml[0][ C[ (int)(seq[i+w-1]) ] ]);    
    score2 = wml[0][ C[ (int)(seq[i+w-1]) ] ];

  }

  // score rest
  for (j=1,k=w-2; j<w; j++,k--) {

    nt_i     = seq[ i+j   ];  // nt at i 
    nt_im1   = seq[ i+j-1 ];  // nt at i-1
    dint_idx = N1[nt_i][nt_im1];

    //printf(" += %f (j=%d) wml[%d][%d] for %c%c\n", wml[j][dint_idx], j, j, dint_idx, seq[ i+j   ], seq[ i+j-1   ]);
    score1  += wml[j][dint_idx];

    if (rna == 0) {
      
      nt_i     = char_compl[ (int)(seq[ i+k   ]) ];  // nt at i 
      nt_im1   = char_compl[ (int)(seq[ i+k+1 ]) ];  // nt at i-1
      dint_idx = N1[nt_i][nt_im1];
      score2  += wml[j][dint_idx];
      
      //printf(" += %f (for nti=%c|ntim1=%c) dint_idx=%d  \n", wml[j][dint_idx], seq[i+k], seq[i+k+1], dint_idx);    

      
    }

  }
  
  //printf("score1 1M = %f, score2 1M = %f\n", score1, score2);

  if ((rna == 1) || (score1 >= score2)) {
    *strand = 1;
    return score1;
  } else {
    *strand = -1;
    return score2;
  }
  

}


void printIntegerWM(int** wm, int w)
{
  int i,j;
  for (j=0; j<4; j++) {
    
    for (i=0; i<w; i++) {
      printf("%4d", wm[i][j]);
    }
    printf("\n");
  }


}

void printIntegerWM_1M(int** wm, int w)
{
  int i,j;
  for (j=0; j<16; j++) {
    
    printf("%4s\t", ntoC1[j]);
    for (i=0; i<w; i++) {
      printf("%4d", wm[i][j]);
    }
    printf("\n");
  }


}


void freeIntegerWM(int** wm, int w) {
  
  int i;
  for(i=0; i<w; i++) {
    free(wm[i]);
  }
  free(wm);
}



void integerWMtoACE(int** intwm, int w, char*** m) 
{
  int i, j, k, l;
  char** mym;

  mym = (char**)malloc(120 * sizeof(char*));
  for (i=0; i<120; i++)
    mym[i] = (char*)malloc(w * sizeof(char));

  for (i=0; i<w; i++) {
    
    l = 0;
    for (j=0; j<4; j++) {      
      for (k=0; k<intwm[i][j]; k++) {	
	mym[l][i] = ntoC[ j ];
	l++;
      }
    }
    
  }

  *m = mym;
}


//  wm2 = log(wm1)
void integerWMtoLog(int** wm1, float*** wm2, int h, int w) 
{
  int lig, col;
  float pseudo = 0.000000001;
  float** mywm2;
  int i;

  // allocate memory
  mywm2 = (float**)calloc( w, sizeof(float*) );
  for (i=0; i<w; i++) {
    mywm2[i] = (float*)calloc(h, sizeof(float));
  }



  for(lig=0; lig < h; lig++) {
    for(col= 0; col < w; col++) {
      mywm2[col][lig] = log( (wm1[col][lig] / 120.0) + pseudo );
    }
  }
  
  *wm2 = mywm2;
  
}

int getRegexpMotifLength(char* motif, int gap) 
{
  

  int      l = strlen(motif);
  int      i = 0;  // cnt in ungapped motif
  int      j = 0;  // cnt in   gapped motif
  int      k = 0;  // cnt characters
  //int      h = 0;
  int      p_b;


  while (i < l) {

    if (motif[i] == '[') {

      p_b = i;
      while (motif[p_b] != ']') {	
	p_b++; 
      } 

      while (i <= p_b) {	
	i++; j++;
      }

      k++;
    } else {
      j++; i++; k++;
    }    
  }
  
  return k + gap;
}

    
float compare_motifs(float *m1, int l1, float* m2, int l2, int ss, int min_inf_cols) 
{

  //int    min_inf_cols = 6;
  int    verbose      = 0;
  int   *st1, *st2, *st3;
  float* in1;
  float* in2;
  int*   in1_idx;
  int*   in2_idx;
  int    cnt_st1, cnt_st2;
  float  lg2;
  int    i, j;
  float  *m3;
  int    i1, i2;
  float  pe, pe_max;
  int    cnt;
  float* V1, *V2;
  int    start_i1, start_i2;
  float  min_inf_t    = 0.8;

  float  last_i;
  float  eps = 1e-10;
  lg2 = log(2.0);

  // minimum number of informative columns is of course bounded by motif width
  min_inf_cols = min(min_inf_cols, l1);
  min_inf_cols = min(min_inf_cols, l2);
  
  st1 = (int*)calloc( l1, sizeof(int));  
  st2 = (int*)calloc( l2, sizeof(int));

  if (verbose == 1)
    printf("min_inf_cols = %d\n", min_inf_cols);

  // calculate information content of each column of motif 1
  in1 = (float*)calloc(l1, sizeof(float));
  for (i=0; i<l1; i++) {
    in1[i] = 2.0;
    for (j=0; j<4; j++) {
      if (fabs(m1[i*4+j]) > eps)
	in1[i] += m1[i*4+j] * log( m1[i*4+j] ) / lg2;
    }
    //printf("pos = %d, I=%4.1f\n", i, in1[i]);
  }


  // flag the 6+ most informative columns
  in1_idx = bubbleSortIndex(in1, l1);
  for (i=0; i<min_inf_cols; i++) {
    //printf("I = %4.3f (%d, 1)\n", in1[ in1_idx[i] ], in1_idx[i]);

    st1[ in1_idx[i] ] = 1;
    last_i            = in1[ in1_idx[i] ];

  }
  i=min_inf_cols;
  //while ((i<l1) && ( ( in1[ in1_idx[min_inf_cols-1] ] - in1[ in1_idx[i] ] ) < 0.001 )) {
  while ((i<l1) && ( in1[ in1_idx[i] ] >= min_inf_t)) {

    //printf("5 - %d = %4.3f\n", i, ( in1[ in1_idx[5] ] - in1[ in1_idx[i] ] ));    
    //printf("I = %4.3f (%d, 1)\n", in1[ in1_idx[i] ], in1_idx[i]);
    st1[ in1_idx[i] ] = 1;
    last_i            = in1[ in1_idx[i] ];

    i++;
  }
  for (j=i; j<l1; j++) {
    //printf("I = %4.3f (%d, 0)\n", in1[ in1_idx[j] ], in1_idx[j]);
    if ( (last_i - in1[ in1_idx[i] ]) < 0.001)
      st1[ in1_idx[j] ] = 1;
    else
      st1[ in1_idx[j] ] = 0;
  }
  
  /*
  for (i=0; i<l1; i++) {
    printf(" %4d", st1[i]);
  }
  printf("\n\n");
  */

  // same for motif 2
  in2 = (float*)calloc(l2, sizeof(float));
  for (i=0; i<l2; i++) {
    in2[i] = 2.0;
    for (j=0; j<4; j++) {
      if (fabs(m2[i*4+j]) > eps)
	in2[i] += m2[i*4+j] * log( m2[i*4+j] ) / lg2;
    }
    //printf("pos = %d, I=%4.1f\n", i, in2[i]);
  }

  // flag the 6+ most informative columns
  in2_idx = bubbleSortIndex(in2, l2);
  for (i=0; i<min_inf_cols; i++) {
    //printf("I = %4.3f (%d, 1)\n", in2[ in2_idx[i] ], in2_idx[i]);    
    st2[ in2_idx[i] ] = 1;
    last_i            = in2[ in2_idx[i] ];

  }

  // if (m
  i=min_inf_cols;
  while ((i<l2) && ( in2[ in2_idx[i] ] >= min_inf_t)) {

  //while ((i<l2) && ( ( in2[ in2_idx[min_inf_cols-1] ] - in2[ in2_idx[i] ] ) < 0.001 )) {    
    //printf("5 - %d = %4.3f\n", i, ( in2[ in2_idx[5] ] - in2[ in2_idx[i] ] ));    
    //printf("I = %4.3f (%d, 1)\n", in2[ in2_idx[i] ], in2_idx[i]);

    st2[ in2_idx[i] ] = 1;
    last_i            = in2[ in2_idx[i] ];

    i++;
  }

  for (j=i; j<l2; j++) {
    //printf("I = %4.3f (%d, 0)\n", in2[ in2_idx[j] ], in2_idx[j]);	
    //st2[ in2_idx[j] ] = 0;
    
    if ( (last_i - in2[ in2_idx[i] ]) < 0.001)
      st2[ in2_idx[j] ] = 1;
    else
      st2[ in2_idx[j] ] = 0;

  }
  
  /*
  for (i=0; i<l2; i++) {
    printf(" %4d", st2[i]);
  }
  printf("\n\n");
  */

  //for (i=0; i<l2; i++) {
  //  printf(" %3.2f", in2[ in2_idx[i] ]);
  //}
    
  if (verbose == 1) {

    printf("l1=%d\n", l1);
    for (i=0; i<l1; i++) {
      for (j=0; j<4; j++) {
	printf(" %4.1f", m1[i*4+j]);
      }
      printf("\t%3.2f\t%4d\n", in1[i], st1[i]);
    }
    
    printf("l2=%d\n", l2);
    for (i=0; i<l2; i++) {
      for (j=0; j<4; j++) {
	printf(" %5.2f", m2[i*4+j]);
      }
      printf(" %5d\n", st2[i]);
    }

  }

  
  V1 = (float*)malloc(4*l1 * sizeof(float));
  V2 = (float*)malloc(4*l2 * sizeof(float));

  //
  // m1 |    x    m2 -
  //
  // diago starting on m1
  pe_max = 0.0;
  for (start_i1=0; start_i1<l1; start_i1++) {

    i1      = start_i1;
    i2      = 0;
    cnt     = 0;						
    cnt_st1 = 0;
    cnt_st2 = 0;

    if (verbose == 1)	
      printf("Aln starting at i1=%d, i2=%d\n", i1, i2);

    while ((i1 < l1) && (i2 < l2)) {
      
      //printf("%d-%d ", i1, i2);
      if (st1[i1] == 1)
	cnt_st1 ++;
      if (st2[i2] == 1)
	cnt_st2 ++;
	
      /*
      printf("%d-%d\t", i1, i2);      
      for (j=0; j<4; j++) {
	printf(" %4.2f", m1[i1*4+j]);
      }
      printf(" and ");
      for (j=0; j<4; j++) {
	printf(" %4.2f", m2[i2*4+j]);
      }   
      printf("\n");
      */

      for (j=0; j<4; j++) {
	V1[cnt] = m1[i1*4+j];
	V2[cnt] = m2[i2*4+j];      
	cnt ++;      
      }
      
      //}

      i1 ++;
      i2 ++;
      
    }
    //printf(", cnt=%d\n", cnt/4);

    if (verbose == 1)
      printf("cnt_st1=%d and cnt_st2=%d\n", cnt_st1, cnt_st2);
	
    
    if ((cnt_st1 >= min_inf_cols) && (cnt_st2 >= min_inf_cols)) {
      
      pe = pearson(V1, V2, cnt);
      
      
      /*
      for (i=0; i<cnt; i++) {
	printf("%f\t%f\n", V1[i], V2[i]);
      }
      */

      if (verbose == 1)
	printf("pe (n=%d) = %4.3f\n", cnt, pe);
      
      if (pe > pe_max) {
	pe_max = pe;
      }
      
    }	
  }

  
  // diago starting on m2
  for (start_i2=0; start_i2<l2; start_i2++) {
    
    i1  = 0;
    i2  = start_i2;
    cnt = 0;
    cnt_st1 = 0;
    cnt_st2 = 0;

    if (verbose == 1)
      printf("Aln starting at i1=%d, i2=%d\n", i1, i2);
  
    while ((i1 < l1) && (i2 < l2)) {
      

      if (st1[i1] == 1)
	cnt_st1 ++;
      if (st2[i2] == 1)
	cnt_st2 ++;

      /*
      printf("%d-%d\t", i1, i2);      
      for (j=0; j<4; j++) {
	printf(" %4.2f", m1[i1*4+j]);
      }
      printf(" and ");
      for (j=0; j<4; j++) {
	printf(" %4.2f", m2[i2*4+j]);
      }   
      printf("\n");
      */
 
      for (j=0; j<4; j++) {
	V1[cnt] = m1[i1*4+j];
	V2[cnt] = m2[i2*4+j];      
	cnt ++;      
      }
      
      i1 ++;
      i2 ++;
      

    }

    //}
    
    //printf(", cnt=%d\n", cnt/4);

    //if ((cnt / 4) >= 6) {

    if (verbose == 1)
      printf("cnt_st1=%d and cnt_st2=%d\n", cnt_st1, cnt_st2);


    if ((cnt_st1 >= min_inf_cols) && (cnt_st2 >= min_inf_cols)) {

      pe = pearson(V1, V2, cnt);

      if (verbose == 1)
	printf("pe = %4.3f\n", pe);
      
      if (pe > pe_max) {
	pe_max = pe;
      }
    
    }	
  }

  if (ss == 0) {

    // take rev comp of motif 2
    m3  = (float*)malloc(l2*4 * sizeof(float));
    st3 = (int*)  malloc(l2 * sizeof(int  ));
    
    for (i=l2-1,j=0; i>=0; i--,j++) {
      
      //printf("i=%d, j=%d\n", i, j);
      
      // A == T
      m3[ j*4 + 0 ] = m2[ i*4 + 3 ];
    
      // T == A
      m3[ j*4 + 3 ] = m2[ i*4 + 0 ];
      
      // C == G
      m3[ j*4 + 1 ] = m2[ i*4 + 2 ];
      
      // G == C
      m3[ j*4 + 2 ] = m2[ i*4 + 1 ];
      
      st3[ j ] = st2[ i ];
      
    }


      if (verbose == 1) {

	printf("l3=%d\n", l2);
	for (i=0; i<l2; i++) {
	  for (j=0; j<4; j++) {
	    printf(" %4.1f", m3[i*4+j]);
	  }
	  printf(" %4d\n", st3[i]);
	}
      }
    
    
  
    
    for (start_i1=0; start_i1<l1; start_i1++) {
      
      i1  = start_i1;
      i2  = 0;
      cnt = 0;
      cnt_st1 = 0;
      cnt_st2 = 0;

      if (verbose == 1)
	printf("Aln starting at i1=%d, i2=%d\n", i1, i2);
      
      
      while ((i1 < l1) && (i2 < l2)) {
      
	//printf("%d-%d ", i1, i2);            
	//if ((st1[i1] == 1) && (st3[i2] == 1)) {

	if (st1[i1] == 1)
	  cnt_st1 ++;
	if (st3[i2] == 1)
	  cnt_st2 ++;
	
	for (j=0; j<4; j++) {
	  V1[cnt] = m1[i1*4+j];
	  V2[cnt] = m3[i2*4+j];      
	  cnt ++;      
	}
      
	//}

	i1 ++;
	i2 ++;
      
      }

      // printf(", cnt=%d\n", cnt/4);
    
      // if ((cnt / 4) >= 6) {

      if (verbose == 1)
	printf("cnt_st1=%d and cnt_st2=%d\n", cnt_st1, cnt_st2);


      if ((cnt_st1 >= min_inf_cols) && (cnt_st2 >= min_inf_cols)) {

	pe = pearson(V1, V2, cnt);
      
	if (verbose == 1)	
	  printf("pe = %4.3f\n", pe);
      
	if (pe > pe_max) {
	  pe_max = pe;
	}
      
      }	
    }

  
    // diago starting on m2
    for (start_i2=0; start_i2<l2; start_i2++) {
    
      i1  = 0;
      i2  = start_i2;
      cnt = 0;
      cnt_st1 = 0;
      cnt_st2 = 0;


      if (verbose == 1)
	printf("Aln starting at i1=%d, i2=%d\n", i1, i2);

  
      while ((i1 < l1) && (i2 < l2)) {
      
	// printf("%d-%d ", i1, i2);      
	// if ((st1[i1] == 1) && (st3[i2] == 1)) {

	if (st1[i1] == 1)
	  cnt_st1 ++;
	if (st3[i2] == 1)
	  cnt_st2 ++;

	for (j=0; j<4; j++) {
	  V1[cnt] = m1[i1*4+j];
	  V2[cnt] = m3[i2*4+j];      
	  cnt ++;      
	}
      
	//}
      
	i1 ++;
	i2 ++;
      
      }
    
      // printf(", cnt=%d\n", cnt/4);

      //if ((cnt / 4) >= 6) {

      if (verbose == 1)
	printf("cnt_st1=%d and cnt_st2=%d\n", cnt_st1, cnt_st2);
    
      if ((cnt_st1 >= min_inf_cols) && (cnt_st2 >= min_inf_cols)) {

	pe = pearson(V1, V2, cnt);

      
	if (verbose == 1)
	  printf("pe = %4.3f\n", pe);

	if (pe > pe_max) {
	  pe_max = pe;
	}

      }  
      //}	
    }

  }

  //printf("%f\n", pe_max);
  
  return pe_max;
}
