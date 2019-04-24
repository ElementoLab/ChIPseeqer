#include <stdio.h>
#include <stdlib.h>
#include <search.h>
#include <string.h>
#include <ctype.h>

#include "hashtable.h"
#include "sequences.h"
#include "dataio.h"
#include "statistics.h"

typedef struct _ScanACEmatch {
  
  int pos;
  int mlen;

  struct _ScanACEmatch* next;

} ScanACEmatch;

void listGeneMatches(ScanACEmatch* m0);

int segmentHasAnyNs(char* seq, int p, int w);
void loadScanACEMatchesInHash(char* file, int numgenes, struct my_hsearch_data* hash_genes, ScanACEmatch*** ma);
int getMedianDistanceBetweenMotifOccurrences(ScanACEmatch** m1, ScanACEmatch** m2, char** seqs, int numgenes, int rand);
double getMedianDistanceZScore(ScanACEmatch** m1, ScanACEmatch** m2, char** seqs, int numgenes, int num, int med, int rand, int* rank);

int main(int argc, char** argv) 
{
  
  char* fastafile = argv[1];
  char** seqs;
  char** seqn;
  int    numgenes = 0;
  ScanACEmatch** m1;
  ScanACEmatch** m2;

  // subset of genes where motifs coocuur
  ScanACEmatch** m1co;
  ScanACEmatch** m2co;
  
  char** seqsco = 0;

  //char buff[1000];
  //int len = 1000;
  //FILE* fm1, *fm2;
  //char** a; 
  //int    m;
  int verbose = 0;
  struct my_hsearch_data* hash_genes;
  //ScanACEmatch* newm;
  //ENTRY e, *ep;
  // read fasta sequences
  
  default_set_seed(1234);

  loadFastaSequencesAndMakeIndex(fastafile, &seqs, &seqn, &numgenes, &hash_genes);
    
  //fprintf(stderr, "%d seqs\n", numgenes);
  
 
  loadScanACEMatchesInHash(argv[2], numgenes, hash_genes, &m1);
  loadScanACEMatchesInHash(argv[3], numgenes, hash_genes, &m2);

  
  m1co  = (ScanACEmatch**)calloc(numgenes, sizeof(ScanACEmatch*));
  m2co  = (ScanACEmatch**)calloc(numgenes, sizeof(ScanACEmatch*));
  seqsco  = (char**)calloc(numgenes, sizeof(char*));

  int ico = 0;
  int i;
  for (i=0; i<numgenes; i++) {
    if ((m1[i] != 0) && (m2[i] != 0)) {
      m1co[ico] = m1[i];
      m2co[ico] = m2[i];
      //listGeneMatches(m2co[ico]);
      seqsco[ico] = seqs[i];
      ico++;
    }
  }

  // update num genes
  numgenes = ico;

  //printf("Number of genes in which the two motifs cooccur = %d\n", numgenes);
  
  int med1 = getMedianDistanceBetweenMotifOccurrences(m1co,m2co, seqsco, numgenes, 0);
  if (verbose == 1) printf("real median = %d\n", med1);

  // shuffle mode 2 where gene labels are shuffled for second motif
  int rank1 = 0;
  double z1 = getMedianDistanceZScore(m1co, m2co, seqsco, numgenes, 1000, med1, 2, &rank1);
  if (verbose == 1) 
    printf("z-score = %f, rank=%d\n", z1, rank1);

  // shuffle mode 1 where positions are drawn at random
  int rank2 = 0;
  double z2 = getMedianDistanceZScore(m1co, m2co, seqsco, numgenes, 1000, med1, 1, &rank2);
  if (verbose == 1) 
    printf("z-score = %f, rank=%d\n", z2, rank2);

  printf("%s\t%s\t%d\t%d\t%4.3f\t%d\t%4.3f\t%d\n", mybasename(argv[2]), mybasename(argv[3]), numgenes, med1, z1, rank1, z2, rank2);

  return 0;
}


void listGeneMatches(ScanACEmatch* m0)
{
  ScanACEmatch* mptr;

  mptr = m0;
  int p;

  while (mptr != 0) {
    p    = mptr->pos;
    printf("p=%d\n", p);
    mptr = mptr->next;
  }
}

double getMedianDistanceZScore(ScanACEmatch** m1, ScanACEmatch** m2, char** seqs, int numgenes, int num, int med, int rand, int* rank)
{  
  int i;
  //double zscore;
  int* a_d = (int*)calloc(num, sizeof(int));
  *rank = 0;
  for (i=0; i<num; i++) {
    int med2 = getMedianDistanceBetweenMotifOccurrences(m1, m2, seqs, numgenes, rand);
    //printf("rand median = %d\n", med2);
    a_d[i] = med2;
    if (med2 < med)
      (*rank) ++;
  }
  

  double s = stddev_int(a_d, num);
  double m = average_int(a_d, num);
  
  return ( med - m ) / s;

}


int getMedianDistanceBetweenMotifOccurrences(ScanACEmatch** m1, ScanACEmatch** m2, char** seqs, int numgenes, int rand)
{
  int i;
  int dmin, d;
  ScanACEmatch* mptr1;
  ScanACEmatch* mptr2;
  
  int* a_d = (int*)calloc(numgenes, sizeof(int));
  int numd = 0;
  int p1, p2;
  int* randidx = 0;

  if (rand == 2)
    randidx = randPermOfIndices(numgenes); 


  for (i=0; i<numgenes; i++) {

    //if (rand == 1)
    //  printf("num=%d\n", i);
   
    if (rand == 2) {
      if ((m1[i] == 0) || (m2[i] == 0)) {
	die("can't have empty matches in that randomization mode\n");
      }
    }
    
    dmin = strlen(seqs[i]);
    
    mptr1 = m1[i];
    while (mptr1 != 0) {

      p1    = mptr1->pos;
      
      
      if (rand == 1) {				
	int trials = 0;
	while (1) {

	  p1 = (int)(0.5+ default_rand() * (strlen(seqs[i]) - mptr1->mlen ));
	  //printf("p1=%d\n", p1);
	  trials ++;
	  if ((!segmentHasAnyNs(seqs[i], p1, mptr1->mlen)) || (trials == 10))
	    break;

	}
      }
      


      mptr2 = m2[i];
      
      // modify
      if (rand == 2)
	mptr2 = m2[ randidx[i] ];

      
      while (mptr2 != 0) {
	
	p2    = mptr2->pos;
	
	if (rand == 1) {
	  int trials = 0;
	  while (1) {

	    p2 = (int)(0.5+ default_rand() * ( strlen(seqs[i]) - mptr2->mlen) );
	    //printf("p2=%d, drawn fron seq of len = %d  (minus %d, orig pos = %d) \n", p2, (int)strlen(seqs[i]), mptr2->mlen, mptr2->pos);
	    //printf("p2=%d\n", p2);

	    trials ++;
	    if ((!segmentHasAnyNs(seqs[i], p2, mptr2->mlen)) || (trials == 10))
	      break;
	  }	  	  
	}
	
	

	d = abs(p1 - p2);
	if (d < dmin)
	  dmin = d;

	mptr2 = mptr2->next;
      }

      mptr1 = mptr1->next;
    }

    //printf("%s\t%d\n", seqn[i], dmin);
    if (dmin < strlen(seqs[i])) {
      a_d[ numd ] = dmin;
      numd++;
    }

  }
  
  if (rand == 2)
    free(randidx);


  int med = median_int(a_d, numd); 
  free(a_d);

  return med;

}

int segmentHasAnyNs(char* seq, int p, int w)
{
  //int hasN = 0;
  int i;

  for (i=p; i<p+w; i++) {
    if (toupper(seq[i]) == 'N')
      return 1;      
  }
  return 0;
}


void loadScanACEMatchesInHash(char* file, int numgenes, struct my_hsearch_data* hash_genes, ScanACEmatch*** ma)
{
  char buff[10000];
  int len = 10000;
  FILE* fm;
  char** a; 
  int    m;
  int   i;
  ScanACEmatch* newm;
  ENTRY e, *ep;
  
  // read scanace matches for motif 1
  *ma  = (ScanACEmatch**)calloc(numgenes, sizeof(ScanACEmatch*));
  for (i=0; i<numgenes; i++) 
    (*ma)[i] = 0;

  fm = fopen(file, "r");
  if (!fm) {
    die("cannot open fm\n");
  }
  while (fgets(buff, len, fm) != 0) {
    chomp(buff);
    ///printf("buff=%s\n", buff);
    split_line_delim_sep(buff, "\t", &a, &m);
    
    // lookup gene
    e.key = a[0];        
    my_hsearch_r(e, FIND, &ep, hash_genes) ;  
    if (!ep) {
      continue;
    }
 
    int idx = (int)(ep->data);
    
    
    // add
    newm       = (ScanACEmatch*)calloc(1, sizeof(ScanACEmatch));
    newm->pos  = atoi(a[1]);
    newm->mlen = strlen(a[4]);
    if (newm->mlen > 100) {
      die("WTF\n");
    }
    newm->next = (*ma)[ idx ];
    (*ma)[idx]    = newm;
    

    free(a);
  }
  
}
