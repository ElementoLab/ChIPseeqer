#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)
#define _GNU_SOURCE


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <search.h>
#include <math.h>

#include <limits.h>

#include "hashtable.h"
#include "statistics.h"
#include "dataio.h"
#include "uthash.h"
#include "sequences.h"



void readChrData(char* file, char*** chrnames, int** chrlens, int* numchroms)
{
  char* buff = 0;
  int m;
  FILE* fp;
  char** p;
  int numlines = 0;
  int i;
  int  mynmax = 10000;

  buff  = (char*)malloc(mynmax * sizeof(char));
  fp    = fopen(file, "r");
  if (!fp) {
    printf("Please enter a valid filename (%s invalid)\n", file);
    exit(0);
  }

  numlines = nbLinesInFile(file); 
  *chrnames = (char**)calloc(numlines, sizeof(char*));
  *chrlens  = (int*)  calloc(numlines, sizeof(int));
  *numchroms = numlines;
  // read lines
  i = 0;
  while (fgets(buff, mynmax, fp) != 0) {
    chomp(buff);
    split_line_delim(buff, "\t", &p, &m);
    (*chrnames)[i] = strdup(p[0]);
    (*chrlens) [i] = atoi(p[1]);
    free(p);
    i++;
  }
  free(buff);
}



void loadFastaSequencesAndMakeIndex(char* fastafile, char*** seqs, char*** seqn, int* numgenes, struct my_hsearch_data** hash_genes)
{

  seqI  si;
  char* seq;
  char* name;
  int   size;
  int   i, m;
  ENTRY e, *ep;
  int   hashret;

  // read sequence data
  int numseqs = numSequencesInFastaFile(fastafile); 
    
  *hash_genes = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(numseqs*2, *hash_genes);
  if (hashret == 0) {
    printf("Could not create hash table ...\n");
    exit(0);
  }
 
  *seqs = (char**)calloc( numseqs, sizeof(char*));
  if (*seqs == 0) {
    die("could not alloc *seqs\n");
  }
  *seqn = (char**)calloc( numseqs,  sizeof(char*));
  if (*seqn == 0) {
    die("could not alloc *seqn\n");
  }

  si.verbose = 0;
  seqI_set_max_seqlen(&si, 5000000);
  seqI_set_seqlen_inc(&si, 3000000);
  seqI_set_fread_chunk(&si, 1000000);
  
  if (seqI_open(&si, fastafile) == 0) {
    die("Error opening file\n");
  }
  
  int idxgene = 0;
  while ((seq = seqI_nextSequence(&si, &name, &size))) {
    
    m = strlen(name);
    for (i=0; i<m; i++) {
      if (name[i] == ' ') {
	name[i] = '\0';
	break;
      }
    }
    
    // lookup sequence name
    e.key  = strdup(name);   
    e.data = (char*)idxgene;

    my_hsearch_r(e, ENTER, &ep, *hash_genes) ;
    if (!ep) {
      die("Couldn't enter data\n");
    } 

    (*seqs)[idxgene] = seq;
    (*seqn)[idxgene] = name;  
    //seql[idxgene] = strlen(seq);

    //fprintf(stderr, "Loaded %s               \r", name);
    
    idxgene++;
  }
  seqI_close(&si);

  *numgenes = idxgene;
}


double aaConsIndex(char* str) {
  
  char* aas = "ARNDCEQGHILKMFPSTWYV";
  int   n   = 20;
 
  int   l = strlen(str);
  int   i, j;
  int  aacnt[20];
  double sum = 0;
  int   m = 0;
  double   px = 0.0;
  for (i=0; i<n; i++)
    aacnt[i] = 0;

  for (j=0; j<l; j++) {
    if ((str[j] != '-') && (str[j] != 'X')) {
      for (i=0; i<n; i++) {
	if (str[j] == aas[i]) {
	  aacnt[i] ++;
	  m++;
	  break;
	}
      }
    }
  }
  
  for (i=0; i<n; i++) {
    px = aacnt[i]/(double)m;
    if (px > 0) {
      px *= log(px)/log(2.0);
      sum += px;
    }    
  }
    
  return -sum;

}


int chr2index(char* chr) 
{
  int l = strlen(chr);
  //die("Using chr2index\n");
  if (l > 5) 
    return -1;  
  else if (chr[3] == 'M')
    return 0;
  else if (chr[3] == 'X')
    return 23;
  else if (chr[3] == 'Y')
    return 24;
  else {
    return atoi(chr+3);    
  }
  
}

char* IUPAC2nt(char c) {
  int    i;
  char*  d_ch = "NACGTMRWSYKVHBD";
  char*  m_ch[] = {"ACGT",
		   "A",
		   "C",
		   "G",
		   "T",
		   "AC",
		   "AG",
		   "AT",
		   "CG",
		   "CT",
		   "GT",
		   "ACG",
		   "ACT",
		   "CGT",
		   "AGT" };
  for (i=0; i<15; i++) {
    if (d_ch[i] == c) {
      return m_ch[i];
    }
  }
  return 0;
}



int numSequencesInFastaFile(char* file) 
{
  
  char* buff;
  FILE* fp;
  int   n;
  int   maxn = 1000000;

  buff = (char*)malloc(maxn * sizeof(char));
  fp = fopen(file, "r");
  if (!fp) {
    printf("could not open seq data %s\n", file);
  } 
  n = 0;
  fgets(buff, maxn, fp);
  if (buff[0] == '>')
    n++;
  while (fgets(buff, maxn, fp) != 0) {
    if (buff[0] == '>')
      n++;

  }  
  fclose(fp);

  free(buff);
  return n;
}



char* strrev(char* str) 
{
  int i, j;
  int l = strlen(str);
  char* rev = (char*)calloc(l+1, 1);
  
  for (i=l-1,j=0; i>=0; i--,j++)
    rev[j] = str[i];

  return rev;
}

void readSNPfile(char* file, char** a_snp, char*** snp_names)
{
  FILE* f;
  // reading file
  char*   buff;
  int     mynmax          = 100000;
  char**  p;
  int     m;
  char*   newc;

  buff   = (char*)calloc(mynmax, sizeof(char)); 


  f = fopen(file, "r");
  while (!feof(f)) {

    fgets(buff, mynmax, f);
    if (feof(f))
	break; 
    chomp(buff);
    split_line_delim(buff, "\t", &p, &m);

    (*a_snp)    [ atoi(p[1]) ] = (char)1;
    if ((*snp_names)[ atoi(p[1]) ] == 0) 
      (*snp_names)[ atoi(p[1]) ] = strdup(p[0]);
    else {
      int l = strlen(p[0]);
      newc  = (char*)calloc(strlen(  (*snp_names)[ atoi(p[1]) ] ) + 1 + l + 1, sizeof(char));
      strcat(newc,  (*snp_names)[ atoi(p[1]) ]);
      strcat(newc, "/");
      strcat(newc, p[0]);			
      free( (*snp_names)[ atoi(p[1]) ]);
      (*snp_names)[ atoi(p[1]) ] = newc;
    }
    
    free(p); // p[x] are just pointers

  }
  
  //fprintf(stderr, "read dbSNP file %s\n", file);
}




void seq2lower(char* s) {
  
  int i;
  int l = strlen(s);
  for (i=0; i<l; i++) 
    s[i] = tolower(s[i]);
}

void seqlower2dots(char* s) 
{
  
  int i;
  int l = strlen(s);
  for (i=0; i<l; i++) {
    if (s[i] == tolower(s[i])) 
      s[i] = '.';
  }
}

/**
 *
 * input:
 * 	wm
 * 	w:	width of weight matrx
 * 	m:	number of sequences to generate
 *
 * output:
 * 	seqs:	random sequences sampled from wm
 *
 */
void generateManyRandomSequences(float** wm, int w, int m, char*** seqs)
{
  int i=0;

  (*seqs) = (char**)malloc(m * sizeof(char*));

  for (i=0; i<m; i++) {
    generateRandomSequence(wm, w, &((*seqs)[i]));
    //printf("%s\n", (*seqs)[i]);
  }

}

//
// sample a sequence from the motif distribution
//
// output: 
// 	s	sequence
//
void generateRandomSequence(float** wm, int w, char** s) 
{

  int i, j = 0, k;
  char* nt = "ACTG";
  float a_cum[5];
  float d = 0;
  char c = '-';

  *s = (char*)calloc(w+1, sizeof(char));
  
  for (j=0; j<w; j++) {
	
    i = 1;
    a_cum[0] = 0.0;
    for (k=0; k<4; k++) {
      a_cum[i] = a_cum[i-1] + wm[j][k];
      
      i++;
    }
    a_cum[i-1] = 1.0;

    d = default_rand();
    /*
    printf("\n");
    for (i=0; i<5; i++)
      printf("%f\n", a_cum[i]);
    printf("d=%f\n", d);
    */

    for (i=1; i<=4; i++) {
      if ((d > a_cum[i-1]) && (d <= a_cum[i])) {
	c = nt[i-1];
	break;
      }
    }
      
    
    (*s)[j] = c;
  }
    
  (*s)[w] = '\0';
}


int basicEditDistance(char* ss1, char* ss2) 
{
  int w1 = strlen(ss1);
  int w2 = strlen(ss2);
  int mmf = 0, mmr = 0, i=0, j=0;
  if (w1 != w2) {
    die("Trying to calc edit distance based on s of unequal length\n");
  }
  
  for (i=0, j=w1-1; i<w1; i++,j--) {
    if (ss1[i] != ss2[i])
      mmf ++;
    if (ss1[i] != char_complement(ss2[j]))
      mmr ++;
  }
  
  return min(mmr, mmf);
}

void sw(char* ss1, char* ss2, int R, int P, int O, int G, float* id, int* alnlen, int* qst, int* qen, int* sst, int* sen, int* sco, char** aln1, char** aln2) 
{

  int**  A;  // scores
  int*** L; // links 
  char*  s1;
  char*  s2;
  int    newsc;
  int    n1 = strlen(ss1);
  int    n2 = strlen(ss2);
  int    j_idx, i_idx;
  int    i, j, k;
  int    myi, myj;

  int    total_max_j = 0;
  int    total_max_i = 0;
  int    total_maxsc = 0;
  int    verbose = 0;
  int    maxsc;
  int    sc;
  int    maxil;
  int    maxjl;
  char*  as1;
  char*  as2;
  int    ii, jj, li, lj;

  s1 = (char*)calloc(n1+2, sizeof(char));
  s2 = (char*)calloc(n2+2, sizeof(char));

  s1[0] = 'N';
  memcpy(s1+1, ss1, n1);
  s1[n1+1] = '\0';

  s2[0] = 'N';
  memcpy(s2+1, ss2, n2);
  s2[n2+1] = '\0';

  n1++;
  n2++;

  // create score array
  A = (int**)calloc(n1, sizeof(int*));  
  for (i=0; i<n1; i++) {    
    A[i] = (int*)calloc(n2, sizeof(int));
  }

  for (i=0; i<n1; i++) {
    A[i][0] = 0;
  }      

  for (j=0; j<n2; j++) {
    A[0][j] = 0;
  }      

  // create link array
  L = (int***)calloc(n1, sizeof(int**));  
  for (i=0; i<n1; i++) {    
    L[i] = (int**)calloc(n2, sizeof(int*));
    for (j=0; j<n2; j++)
      L[i][j] = (int*)calloc(2, sizeof(int));
  }
  
  L[0][0][0] = -1;
  L[0][0][0] = -1;

  for (i=0; i<n1; i++) {    
    for (j=0; j<n2; j++) {

      if ((i == 0) && (j == 0))
	continue;
      
      if (s1[i] == s2[j])
	sc = R;
      else 
	sc = G;
      
      maxsc =  0;
      maxil = -1;
      maxjl = -1;

      // diag score
      if ((i>0) && (j>0)) {
	newsc = A[i-1][j-1] + sc;
	if (newsc > maxsc) {
	  maxsc = newsc;
	  maxil = i-1;
	  maxjl = j-1;
	}
      } // if
      
      // up
      for (ii = i-1; ii>=0; ii--) {
	newsc = A[ii][j] + (i-ii) * G + O;
	if (newsc > maxsc) {
	  maxsc = newsc;
	  maxil = ii;
	  maxjl = j;
	}
      }
      
      // left
      for (jj = j-1; jj>=0; jj--) {
	newsc = A[i][jj] + (j-jj) * G + O;
	if (newsc > maxsc) {
	  maxsc = newsc;
	  maxil = i;
	  maxjl = jj;
	}
      }
      
      // store links
      L[i][j][0] = maxil;
      L[i][j][1] = maxjl;
  
      // store max score
      A[i][j]    = maxsc;
      
      // keep track of max score
      if (maxsc > total_maxsc) {
	total_maxsc = maxsc;
	total_max_i = i;
	total_max_j = j;
      }

    }
  }
  

  
  if (verbose == 1) {
    printf("    ");
    for (j=0; j<n2; j++) {
      printf("%4d", j);
    }
    printf("\n");
  
    for (i=0; i<n1; i++) { 
      printf("%4d", i);
      for (j=0; j<n2; j++) {
	printf("%4d", A[i][j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  // done with filling array
  i = total_max_i;
  j = total_max_j;
  
  *qen =  i-1;  // because added N to sequence ...
  *sen =  j-1;
  *qst = -1;
  *sst = -1;
  
  i_idx = 0;
  j_idx = 0;

  // string
  int maxal = max(n1,n2)*2;
  as1 = (char*)calloc(maxal, sizeof(char));
  as2 = (char*)calloc(maxal, sizeof(char));

  while (A[i][j] > 0) {

    *qst = i;
    *sst = j;
        
    myi = L[i][j][0];
    myj = L[i][j][1];
    
    //printf("Go to myi=%d, myj=%d\n", myi, myj);
    
    //if ((myi == -1) && (myj == -1))
    //  break;

    li = i - myi;
    lj = j - myj;


    if (li == 0) {
      for (k=j; k>myj; k--) {
	as1[i_idx++] = '-';
	as2[j_idx++] = s2[k];
      }
    } else if (lj == 0) {
       for (k=i; k>myi; k--) {
	 as1[i_idx++] = s1[k];
	 as2[j_idx++] = '-';
       }
    } else {
      as1[i_idx++] = s1[i];
      as2[j_idx++] = s2[j];
    }
    
    i = myi;
    j = myj;


  }

  *qst = *qst - 1;
  *sst = *sst - 1;

  // reverse strings  
  char* rev_as1 = (char*)calloc(i_idx+1, 1);
  char* rev_as2 = (char*)calloc(j_idx+1, 1);
  
  for (i=i_idx-1,j=0; i>=0; i--,j++)
    rev_as1[j] = as1[i];

  for (i=j_idx-1,j=0; i>=0; i--,j++)
    rev_as2[j] = as2[i];

  rev_as1[i_idx] = '\0';
  rev_as2[j_idx] = '\0';
  
  //printf("rev_sa1=%s\n", rev_as1);
  //printf("rev_sa2=%s\n", rev_as2);

  free(s1);
  free(s2);
  free(as1);
  free(as2);

  
  for (i=0; i<n1; i++)
    free(A[i]);
  free(A);
  
  
  for (i=0; i<n1; i++) {    
    for (j=0; j<n2; j++)
      free(L[i][j]);
    free(L[i]);
  }
  
  free(L);
  
  // sim
  for (i=0; i<i_idx; i++) {    
    if (rev_as1[i] == rev_as2[i]) 
      (*id) += 1;
  }
  *id     = *id / i_idx;

  *aln1   = rev_as1;
  *aln2   = rev_as2;
  *alnlen = i_idx;
  *sco    = total_maxsc;
  

}


int rand_location_around_pos_in_chrlen(int pos, int range, int chrlen)
{
  int randpos = -10;
  while ((randpos < 0) || (randpos >= chrlen) || (abs(pos-randpos)<5)) {
    randpos = pos + (int)(0.5 + default_rand() * range - range / 2);
  }
  return randpos;
}

int rand_location_around_pos(int pos, int range, char* c)
{
  int randpos = -10;
  
  int chrl = hg18_chrlen(c);
  if (chrl == -1) {
    printf("Not recognizing %s\n", c);
    exit(0);
  }

  while ((randpos < 0) || (randpos >= chrl) || (abs(pos-randpos)<5)) {
    randpos = pos + (int)(0.5 + default_rand() * range - range / 2);
  }

  //printf("oldpos = %d, newpos = %d\n", pos, randpos);

  return randpos;
}

int rand_chr_location(char* c)
{
  int chrl = hg18_chrlen(c);
  if (chrl == -1) {
    printf("Not recognizing %s\n", c);
    exit(0);
  }

  int i    = (int)(0.5 + default_rand() * chrl);

  return i;
}


int rand_chr_interval_chrlen(int chrl, int l)
{

  int i    = (int)(0.5 + default_rand() * chrl - l);

  return i;
}


int rand_chr_interval(char* c, int l)
{
  int chrl = hg18_chrlen(c);
  //printf("chrl=%d\n", chrl);
  if (chrl == -1) {
    printf("Not recognizing %s\n", c);
    exit(0);
  }

  int i    = (int)(0.5 + default_rand() * chrl - l);

  return i;
}


float hg18_30mer_mappability(char* c) 
{
  int i;
  char* chrnames[] = 
    {"chr18",
     "chr2",
     "chr3",
     "chr4",
     "chr6",
     "chr8",
     "chr5",
     "chr12",
     "chr11",
     "chr20",
     "chr7",
     "chr10",
     "chr17",
     "chrX",
     "chr1",
     "chr13",
     "chr16",
     "chr14",
     "chrM",
     "chr9",
     "chr19",
     "chr15",
     "chr21",
     "chr22",
     "chrY"};
  
  float chrmap[] = {
    0.123087459668913,
    0.139325292921335,
    0.139622164963933,
    0.141630052737745,
    0.142365618132972,
    0.143106107677065,
    0.143924633059643,
    0.14645199279659,
    0.149044906485258,
    0.160449000194824,
    0.164855517225434,
    0.169669404417753,
    0.184021980040252,
    0.201849540099583,
    0.223202247602959,
    0.251159230291692,
    0.265147676410215,
    0.279122120502026,
    0.279283084907368,
    0.307441635415995,
    0.317359834491667,
    0.317702957023205,
    0.359312392256674,
    0.426161556382597,
    0.747921665906161,
  };
  
  for (i=0; i<25; i++) {
    if (strcmp(chrnames[i], c) == 0) {
      return 1.0 - chrmap[i];
    }
  }	
  return -1;
}


void get_hg18_chroms(char*** chroms, int* numchroms, int** chrlens)
{
  int i;
  *numchroms   = 24;
  *chroms      = (char**)calloc(*numchroms, sizeof(char*));
  (*chroms)[0] =  strdup("chrY");		//array that holds the chromosomes names
  (*chroms)[1] =  strdup("chrX");
  (*chroms)[2] =  strdup("chr9");
  (*chroms)[3] =  strdup("chr8");
  (*chroms)[4] =  strdup("chr7");
  (*chroms)[5] =  strdup("chr6");
  (*chroms)[6] =  strdup("chr5");
  (*chroms)[7] =  strdup("chr4");
  (*chroms)[8] =  strdup("chr3");
  (*chroms)[9] =  strdup("chr22");
  (*chroms)[10] = strdup("chr21");
  (*chroms)[11] = strdup("chr20");
  (*chroms)[12] = strdup("chr2");
  (*chroms)[13] = strdup("chr19");
  (*chroms)[14] = strdup("chr18");
  (*chroms)[15] = strdup("chr17");
  (*chroms)[16] = strdup("chr16");
  (*chroms)[17] = strdup("chr15");
  (*chroms)[18] = strdup("chr14");
  (*chroms)[19] = strdup("chr13");
  (*chroms)[20] = strdup("chr12");
  (*chroms)[21] = strdup("chr11");
  (*chroms)[22] = strdup("chr10");
  (*chroms)[23] = strdup("chr1");

  *chrlens      = (int*)calloc(*numchroms, sizeof(int));
  for (i=0; i<*numchroms; i++) {
    (*chrlens)[i] = hg18_chrlen( (*chroms)[i] );
  }

}

int hg18_chrlen(char* c) 
{
  
  int i;
  char* chrnames[] = 
    {"chr1",
     "chr2",
     "chr3",
     "chr4",
     "chr5",
     "chr6",
     "chr7",
     "chrX",
     "chr8",
     "chr9",
     "chr10",
     "chr11",
     "chr12",
     "chr13",
     "chr14",
     "chr15",
     "chr16",
     "chr17",
     "chr18",
     "chr19",
     "chr20",
     "chrY",
     "chr22",
     "chr21",
     "chrM"};

  int chrlens[] = {247249719,
		   242951149,
		   199501827,
		   191273063,
		   180857866,
		   170899992,
		   158821424,
		   154913754,
		   146274826,
		   140273252,
		   135374737,
		   134452384,
		   132349534,
		   114142980,
		   106368585,
		   100338915,
		   88827254,
		   78774742,
		   76117153,
		   63811651,
		   62435964,
		   57772954,
		   49691432,
		   46944323,
		   16571};
  
  for (i=0; i<25; i++) {
    if (strcmp(chrnames[i], c) == 0) {
      return chrlens[i];
    }
  }	
  return -1;
}



long sequencesOverlap(long s1, long e1, long s2, long e2) {

  long p1 = max(s1, s2);    
  long p2 = min(e1, e2);
  


  long ll = p2 - p1 + 1;
  
  //printf("%ld\t%ld\tll=%ld\n", p1,p2,ll);

  if ( ll > 0  ) {
    return ll;
  } else {
    return 0;
  }
    
}


//  P(A1,A2,A3) = P(A1)P(A2)P(A3)  (all col indep)
//  P(A1,A2,A3) = P(A1,A2)P(A3|A1,A2)   = P(A1)P(A2,A3|A1) = P(A1)
//                        P(A3|A1,A2) = 
//              = P(A1)P(A2|A1)P(A3|A2)

//  comes down to P(A3|A1,A2) = P(A3|A2) if A1 and A3 are indep, ie P(A1,A3)=0 
 

// input a file, returns an array
void readScanACEMatches(char* file, ScanACE_match** scm, int* nmp) 
{

  char* buff;
  int   mynmax = 1000;
  FILE* fp;
  char* genename;
  int   position;
  int   strand;  
  float score;   
  int   maxscanacematches = 100000;
  int   nm = 0;
  //
  // read scanace matches
  //  
  *scm = (ScanACE_match*)malloc(maxscanacematches * sizeof(ScanACE_match));
    
  buff  = (char*)malloc(mynmax * sizeof(char));
  fp    = fopen(file, "r");
  if (!fp) {
    printf("Please enter a valid filename (%s invalid)\n", file);
    exit(0);
  }
   
  
  // read matches
  while (!feof(fp)) {
    
    fgets(buff, mynmax, fp);
    
    if (feof(fp))
      break; 
    
    chomp(buff);
    
    //
    // get the row name 
    //
    genename =    mystrtok(buff, '\t');
    position = atoi(mystrtok(0, '\t'));
    strand   = atoi(mystrtok(0, '\t'));
    score    = atof(mystrtok(0, '\t'));
	
    (*scm)[nm].gene     = genename;
    (*scm)[nm].strand   = strand;
    (*scm)[nm].position = position;
    (*scm)[nm].score    = score;
	
    nm ++;
	
    if (nm == maxscanacematches) {
      die("maxscanacematches reached, please modify source and recompile.\n");
    }
    
  }
  

  *nmp = nm;
  
  fclose(fp);

  
}

/*
 * 
 * TO DO:
 * 	dangerous calling hidden global variables!
 *
 */
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
  
  C['A']  = 2;
  C['C']  = 3;
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


void getWMConsensus(float** wm, int w, char** consensus)
{
  int i,j;
  float max_pi;
  //int max_pi_idx;
  char *nt = "ACTG";

  *consensus = (char*)calloc((w+1), sizeof(char));

  for (i=0; i<w; i++) {
    // nt at pos i
    max_pi = -1000.0;
    for (j=0; j<4; j++) {
      //printf("max_pi = %f\n", max_pi);
      if (wm[i][j] > max_pi) {
	max_pi = wm[i][j];
	(*consensus)[i] = nt[j];
      }
    }
  }
  
  (*consensus)[w] = '\0';

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

void printBckd(float* m)
{
  int i;
  printf("BKG");
  for (i=0; i<4; i++) {
    printf("\t%4.3f", m[i]);
  }
  printf("\n");
  
}

//
// read PBM matrix 
//
void read_BULYK_WM(char* wmfile, float*** wm, int* w)
{
  char* buff;
  int   len = 100000;
  FILE* fp = 0;
  char** p;
  int    m, i, col;

  buff = (char*)malloc(1000*sizeof(char));

  fp = fopen(wmfile, "r");

  // order in matrices is ACGT
  // order in outout is ACTG

  // rid of row 1
  fgets(buff, len, fp);
  
  //
  // Motif line 1
  //
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, "\t", &p, &m);

  // must count nt
  *w = m-1;

  //
  // allocate memory for matrix (attention matrix is inverted)
  //
  *wm = (float**)malloc((m-1) * sizeof(float*));
  for (col= 0; col < m-1; col++) {
    (*wm)[col] = (float*)malloc(4 * sizeof(float));
    for (i=0; i<4; i++) {
      (*wm)[col][i] = 0;
    }
  }


  // enter data from line 1 (A)
  for (i=1; i<m; i++) {
    //printf("'%s' ", p[i]);
    (*wm)[i-1][0] = atof(p[i]);
  }
  //printf("\n");
  free(p);


  // read data from line 2 (C)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, "\t", &p, &m); 
  for (i=1; i<m; i++) {
    (*wm)[i-1][1] = atof(p[i]);
  }
  free(p);

  // read data from line 3 (T, is line 3 in 0..3 scale)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, "\t", &p, &m);   
  for (i=1; i<m; i++) {
    (*wm)[i-1][3] = atof(p[i]);
  }
  free(p);

  // read data from line 3 (G, is line 2 in 0..3 scale)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, "\t", &p, &m);   
  for (i=1; i<m; i++) {
    (*wm)[i-1][2] = atof(p[i]);
  }
  free(p);

  
}



void read_JASPAR_WM(char* wmfile, int*** wm, int* w, int** nsites_j)
{
  char* buff;
  int   len = 100000;
  FILE* fp = 0;
  char** p;
  int    m, i, j, col;

  buff = (char*)malloc(1000*sizeof(char));

  fp = fopen(wmfile, "r");

  // order in matrices is ACGT
  // order in outout is ACTG

  //
  // Motif line 1
  //
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, " ", &p, &m);

  //for (col=0; col<m; col++) {
  //  printf("%s\n", p[col]);
  //}

  *w = m;

  // note : first spaces are not counted in split_line_delim
  //
  // allocate memory for matrix (attention matrix is inverted)
  //
  *wm = (int**)malloc(m * sizeof(int*));
  for (col= 0; col < m; col++) {
    (*wm)[col] = (int*)malloc(4 * sizeof(int));
    for (i=0; i<4; i++) {
      (*wm)[col][i] = 0;
    }
  }

  *nsites_j = (int*)calloc(m, sizeof(int));

  // enter data from line 1 (A)
  for (i=0; i<m; i++) {
    //printf("'%s' ", p[i]);
    (*wm)[i][0] = atoi(p[i]);
  }
  //printf("\n");
  free(p);


  // read data from line 2 (C)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, " ", &p, &m); 
  
  for (i=0; i<m; i++) {
    //printf("%s ", p[i]);
    (*wm)[i][1] = atoi(p[i]);
  }
  //printf("\n");
  free(p);

  // read data from line 3 (G, is line 3 in 0-3)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, " ", &p, &m); 
  
  for (i=0; i<m; i++) {
    //printf("%s ", p[i]);
    (*wm)[i][3] = atoi(p[i]);
  }
  //printf("\n");
  free(p);

  // read data from line 3 (G, is line 3 in 0-3)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, " ", &p, &m); 
  
  for (i=0; i<m; i++) {
    //printf("%s ", p[i]);
    (*wm)[i][2] = atoi(p[i]);
  }
  //printf("\n");
  free(p);

  
  // calculate nsites
  //int oldsum = 0;
  //int maxsum = -1;
  for (i=0; i<m; i++) {
    int sum = 0;
    for (j=0; j<4; j++) {
      sum += (*wm)[i][j];
    }
    (*nsites_j)[i] = sum;
    //printf("sum=%d\n", sum);
    //if ((i > 0) && (sum != oldsum)) {
    //fprintf(stderr, "Warning: Problem at col %d: obs=%d != %d\n", i, sum, oldsum);
    //exit(0);
    //}
    //oldsum = sum;     
    //if (sum > maxsum)
    //  maxsum = sum;
  }

  //*nsites = maxsum;
  
}




void read_JASPAR_WM_OLD(char* wmfile, int*** wm, int* w, int* nsites)
{
  char* buff;
  int   len = 100000;
  FILE* fp = 0;
  char** p;
  int    m, i, j, col;

  buff = (char*)malloc(1000*sizeof(char));

  fp = fopen(wmfile, "r");

  // order in matrices is ACGT
  // order in outout is ACTG

  //
  // Motif line 1
  //
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, " ", &p, &m);

  //for (col=0; col<m; col++) {
  //  printf("%s\n", p[col]);
  //}

  *w = m;

  // note : first spaces are not counted in split_line_delim
  //
  // allocate memory for matrix (attention matrix is inverted)
  //
  *wm = (int**)malloc(m * sizeof(int*));
  for (col= 0; col < m; col++) {
    (*wm)[col] = (int*)malloc(4 * sizeof(int));
    for (i=0; i<4; i++) {
      (*wm)[col][i] = 0;
    }
  }


  // enter data from line 1 (A)
  for (i=0; i<m; i++) {
    //printf("'%s' ", p[i]);
    (*wm)[i][0] = atoi(p[i]);
  }
  //printf("\n");
  free(p);


  // read data from line 2 (C)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, " ", &p, &m); 
  
  for (i=0; i<m; i++) {
    //printf("%s ", p[i]);
    (*wm)[i][1] = atoi(p[i]);
  }
  //printf("\n");
  free(p);

  // read data from line 3 (G, is line 3 in 0-3)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, " ", &p, &m); 
  
  for (i=0; i<m; i++) {
    //printf("%s ", p[i]);
    (*wm)[i][3] = atoi(p[i]);
  }
  //printf("\n");
  free(p);

  // read data from line 3 (G, is line 3 in 0-3)
  fgets(buff, len, fp);
  chomp(buff);
  split_line_delim(buff, " ", &p, &m); 
  
  for (i=0; i<m; i++) {
    //printf("%s ", p[i]);
    (*wm)[i][2] = atoi(p[i]);
  }
  //printf("\n");
  free(p);

  
  // calculate nsites
  int oldsum = 0;
  int maxsum = -1;
  for (i=0; i<m; i++) {
    int sum = 0;
    for (j=0; j<4; j++) {
      sum += (*wm)[i][j];
    }
    //(*nsites_j)[i] = sum;
    //printf("sum=%d\n", sum);
    if ((i > 0) && (sum != oldsum)) {
      fprintf(stderr, "Warning: Problem at col %d: obs=%d != %d\n", i, sum, oldsum);
      //exit(0);
    }
    oldsum = sum;     
    if (sum > maxsum)
      maxsum = sum;
  }

  *nsites = maxsum;
  
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
      mywm[col][i] = 0;
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
// convert int matric into background-corrected weight matrix
//
void intWMtofloatWM(int** intwm, int w, int nsites, float*** fwm) 
{

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
      mywm[col][i]  = intwm[col][i] / (float)nsites;
    }
  }
  
  *fwm = mywm;
  
}


//
// convert int matric into background-corrected weight matrix
//
void intWMtofloatWM_J(int** intwm, int w, int* nsites, float*** fwm) 
{

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
      mywm[col][i]  = intwm[col][i] / (float)(nsites[col]);
    }
  }
  
  *fwm = mywm;
  
}



//
// convert int matric into background-corrected weight matrix
//
// input:
// 	intwm
// 	w
// 	nsites
// 	rawbkg:		default= .52/2, (1-.52)/2, ...
//
// output:
// 	logwm
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
// convert int matric into background-corrected weight matrix (JASPAR) 
//
// input:
// 	intwm
// 	w
// 	nsites
// 	rawbkg:		default= .52/2, (1-.52)/2, ...
//
// output:
// 	logwm
//
void ACEintWMtologWM_J(int** intwm, int w, int* nsites, float* rawbkg, float*** logwm) {

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
      mywm[col][i]  = log2( ((intwm[col][i] + rawbkg[i]) / (float)(nsites[col] + 1)) / rawbkg[i] ); 
    }
  }
  
  *logwm = mywm;
  
}




//
// 
//
void WMtologWMbkg(float** fwm, int w, float* rawbkg, float*** logwm) 
{

  int col, i;
  float** mywm;
  float eps = 0.000001;

  mywm = (float**)malloc(w * sizeof(float*));
  for (col= 0; col<w; col++) {
    mywm[col] = (float*)malloc(5 * sizeof(float));
    for (i=0; i<5; i++) {
      mywm[col][i] = 0.0;
    }
  }

  for (col= 0; col < w; col++) {
    for (i=0; i<4; i++) {	
      mywm[col][i]  = log2( 
			   (fwm[col][i] + eps) / rawbkg[i] 
			   ); 
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

/*
 * input: 
 * 	wml:		weight matrix?
 * 	w:		width?
 *	stars:		??
 * 	bkgl:		default= 0 (if weight matrix has been background corrected)
 * 	t:		average score
 * 	seq:		wild type sequence
 * 	rna:		is RNA
 * 	matches_pos
 * 	num_matches
 * 	matches_ori:	orienation
 * 	matches_sco:
 * 	max_num_matches:
 *
 * example usage:
 *
 * findAllWeightMatrixMatches(jm[j].wm, jm[j].w, jm[j].stars, 0, jm[j].t0, wtseq, 
									0, &wtseq_matches_pos, &wtseq_np, &wtseq_matches_ori, 
									&wtseq_matches_sco, 10000);

 */
void findAllWeightMatrixMatches(float** wml, int w, int* stars, float* bkgl, float t, char* seq, int rna, int** matches_pos, int* num_matches, char** matches_ori, float** matches_sco, int max_num_matches)
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
    printf("findAllWeightMatrixMatches: not enough memory for matches_pos (%ld)\n", my_max_num_matches * sizeof(int));    
    exit(1);
  }
  *matches_ori = (char*)malloc(my_max_num_matches * sizeof(char));
  if (!*matches_ori) {
    die("findAllWeightMatrixMatches: not enough memory for matches_ori\n");    
  }
  *matches_sco = (float*)malloc(my_max_num_matches * sizeof(float));
  if (!*matches_sco) {
    die("findAllWeightMatrixMatches: not enough memory for matches_sco\n");    
  }

  for (i=0; i<l; i++) {

    s = getScoreOnSeq(wml, w, stars, bkgl, seq, i, rna, &strand, &hasn);
    
    if (!hasn && (s >= t)) {

      (*matches_pos)[ my_num_matches ] = i;
      (*matches_ori)[ my_num_matches ] = (char)strand;
      (*matches_sco)[ my_num_matches ] = s;
      
      my_num_matches ++;

      if (my_num_matches == my_max_num_matches) {
	my_max_num_matches += max_num_matches;
	ptr = realloc( *matches_pos, my_max_num_matches * sizeof(int));	
	if (ptr == 0)
	  die("findAllWeightMatrixMatches: not enough memory for matches_pos (realloc)\n");
	ptr = realloc( *matches_ori, my_max_num_matches * sizeof(char));		
	if (ptr == 0)
	  die("findAllWeightMatrixMatches: not enough memory for matches_ori (realloc)\n");
	ptr = realloc( *matches_sco, my_max_num_matches * sizeof(float));		
	if (ptr == 0)
	  die("findAllWeightMatrixMatches: not enough memory for matches_sco (realloc)\n");
	
      }
      
    }    
  } 

  *num_matches = my_num_matches;
  
}

//
// find all motif matches in intervals
//
void findAllWeightMatrixMatchesInIntervals(float** wml, int w, int* stars, float* bkgl, float t, char* seq, int rna, int** matches_pos, int* num_matches, char** matches_ori, float** matches_sco, int max_num_matches, int** intervals, int numint)
{

  int i;
  int l = strlen(seq);
  float s;
  int strand;
  int my_num_matches = 0;
  int my_max_num_matches = max_num_matches;
  void* ptr;
  int  hasn;
  int  j;

  *matches_pos = (int*) malloc(my_max_num_matches * sizeof(int));  
  if (!*matches_pos) {
    printf("findAllWeightMatrixMatches: not enough memory for matches_pos (%ld)\n", my_max_num_matches * sizeof(int));    
    exit(1);
  }
  *matches_ori = (char*)malloc(my_max_num_matches * sizeof(char));
  if (!*matches_ori) {
    die("findAllWeightMatrixMatches: not enough memory for matches_ori\n");    
  }
  *matches_sco = (float*)malloc(my_max_num_matches * sizeof(float));
  if (!*matches_sco) {
    die("findAllWeightMatrixMatches: not enough memory for matches_sco\n");    
  }

  // go thru all intervals
  for (j=0; j<numint; j++) {
    
    for (i=intervals[j][0]; i<=min(l-1,intervals[j][1]); i++) {
      
      s = getScoreOnSeq(wml, w, stars, bkgl, seq, i, rna, &strand, &hasn);
      
      if (!hasn && (s >= t)) {
	
	(*matches_pos)[ my_num_matches ] = i;
	(*matches_ori)[ my_num_matches ] = (char)strand;
	(*matches_sco)[ my_num_matches ] = s;
	
	my_num_matches ++;
	
	if (my_num_matches == my_max_num_matches) {
	  my_max_num_matches += max_num_matches;
	  ptr = realloc( *matches_pos, my_max_num_matches * sizeof(int));	
	  if (ptr == 0)
	    die("findAllWeightMatrixMatches: not enough memory for matches_pos (realloc)\n");
	  ptr = realloc( *matches_ori, my_max_num_matches * sizeof(char));		
	  if (ptr == 0)
	    die("findAllWeightMatrixMatches: not enough memory for matches_ori (realloc)\n");
	  ptr = realloc( *matches_sco, my_max_num_matches * sizeof(float));		
	  if (ptr == 0)
	    die("findAllWeightMatrixMatches: not enough memory for matches_sco (realloc)\n");
	  
	} // if
	
      } // if     
    } // for (within int)
  } // for each int

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

  int*  my_a_idx = 0;
  
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
  float s = -100000;
  int   hasn, strand;

  l = strlen(seq);
  *idx = -1;

  //printWM(logwm, w);
      
  for (i=0; i<l-w; i++) {    

    if (mt == 0) 
      s = getScoreOnSeq(logwm, w, stars, bkg, seq, i, rna, &strand, &hasn);	    

    else if (mt == 1)
      s = getScoreOnSeq_1M(logwm, w, stars, bkg, seq, i, rna, &strand, &hasn);	 

    //printf("new max for %s at i=%d, score = %f\n", substr(seq, i, w), i, s);

    if ((hasn == 0) && (s > maxs)) {      
      
      //printf("new max for %s, score = %f\n", substr(seq, i, w), s);
      
      maxs = s;
      *idx = i;
    } 
     
  }

  return maxs;
}



//
//  get the score threshold out of a WM (log or not)
//
//  input:
//  	stars:		filter indicating if we compute for the position??
//  	transfo:	wm already background corrected???
//
//  output:
//  	minscore
//  	maxscore
//
float getWMThreshold(float** wm, int w, int* stars, float* bkg, int transfo, char** sites, int n, float c, float* minscore, float* maxscore) {

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

	*minscore = 10000.0;				
	*maxscore = -10000.0;
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
				score += wm[j][ N[ (int)(sites[i][j]) ] ]; 			// - bkg[ N[ (int)(sites[i][j]) ] ];
				//printf("j=%d %c mult pby %4.3f (%c) - %4.3f\n", j, sites[i][j], wm[j][ N[sites[i][j]] ], sites[i][j], bkg[ N[ sites[i][j]] ]);

			}
		}

		// min max score
		if (score > *maxscore)
			*maxscore = score;
		if (score < *minscore)
			*minscore = score;

		///printf("%s\t%4.3f\n", sites[i], score);

		SX  += score;
		SX2 += score * score;

		scores[i] = score;

	}

	AVG =  SX/n;
	STD =  (float)sqrt(( SX2 - SX*SX / n )  / (n - 1));
	//printf("AVG = %3.2f, S = %3.2f\n", AVG, STD);

	t = AVG - c * STD;

	cnt_above = 0;
	for (i=0; i<n; i++) {
		if (scores[i] >= t) {
			cnt_above ++;
		}
	}

	//printf("%d/%d sites score better than t\n", cnt_above, n);

	return t;
}





//
//  write WM scores to a file for distribution
//
void writeWMScoresToFile(float** wm, int w, int* stars, float* bkg, int transfo, char** sites, int n, char* file) 
{
  
  float score; 
  int   i, j;
  //char* T    = "ACTG";
  //float SX   = 0.0;
  //float SX2  = 0.0;
  //float AVG;
  //float STD  = 0.0;
  //float t;
  float* scores;
  //int   cnt_above;
  FILE* fp;

  fp = fopen(file, "w");
  if (!fp) {
    die("Cannot open score file for writing.\n");
  }

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

    fprintf(fp, "%s\t%4.3f\n", sites[i], score);
    

  }
  
  fclose(fp);

  return;

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


// transform a dinurecleotide count vector into a P(Xi|Xi-1) matrix to be used for scoring.
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

void shuffleRegexpMotif(char* motif, int gap, char** shuffled_motif)
{
  int  nc = 0;
  int* v  = 0;
  int  l,i;
  char** a_re;
  char*  newmotif;
  int  lg = 0, rg = 0;
  int  lnc;

  // motif length
  nc  = getRegexpMotifLength(motif, gap); 
  lnc = nc;

  //
  // printf("Motif has %d characters.\n", nc);
  // getchar();
  //

  l  = strlen(motif);

  getCharArrayFromRegexpMotif(motif, gap, &a_re); 

  while (a_re[lg][0] == '.') {
    lg++;
  }

  i = nc-1;
  while (a_re[i][0] == '.') {
    i--;
    rg++;
  }
  

  //printf("%d gaps to the left, %d to the right.\n", lg, rg);

  //for (i=0; i<nc; i++) {
  //  printf(" %s\n", a_re[i]);
  //}
  
  nc = nc - lg - rg;

  v  = randPermOfIndices(nc); 

  newmotif = (char*)calloc(l+1, sizeof(char));
  for (i=0; i<lg; i++) 
    strcat(newmotif, ".");

  for (i=0; i<nc; i++) {
    strcat(newmotif, a_re[ lg + v[i] ]);
  }

  for (i=0; i<rg; i++) 
    strcat(newmotif, ".");
  
  newmotif[l] = '\0';
  free(v);

  for (i=0; i<lnc; i++)
    free(a_re[i]);

  free(a_re);

  *shuffled_motif = newmotif;
}


void getCharArrayFromRegexpMotif(char* motif, int gap, char*** a_re) 
{
  

  int      l = strlen(motif);
  int      i = 0;  // cnt in ungapped motif
  int      j = 0;  // cnt in   gapped motif
  int      k = 0;  // cnt characters
  int      p_b;
  int      nc;
  char**   my_a_re = 0;
  char*    str;
  int      m = 0;

  nc = getRegexpMotifLength(motif, gap); 

  my_a_re = (char**)malloc(nc * sizeof(char*));

  while (i < l) {

    str = (char*)calloc(5, sizeof(char));
    m   = 0;

    if (motif[i] == '[') {

      p_b = i;
      while (motif[p_b] != ']') {	
	p_b++; 
      } 


      while (i <= p_b) {	
	str[m] = motif[i];
	i++; j++; m++;
      }


      
    } else {
      str[m] = motif[i];
      j++; i++; m++;
    }  

    str[m] = '\0';

    //printf("str=%s\n", str);
    my_a_re[k] = str;
    
    k++;
  }
  
  *a_re = my_a_re;
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

int BLOSUM62(char aa1, char aa2) {
  int** a;
  int i;
  a = (int**)malloc(100 * sizeof(int*));
  for (i=0; i<100; i++) {
    a[i] = (int*)malloc(100 * sizeof(int));
  }

  /* A => ... */
  a[65][65] = 4;	a[65][82] = -1;	a[65][78] = -2;	a[65][68] = -2;	a[65][67] = 0;	a[65][81] = -1;	a[65][69] = -1;	a[65][71] = 0;	a[65][72] = -2;	a[65][73] = -1;	a[65][76] = -1;	a[65][75] = -1;	a[65][77] = -1;	a[65][70] = -2;	a[65][80] = -1;	a[65][83] = 1;	a[65][84] = 0;	a[65][87] = -3;	a[65][89] = -2;	a[65][86] = 0;	a[65][66] = -2;	a[65][74] = -1;	a[65][90] = -1;	a[65][88] = -1;	a[65][42] = -4;	
  /* R => ... */
  a[82][65] = -1;	a[82][82] = 5;	a[82][78] = 0;	a[82][68] = -2;	a[82][67] = -3;	a[82][81] = 1;	a[82][69] = 0;	a[82][71] = -2;	a[82][72] = 0;	a[82][73] = -3;	a[82][76] = -2;	a[82][75] = 2;	a[82][77] = -1;	a[82][70] = -3;	a[82][80] = -2;	a[82][83] = -1;	a[82][84] = -1;	a[82][87] = -3;	a[82][89] = -2;	a[82][86] = -3;	a[82][66] = -1;	a[82][74] = -2;	a[82][90] = 0;	a[82][88] = -1;	a[82][42] = -4;	
  /* N => ... */
  a[78][65] = -2;	a[78][82] = 0;	a[78][78] = 6;	a[78][68] = 1;	a[78][67] = -3;	a[78][81] = 0;	a[78][69] = 0;	a[78][71] = 0;	a[78][72] = 1;	a[78][73] = -3;	a[78][76] = -3;	a[78][75] = 0;	a[78][77] = -2;	a[78][70] = -3;	a[78][80] = -2;	a[78][83] = 1;	a[78][84] = 0;	a[78][87] = -4;	a[78][89] = -2;	a[78][86] = -3;	a[78][66] = 4;	a[78][74] = -3;	a[78][90] = 0;	a[78][88] = -1;	a[78][42] = -4;	
  /* D => ... */
  a[68][65] = -2;	a[68][82] = -2;	a[68][78] = 1;	a[68][68] = 6;	a[68][67] = -3;	a[68][81] = 0;	a[68][69] = 2;	a[68][71] = -1;	a[68][72] = -1;	a[68][73] = -3;	a[68][76] = -4;	a[68][75] = -1;	a[68][77] = -3;	a[68][70] = -3;	a[68][80] = -1;	a[68][83] = 0;	a[68][84] = -1;	a[68][87] = -4;	a[68][89] = -3;	a[68][86] = -3;	a[68][66] = 4;	a[68][74] = -3;	a[68][90] = 1;	a[68][88] = -1;	a[68][42] = -4;	
  /* C => ... */
  a[67][65] = 0;	a[67][82] = -3;	a[67][78] = -3;	a[67][68] = -3;	a[67][67] = 9;	a[67][81] = -3;	a[67][69] = -4;	a[67][71] = -3;	a[67][72] = -3;	a[67][73] = -1;	a[67][76] = -1;	a[67][75] = -3;	a[67][77] = -1;	a[67][70] = -2;	a[67][80] = -3;	a[67][83] = -1;	a[67][84] = -1;	a[67][87] = -2;	a[67][89] = -2;	a[67][86] = -1;	a[67][66] = -3;	a[67][74] = -1;	a[67][90] = -3;	a[67][88] = -1;	a[67][42] = -4;	
  /* Q => ... */
  a[81][65] = -1;	a[81][82] = 1;	a[81][78] = 0;	a[81][68] = 0;	a[81][67] = -3;	a[81][81] = 5;	a[81][69] = 2;	a[81][71] = -2;	a[81][72] = 0;	a[81][73] = -3;	a[81][76] = -2;	a[81][75] = 1;	a[81][77] = 0;	a[81][70] = -3;	a[81][80] = -1;	a[81][83] = 0;	a[81][84] = -1;	a[81][87] = -2;	a[81][89] = -1;	a[81][86] = -2;	a[81][66] = 0;	a[81][74] = -2;	a[81][90] = 4;	a[81][88] = -1;	a[81][42] = -4;	
  /* E => ... */
  a[69][65] = -1;	a[69][82] = 0;	a[69][78] = 0;	a[69][68] = 2;	a[69][67] = -4;	a[69][81] = 2;	a[69][69] = 5;	a[69][71] = -2;	a[69][72] = 0;	a[69][73] = -3;	a[69][76] = -3;	a[69][75] = 1;	a[69][77] = -2;	a[69][70] = -3;	a[69][80] = -1;	a[69][83] = 0;	a[69][84] = -1;	a[69][87] = -3;	a[69][89] = -2;	a[69][86] = -2;	a[69][66] = 1;	a[69][74] = -3;	a[69][90] = 4;	a[69][88] = -1;	a[69][42] = -4;	
  /* G => ... */
  a[71][65] = 0;	a[71][82] = -2;	a[71][78] = 0;	a[71][68] = -1;	a[71][67] = -3;	a[71][81] = -2;	a[71][69] = -2;	a[71][71] = 6;	a[71][72] = -2;	a[71][73] = -4;	a[71][76] = -4;	a[71][75] = -2;	a[71][77] = -3;	a[71][70] = -3;	a[71][80] = -2;	a[71][83] = 0;	a[71][84] = -2;	a[71][87] = -2;	a[71][89] = -3;	a[71][86] = -3;	a[71][66] = -1;	a[71][74] = -4;	a[71][90] = -2;	a[71][88] = -1;	a[71][42] = -4;	
  /* H => ... */
  a[72][65] = -2;	a[72][82] = 0;	a[72][78] = 1;	a[72][68] = -1;	a[72][67] = -3;	a[72][81] = 0;	a[72][69] = 0;	a[72][71] = -2;	a[72][72] = 8;	a[72][73] = -3;	a[72][76] = -3;	a[72][75] = -1;	a[72][77] = -2;	a[72][70] = -1;	a[72][80] = -2;	a[72][83] = -1;	a[72][84] = -2;	a[72][87] = -2;	a[72][89] = 2;	a[72][86] = -3;	a[72][66] = 0;	a[72][74] = -3;	a[72][90] = 0;	a[72][88] = -1;	a[72][42] = -4;	
  /* I => ... */
  a[73][65] = -1;	a[73][82] = -3;	a[73][78] = -3;	a[73][68] = -3;	a[73][67] = -1;	a[73][81] = -3;	a[73][69] = -3;	a[73][71] = -4;	a[73][72] = -3;	a[73][73] = 4;	a[73][76] = 2;	a[73][75] = -3;	a[73][77] = 1;	a[73][70] = 0;	a[73][80] = -3;	a[73][83] = -2;	a[73][84] = -1;	a[73][87] = -3;	a[73][89] = -1;	a[73][86] = 3;	a[73][66] = -3;	a[73][74] = 3;	a[73][90] = -3;	a[73][88] = -1;	a[73][42] = -4;	
  /* L => ... */
  a[76][65] = -1;	a[76][82] = -2;	a[76][78] = -3;	a[76][68] = -4;	a[76][67] = -1;	a[76][81] = -2;	a[76][69] = -3;	a[76][71] = -4;	a[76][72] = -3;	a[76][73] = 2;	a[76][76] = 4;	a[76][75] = -2;	a[76][77] = 2;	a[76][70] = 0;	a[76][80] = -3;	a[76][83] = -2;	a[76][84] = -1;	a[76][87] = -2;	a[76][89] = -1;	a[76][86] = 1;	a[76][66] = -4;	a[76][74] = 3;	a[76][90] = -3;	a[76][88] = -1;	a[76][42] = -4;	
  /* K => ... */
  a[75][65] = -1;	a[75][82] = 2;	a[75][78] = 0;	a[75][68] = -1;	a[75][67] = -3;	a[75][81] = 1;	a[75][69] = 1;	a[75][71] = -2;	a[75][72] = -1;	a[75][73] = -3;	a[75][76] = -2;	a[75][75] = 5;	a[75][77] = -1;	a[75][70] = -3;	a[75][80] = -1;	a[75][83] = 0;	a[75][84] = -1;	a[75][87] = -3;	a[75][89] = -2;	a[75][86] = -2;	a[75][66] = 0;	a[75][74] = -3;	a[75][90] = 1;	a[75][88] = -1;	a[75][42] = -4;	
  /* M => ... */
  a[77][65] = -1;	a[77][82] = -1;	a[77][78] = -2;	a[77][68] = -3;	a[77][67] = -1;	a[77][81] = 0;	a[77][69] = -2;	a[77][71] = -3;	a[77][72] = -2;	a[77][73] = 1;	a[77][76] = 2;	a[77][75] = -1;	a[77][77] = 5;	a[77][70] = 0;	a[77][80] = -2;	a[77][83] = -1;	a[77][84] = -1;	a[77][87] = -1;	a[77][89] = -1;	a[77][86] = 1;	a[77][66] = -3;	a[77][74] = 2;	a[77][90] = -1;	a[77][88] = -1;	a[77][42] = -4;	
  /* F => ... */
  a[70][65] = -2;	a[70][82] = -3;	a[70][78] = -3;	a[70][68] = -3;	a[70][67] = -2;	a[70][81] = -3;	a[70][69] = -3;	a[70][71] = -3;	a[70][72] = -1;	a[70][73] = 0;	a[70][76] = 0;	a[70][75] = -3;	a[70][77] = 0;	a[70][70] = 6;	a[70][80] = -4;	a[70][83] = -2;	a[70][84] = -2;	a[70][87] = 1;	a[70][89] = 3;	a[70][86] = -1;	a[70][66] = -3;	a[70][74] = 0;	a[70][90] = -3;	a[70][88] = -1;	a[70][42] = -4;	
  /* P => ... */
  a[80][65] = -1;	a[80][82] = -2;	a[80][78] = -2;	a[80][68] = -1;	a[80][67] = -3;	a[80][81] = -1;	a[80][69] = -1;	a[80][71] = -2;	a[80][72] = -2;	a[80][73] = -3;	a[80][76] = -3;	a[80][75] = -1;	a[80][77] = -2;	a[80][70] = -4;	a[80][80] = 7;	a[80][83] = -1;	a[80][84] = -1;	a[80][87] = -4;	a[80][89] = -3;	a[80][86] = -2;	a[80][66] = -2;	a[80][74] = -3;	a[80][90] = -1;	a[80][88] = -1;	a[80][42] = -4;	
  /* S => ... */
  a[83][65] = 1;	a[83][82] = -1;	a[83][78] = 1;	a[83][68] = 0;	a[83][67] = -1;	a[83][81] = 0;	a[83][69] = 0;	a[83][71] = 0;	a[83][72] = -1;	a[83][73] = -2;	a[83][76] = -2;	a[83][75] = 0;	a[83][77] = -1;	a[83][70] = -2;	a[83][80] = -1;	a[83][83] = 4;	a[83][84] = 1;	a[83][87] = -3;	a[83][89] = -2;	a[83][86] = -2;	a[83][66] = 0;	a[83][74] = -2;	a[83][90] = 0;	a[83][88] = -1;	a[83][42] = -4;	
  /* T => ... */
  a[84][65] = 0;	a[84][82] = -1;	a[84][78] = 0;	a[84][68] = -1;	a[84][67] = -1;	a[84][81] = -1;	a[84][69] = -1;	a[84][71] = -2;	a[84][72] = -2;	a[84][73] = -1;	a[84][76] = -1;	a[84][75] = -1;	a[84][77] = -1;	a[84][70] = -2;	a[84][80] = -1;	a[84][83] = 1;	a[84][84] = 5;	a[84][87] = -2;	a[84][89] = -2;	a[84][86] = 0;	a[84][66] = -1;	a[84][74] = -1;	a[84][90] = -1;	a[84][88] = -1;	a[84][42] = -4;	
  /* W => ... */
  a[87][65] = -3;	a[87][82] = -3;	a[87][78] = -4;	a[87][68] = -4;	a[87][67] = -2;	a[87][81] = -2;	a[87][69] = -3;	a[87][71] = -2;	a[87][72] = -2;	a[87][73] = -3;	a[87][76] = -2;	a[87][75] = -3;	a[87][77] = -1;	a[87][70] = 1;	a[87][80] = -4;	a[87][83] = -3;	a[87][84] = -2;	a[87][87] = 11;	a[87][89] = 2;	a[87][86] = -3;	a[87][66] = -4;	a[87][74] = -2;	a[87][90] = -2;	a[87][88] = -1;	a[87][42] = -4;	
  /* Y => ... */
  a[89][65] = -2;	a[89][82] = -2;	a[89][78] = -2;	a[89][68] = -3;	a[89][67] = -2;	a[89][81] = -1;	a[89][69] = -2;	a[89][71] = -3;	a[89][72] = 2;	a[89][73] = -1;	a[89][76] = -1;	a[89][75] = -2;	a[89][77] = -1;	a[89][70] = 3;	a[89][80] = -3;	a[89][83] = -2;	a[89][84] = -2;	a[89][87] = 2;	a[89][89] = 7;	a[89][86] = -1;	a[89][66] = -3;	a[89][74] = -1;	a[89][90] = -2;	a[89][88] = -1;	a[89][42] = -4;	
  /* V => ... */
a[86][65] = 0;	a[86][82] = -3;	a[86][78] = -3;	a[86][68] = -3;	a[86][67] = -1;	a[86][81] = -2;	a[86][69] = -2;	a[86][71] = -3;	a[86][72] = -3;	a[86][73] = 3;	a[86][76] = 1;	a[86][75] = -2;	a[86][77] = 1;	a[86][70] = -1;	a[86][80] = -2;	a[86][83] = -2;	a[86][84] = 0;	a[86][87] = -3;	a[86][89] = -1;	a[86][86] = 4;	a[86][66] = -3;	a[86][74] = 2;	a[86][90] = -2;	a[86][88] = -1;	a[86][42] = -4;	
/* B => ... */
 a[66][65] = -2;	a[66][82] = -1;	a[66][78] = 4;	a[66][68] = 4;	a[66][67] = -3;	a[66][81] = 0;	a[66][69] = 1;	a[66][71] = -1;	a[66][72] = 0;	a[66][73] = -3;	a[66][76] = -4;	a[66][75] = 0;	a[66][77] = -3;	a[66][70] = -3;	a[66][80] = -2;	a[66][83] = 0;	a[66][84] = -1;	a[66][87] = -4;	a[66][89] = -3;	a[66][86] = -3;	a[66][66] = 4;	a[66][74] = -3;	a[66][90] = 0;	a[66][88] = -1;	a[66][42] = -4;	
 /* J => ... */
 a[74][65] = -1;	a[74][82] = -2;	a[74][78] = -3;	a[74][68] = -3;	a[74][67] = -1;	a[74][81] = -2;	a[74][69] = -3;	a[74][71] = -4;	a[74][72] = -3;	a[74][73] = 3;	a[74][76] = 3;	a[74][75] = -3;	a[74][77] = 2;	a[74][70] = 0;	a[74][80] = -3;	a[74][83] = -2;	a[74][84] = -1;	a[74][87] = -2;	a[74][89] = -1;	a[74][86] = 2;	a[74][66] = -3;	a[74][74] = 3;	a[74][90] = -3;	a[74][88] = -1;	a[74][42] = -4;	
 /* Z => ... */
 a[90][65] = -1;	a[90][82] = 0;	a[90][78] = 0;	a[90][68] = 1;	a[90][67] = -3;	a[90][81] = 4;	a[90][69] = 4;	a[90][71] = -2;	a[90][72] = 0;	a[90][73] = -3;	a[90][76] = -3;	a[90][75] = 1;	a[90][77] = -1;	a[90][70] = -3;	a[90][80] = -1;	a[90][83] = 0;	a[90][84] = -1;	a[90][87] = -2;	a[90][89] = -2;	a[90][86] = -2;	a[90][66] = 0;	a[90][74] = -3;	a[90][90] = 4;	a[90][88] = -1;	a[90][42] = -4;	
 /* X => ... */
 a[88][65] = -1;	a[88][82] = -1;	a[88][78] = -1;	a[88][68] = -1;	a[88][67] = -1;	a[88][81] = -1;	a[88][69] = -1;	a[88][71] = -1;	a[88][72] = -1;	a[88][73] = -1;	a[88][76] = -1;	a[88][75] = -1;	a[88][77] = -1;	a[88][70] = -1;	a[88][80] = -1;	a[88][83] = -1;	a[88][84] = -1;	a[88][87] = -1;	a[88][89] = -1;	a[88][86] = -1;	a[88][66] = -1;	a[88][74] = -1;	a[88][90] = -1;	a[88][88] = -1;	a[88][42] = -4;	
/* * => ... */
 a[42][65] = -4;	a[42][82] = -4;	a[42][78] = -4;	a[42][68] = -4;	a[42][67] = -4;	a[42][81] = -4;	a[42][69] = -4;	a[42][71] = -4;	a[42][72] = -4;	a[42][73] = -4;	a[42][76] = -4;	a[42][75] = -4;	a[42][77] = -4;	a[42][70] = -4;	a[42][80] = -4;	a[42][83] = -4;	a[42][84] = -4;	a[42][87] = -4;	a[42][89] = -4;	a[42][86] = -4;	a[42][66] = -4;	a[42][74] = -4;	a[42][90] = -4;	a[42][88] = -4;	a[42][42] = 1;	
 
 int b = a[(int)aa1][(int)aa2];

 for (i=0; i<100; i++) {
   free(a[i]);
 }
 free(a); 
 
 return b;
}
