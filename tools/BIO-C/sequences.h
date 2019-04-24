#include <search.h>
#include "hashtable.h"
#include "uthash.h"



typedef struct _ScanACE_match {
  
  float score;
  int   position;
  char* gene;
  int   strand;

} ScanACE_match;


// 
// hidden variables. dangerous!
//
int* C;
int* N;
char* char_compl;

int  ntoC[6];
char ntoC1[16][3];

int** N1;

void readChrData(char* file, char*** chrnames, int** chrlens, int* numchroms);

void readSNPfile(char* file, char** a_snp, char*** snp_names);

void loadFastaSequencesAndMakeIndex(char* fastafile, char*** seqs, char*** seqn, int* numgenes,  struct my_hsearch_data** hash_genes);

int  basicEditDistance(char* ss1, char* ss2); 

void getWMConsensus(float** wm, int w, char** consensus);


void get_hg18_chroms(char*** chroms, int* numchroms, int** chrlens);

void writeWMScoresToFile(float** wm, int w, int* stars, float* bkg, int transfo, char** sites, int n, char* file);

double aaConsIndex(char* str);

int   chr2index(char* chr);
char* strrev (char* str); 
int   numSequencesInFastaFile(char* file); 
char* IUPAC2nt(char c);
int   BLOSUM62(char aa1, char aa2);

int   rand_location_around_pos_in_chrlen(int pos, int range, int chrlen);

int   rand_location_around_pos(int pos, int range, char* c);
int   rand_chr_location(char* c);

void  readSNPfile(char* file, char** a_snp, char*** snp_names);

void  seq2lower(char* s);
void  seqlower2dots(char* s); 

void  generateManyRandomSequences(float** wm, int w, int m, char*** seqs);
void  generateRandomSequence(float** wm, int w, char** s); 

void  sw(char* ss1, char* ss2, int R, int P, int O, int G, float* id, int* alnlen, int* qst, int* qen, int* sst, int* sen, int* sco, char** aln1, char** aln2); 

void  intWMtofloatWM(int** intwm, int w, int nsites, float*** fwm); 
void intWMtofloatWM_J(int** intwm, int w, int* nsites, float*** fwm) ;

void  WMtologWMbkg(float** fwm, int w, float* rawbkg, float*** logwm); 
void  read_BULYK_WM(char* wmfile, float*** wm, int* w);
 
int   hg18_chrlen(char* c); 

int   rand_chr_interval(char* c, int l);

long  sequencesOverlap(long s1, long e1, long s2, long e2);

void  read_JASPAR_WM(char* wmfile, int*** wm, int* w, int** nsites);
void read_JASPAR_WM_OLD(char* wmfile, int*** wm, int* w, int* nsites_j);

void  readScanACEMatches(char* file, ScanACE_match** scm, int* nm); 

void  shuffleRegexpMotif(char* motif, int gap, char** shuffled_motif);
void  getCharArrayFromRegexpMotif(char* motif, int gap, char*** a_re); 

char* getGappedKmer(char* kmer, int gap); 
char* getGappedMotif(char* motif, int motifsize, int gap); 
int   getRegexpMotifLength(char* motif, int gap); 
void  getIntegerWMfromRegexp(char* motif, int gap, int*** wm, int* w);
void  getFirstOrderMarkovIntegerWMfromRegexp(char* motif, int gap, int*** wm, int* w);
float getScoreOnSeq_1M(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn); 
void  convert_wm_0m_to_1m(int** intwm_0m, int w, int*** intwm_1m);
void  printIntegerWM(int** wm, int w);
void  printIntegerWM_1M(int** wm, int w);
void printBckd(float* m);

void  integerWMtoLog(int** wm1, float*** wm2, int h, int w); 
void  integerWMtoACE(int** intwm, int w, char*** m); 
void  findAllSeqsMaxWMScores(char** seqs, int numseqs, float** logwm, int mt, int w, int* stars, float* bkg, int rna, float** scores, int** a_idx);
float findMaxWMScore(char* seq, float** logwm, int mt, int w, int* stars, float* bkg, int rna, int* idx);

void  initialize_nt(); 
float getWMThreshold(float** wm, int w, int* stars, float* bkg, int log, char** sites, int n, float c, float* minscore, float* maxscore);
float getScoreOnSeq(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn); 
//void  readWM(char* wmfile, float*** wm, int* w, char*** sites, int* n, int** stars, int transfo, float* bkgraw);
void  readACEintWM(char* wmfile, int*** wm, int* w, char*** sites, int* n, int** stars);
void  ACEintWMtologWM(int** intwm, int w, int nsites, float* rawbkg, float*** logwm);
void ACEintWMtologWM_J(int** intwm, int w, int* nsites, float* rawbkg, float*** logwm);

void  findAllWeightMatrixMatches(float** wml, int w, int* stars, float* bkgl, float t, char* seq, int rna, int** matches_pos, int* num_matches, char** matches_ori, float** matches_sco, int max_num_matches);
void findAllWeightMatrixMatchesInIntervals(float** wml, int w, int* stars, float* bkgl, float t, char* seq, int rna, int** matches_pos, int* num_matches, char** matches_ori, float** matches_sco, int max_num_matches, int** intervals, int numint);

void  printWM(float** wm, int w);
void  dint_to_scoringWM(int** diwm, int w, int*** scwm); 





