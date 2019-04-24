// CS
#define ELAND 0
#define MIT 1
#define BED 2
#define SAM 3
#define BOWTIESAM 4
#define BAM 5


typedef struct _ScanACE_match {
  
  float score;
  int   position;
  char* gene;
  int   strand;

} ScanACE_match;

int* C;
int* N;
char* char_compl;

int  ntoC[6];
char ntoC1[16][3];

int** N1;
int formatToId(char* formattext);

int readMapData(char* file, unsigned char** map, int t);
int chr2index(char* chr) ;

int SAMbitflagNum(int bitflag, int num);
int SAMbitflag_strand(int bitflag);
int SAMbitflag_isPairedProper(int bitflag);
int SAMbitflag_ismapped(int bitflag); 

void readChrLen(char* file, int** chrlens);
void readChrData(char* file, char*** chrnames, int** chrlens, float** chrmap, int* numchroms);
void readChrData2(char* file1, char* file2, int hasid1, int hasid2, char*** chrnames, int* numchroms);

long int CountAlignedNucleotides(char* file, int format);

int CountNonClonalReads(char* file, int format, int chrlen);

int CountReads(char* file, int format);

int hg18_chrlen(char* c); 
void get_hg18_chroms(char*** chroms, int* numchroms, int** chrlens);

int hg19_chrlen(char* c); 
void getCountFromReadFile(char* file, int format, char* chr, long chrlen, int readlen, int fraglen, int uniquereads, int uniqmap, unsigned short** counts, int* numreads, int* numclonalreads);
void getCountFromReadFileOLD(char* file, int format, char* chr, long chrlen, int readlen, int fraglen, unsigned short** counts, int* numreads);
float hg18_30mer_mappability(char* c);
float hg19_30mer_mappability(char* c);

void generateRandomSequence(float** wm, int w, char** s); 


int rand_chr_interval(char* c, int l);
int rand_chr_interval2(char* c, int l, int m);
int rand_chr_interval_startpoint_around_interval(char* c, int st, int en, int d);

long sequencesOverlap(long s1, long e1, long s2, long e2);
long sequencesDistance(long s1, long e1, long s2, long e2);

void read_JASPAR_WM(char* wmfile, int*** wm, int* w, int* nsites);

void readScanACEMatches(char* file, ScanACE_match** scm, int* nm); 

void shuffleRegexpMotif(char* motif, int gap, char** shuffled_motif);
void getCharArrayFromRegexpMotif(char* motif, int gap, char*** a_re); 

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
float getWMThreshold(float** wm, int w, int* stars, float* bkg, int log, char** sites, int n, float c);
float getScoreOnSeq(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn); 
//void  readWM(char* wmfile, float*** wm, int* w, char*** sites, int* n, int** stars, int transfo, float* bkgraw);
void  readACEintWM(char* wmfile, int*** wm, int* w, char*** sites, int* n, int** stars);
void  ACEintWMtologWM(int** intwm, int w, int nsites, float* rawbkg, float*** logwm);

void  findAllWeightMatrixMatches(float** wml, int w, int* stars, float* bkgl, float t, char* seq, int rna, int** matches_pos, int* num_matches, char** matches_ori, float** matches_sco, int max_num_matches);

void  printWM(float** wm, int w);
void  dint_to_scoringWM(int** diwm, int w, int*** scwm); 





