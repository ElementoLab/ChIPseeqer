int* C;
int* N;
char* char_compl;

int  ntoC[6];
char ntoC1[16][3];

int** N1;



char* getGappedKmer(char* kmer, int gap); 
char* getGappedMotif(char* motif, int motifsize, int gap); 
int   getRegexpMotifLength(char* motif, int gap); 
void  getIntegerWMfromRegexp(char* motif, int gap, int*** wm, int* w);
void  getFirstOrderMarkovIntegerWMfromRegexp(char* motif, int gap, int*** wm, int* w);
float getScoreOnSeq_1M(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn); 
void  convert_wm_0m_to_1m(int** intwm_0m, int w, int*** intwm_1m);
void  printIntegerWM(int** wm, int w);
void  printIntegerWM_1M(int** wm, int w);
void  integerWMtoLog(int** wm1, float*** wm2, int h, int w); 
void  integerWMtoACE(int** intwm, int w, char*** m); 
void  findAllSeqsMaxWMScores(char** seqs, int numseqs, float** logwm, int mt, int w, int* stars, float* bkg, int rna, float** scores, int** a_idx);
float findMaxWMScore(char* seq, float** logwm, int mt, int w, int* stars, float* bkg, int rna, int* idx);

void  initialize_nt(); 
float getWMThreshold(float** wm, int w, int* stars, float* bkg, int log, char** sites, int n);
float getScoreOnSeq(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn); 
//void  readWM(char* wmfile, float*** wm, int* w, char*** sites, int* n, int** stars, int transfo, float* bkgraw);
void  readACEintWM(char* wmfile, int*** wm, int* w, char*** sites, int* n, int** stars);
void  ACEintWMtologWM(int** intwm, int w, int nsites, float* rawbkg, float*** logwm);

void  findAllWeightMatrixMatches(float** wml, int w, int* stars, float* bkgl, float t, char* seq, int rna, int** matches_pos, int* num_matches, char** matches_ori, int max_num_matches);

void  printWM(float** wm, int w);
void  dint_to_scoringWM(int** diwm, int w, int*** scwm); 
float compare_motifs(float *m1, int l1, float* m2, int l2, int ss, int min_inf_cols); 
void getLinearFloatWMfromRegexp(char* motif, float** lwm, int* lw);




