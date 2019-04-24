#ifndef MI_LIBRARY_H
#define MI_LIBRARY_H

#ifndef NAN
#define NAN (0.0 / 0.0)
#endif


#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)


typedef struct _Params {
  
  char* expfile;
  char* kmerfile;
  char* fastafile;
  int   kmersize;
  int   quantized;
  int   singlestrand;
  int   verbose;
  int   shuffle;
  int   shuffle_rank;
  int   mbins;
  char* outfile;
  int   report;
  int   mbins_dist;

} Params;


void encode_kmer(char* kmer, int kmersize, int add5, int add3, char** d_ch, char** ch, int nbchars, int** encoded_kmer, char** encoded_re);

char* getGappedKmer(char* kmer, int gap); 
void add_small_values_to_identical_floats(float* E, int n);
void bubbleSortFI(FloatAndIndex* fi, int n);


int max_rank_test(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int do_val, int* val, int do_z, double* z);
int jacknife_max_rank_test(int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int jn, int jn_f, int jn_t, int do_val, int* val); 
 

double get_zscore_and_rank_value(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int shuffle_rank, double* c); 
void get_rank_and_zscore(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int* rank, double* zscore); 

float calc_gc_content(char* seq, int l);
float calc_CpG_content(char* seq, int l);

void get_Params(int argc, char** argv, Params* p); 

void readKmers_general (char* kmerfile, int max_nb_kmers, int inc_nb_kmers, char*** kmers, int* nbkmers, int* kmersize); 
void readKmers_general_special_optim (char* kmerfile, int max_nb_kmers, int inc_nb_kmers, char*** motifs, char*** seeds, int* nbkmers, int* kmersize); 


void quantize_M_counts( short* M, int nborfs, int* mbins, int** M_q); 
void quantize_M_zero_eqpop( float* M, int nborfs, int mbins, int** M_q); 
void add_to_report(FILE* fpr, char* kmer, int kmersize, int gap, int* M_q, int mbins, int* E_q, int ebins, float* E_q_bins, int nborfs);
double get_zscore(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle); 
double get_zscore_and_do_interval(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, int nbkmers, double* c); 

void quantize_E(float* E, int nborfs, int quantized, int* ebins, int** E_q, float** E_q_bins); 

void quantize_M( float* M, int nborfs, int mbinary, int m3bins, int* mbins, int** M_q); 
char** allocate_more_kmers(char** kmers, int n, int max_nb_kmers, int inc_nb_kmers); 
void readKmers (char* kmer, char* kmerfile, int max_nb_kmers, int inc_nb_kmers, char*** kmers, float** oldmis, int* nbkmers, int* kmersize); 

#endif
