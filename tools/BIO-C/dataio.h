#ifndef DATAIO_H
#define DATAIO_H



/*
 *  structure that keeps running info
 *
 */
typedef struct _seqI {
  
  FILE* fp;
  int   nextSequence_started; 
  int   nextSequence_ended;   
  char* nextSequence_currentLine;
  int   seqlen_inc;
  int   verbose;

  char* seq;
  unsigned long int   pos;
  unsigned long int   max_seqlen;
  unsigned long int   fread_chunksize;
  char* chunk;

  char* name;
  unsigned long int pos_in_chunk;
  unsigned long int enter_new_seq;
  unsigned long int pos_in_seq;
  int inname;
  int inseq;

  int pos_in_name;
  unsigned long int chunk_len;
  
  int eof;
  int freadeof;
} seqI;


typedef struct _LineI {

  //int chumkum;
  int bufflen;
  FILE* fp;
  char* buff;
  //char* prevline;
  char* buffptr;
  int   eof;

} LineI;

char * readFile(char *);

void  LineIclose(LineI* li);
char* nextLine(LineI* li);
int   LineIopen(LineI* li, char* file);

void  readmiRNASeeds(char* seedfile, char*** seeds, char*** mirnames, int* numseeds);
int   file_exists (char * fileName);

char*          decodeDNAseq(unsigned char* c, int l); 
unsigned char* encodeDNAseq(char* s); 

void split_line_delim_sep(char* l, char* delim, char*** p, int* m); 

void split_line_delim(char* l, char* delim, char*** p, int* m);

char char_complement(char s); 
int nbSequencesInFastaFile(char* file); 

int seqI_open(seqI* s, char* file);
char* seqI_nextSequence(seqI* s, char** name, int* size);
void seqI_close(seqI* s);
void seqI_set_seqlen_inc(seqI* s, int n);
void seqI_set_verbose(seqI* s, int v);

char* seqI_nextSequence_fast(seqI* s, char** name, int* size);
char* seqI_nextSequence_fast_light(seqI* s, char** name, int* size);

int seqI_open_fast_light(seqI* s, char* file);

int seqI_open_fast(seqI* s, char* file);
void seqI_set_max_seqlen(seqI* s, unsigned long int n);
void seqI_set_fread_chunk(seqI* s, unsigned long int n);

int getNumPositiveEntries(int* v, int n);

unsigned char* create_binarized_array(int n);
void     set_entry_in_binarized_array(unsigned char* p, int i); 
int      get_entry_in_binarized_array(unsigned char* p, int i);
void     unset_entry_in_binarized_array(unsigned char* p, int i); 
int*     c_matrix_column_binarized      (unsigned char** data, int j, int n); 

void readAASequences(char* filename, int* m, int* n, int*** data, char*** rownames, char** symbols, int* nbsymbols);
void readAlignACEMotif(char* filename, int* m, int* n, int*** data); 
int data_nucl(char c);

int exist_parameter(int argc, char** argv, char* param);
int exist_flag(int argc, char** argv, char* param);


char*** initialize_aacode(); 
char*   substr(char *string, int start, int length); 
char*   translate(char* seq, char*** aacode);
int     indexOfKmer(char* kmer); 
int     powint(int x, int y); 
char    codon2aa(char* codon); 
void    readFloatProfile(char* filename, int* n, float** data, char*** rownames); 
void    readFloatTable(char *filename, int *m, int *n, float ***data, char ***rownames, char ***colnames, int logt, int header);
void    readIntTable(char* filename, int* m, int* n, int*** data, char*** rownames, char*** colnames); 
void    readIntProfile(char* filename, int* n, int** data, char*** rownames); 
void    readStringTable(char* filename, char**** data, int* n, int* m); 
char**  readStringSet(char* setfile, int* n); 
float*  readFloatSet(char *setfile, int *n);
char*   mystrtok(char *s, char delim);
void    chomp(char *s);
char*   complement(char *s);
char*   mybasename(char *filename);
char*   get_parameter(int argc, char **argv, char *param);
void    die(char* s); 
char*   nextSequence(FILE* fp, char** name, int* size, int* nextSequence_started, int* nextSequence_ended, char* nextSequence_currentLine); 
void    uctransform(char* seq);
void    lctransform(char* seq);
char*   uc(char* seq); 
char*   lc(char* seq);
int     nbLinesInFile(char* file); 
void    showIntVector(int* v, int m); 
void    showFloatVector(float* v, int m); 
float*  itof_vector(int* v, int n);
int     nbLinesInFile(char* file); 
float*  f_matrix_column(float** data, int j, int n);
float*  stof_matrix_column(short** data, int j, int n) ;
float*  itof_matrix_column(int** data, int j, int n); 

short*  s_matrix_column(short** data, int j, int n); 

short*  ftos_vector(float* v, int n); 
int*   ctoi_vector(char* v, int n); 

int* stoi_matrix_column(short** data, int j, int n); 

int* i_matrix_column(int** data, int j, int n); 
int* c_matrix_column(char** data, int j, int n); 

void showTwoIntVectors(int* v1, int* v2, int m); 
 
#endif

