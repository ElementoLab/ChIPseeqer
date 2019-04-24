#ifndef DATAIO_H
#define DATAIO_H



/*
 *  structure that keeps running info
 *
 */
typedef struct _seqI {
  
  FILE* fp;
  int   nextSequence_started;  // boolean
  int   nextSequence_ended;    // boolean
  char* nextSequence_currentLine;
  
} seqI;

int seqI_open(seqI* s, char* file);
char* seqI_nextSequence(seqI* s, char** name, int* size);
void seqI_close(seqI* s);


unsigned char* create_binarized_array(int n);
void     set_entry_in_binarized_array(unsigned char* p, int i); 
int      get_entry_in_binarized_array(unsigned char* p, int i);
void     unset_entry_in_binarized_array(unsigned char* p, int i); 
int*     c_matrix_column_binarized      (unsigned char** data, int j, int n); 

void readAASequences(char* filename, int* m, int* n, int*** data, char*** rownames, char** symbols, int* nbsymbols);
int exist_parameter(int argc, char** argv, char* param);

char*** initialize_aacode(); 
char*   substr(char *string, int start, int length); 
char*   translate(char* seq, char*** aacode);
int     indexOfKmer(char* kmer); 
int     powint(int x, int y); 
char    codon2aa(char* codon); 


void    readFloatTable(char *filename, int *m, int *n, float ***data, char ***rownames, char ***colnames, int logt, int header);
void    readIntTable(char* filename, int* m, int* n, int*** data, char*** rownames, char*** colnames); 
void    readStringTable(char* filename, char**** data, int* n, int* m); 

float*  readFloatSet(char *setfile, int *n);
char*   mystrtok(char *s, char delim);
void    chomp(char *s);
char*   complement(char *s);
char*   mybasename(char *filename);

char  *get_parameter(int argc, char **argv, char *param);
void   die(char* s); 
char  *nextSequence(FILE* fp, char** name, int* size, int* nextSequence_started, int* nextSequence_ended, char* nextSequence_currentLine); 
char  *uc(char* seq); 
int    nbLinesInFile(char* file); 
void   showIntVector(int* v, int m); 
void showFloatVector(float* v, int m); 

float* itof_vector(int* v, int n);
int nbLinesInFile(char* file); 
float* f_matrix_column(float** data, int j, int n);
float* stof_matrix_column(short** data, int j, int n) ;
float* itof_matrix_column(int** data, int j, int n); 

short* s_matrix_column(short** data, int j, int n); 
short* ftos_vector(float* v, int n); 
int* stoi_matrix_column(short** data, int j, int n); 
int* i_matrix_column(int** data, int j, int n); 
int* c_matrix_column(char** data, int j, int n); 

 
#endif

