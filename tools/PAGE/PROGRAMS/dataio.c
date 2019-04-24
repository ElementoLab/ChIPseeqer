#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>

#ifdef USEMPATROL
#include <mpatrol.h>
#endif

#include "dataio.h"
#include "statistics.h"

int verbose;

#ifndef NAN
static const double NAN = (0.0 / 0.0);
#endif

unsigned char* create_binarized_array(int n) {

  unsigned char* p;
  int o = n / 8 + 1;
  
  p = (unsigned char*)calloc(o, sizeof(unsigned char));

  return p;
}



void unset_entry_in_binarized_array(unsigned char* p, int i) 
{

  int  o = i / 8;
  int  n = i % 8;

  p[o] = p[o] & ~(1 << n); 

}

void set_entry_in_binarized_array(unsigned char* p, int i) 
{

  int  o = i / 8;
  int  n = i % 8;

  p[o] = p[o] | (1 << n); 

}


int get_entry_in_binarized_array(unsigned char* p, int i) {
 
 
  int  o = i / 8;
  int  n = i % 8;

  return ((p[o] & (1 << n)) > 0); 

}




int* c_matrix_column_binarized(unsigned char** data, int j, int n) 
{
  int* v;
  int    i;

  v = (int*) malloc ( n * sizeof( int ));
  for (i=0; i<n; i++) {
    v[i] = get_entry_in_binarized_array(data[i], j); //(int)(data[i][j]);
  }
  
  return v;
}



void writeFloatSet(char* file, float* f, int n) 
{
  
  FILE* fp;
  int   i;

  fp = fopen(file, "w");
  if (!fp) {
    die("writeFloatSet: cannot open file for writing\n");
  }
  for (i=0; i<n; i++) {
    fprintf(fp, "%f\n", f[i]);
  }
  fclose(fp);

}


int nbLinesInFile(char* file) 
{
  
  char* buff;
  FILE* fp;
  int   n;

  buff = (char*)malloc(10000 * sizeof(char));
  fp = fopen(file, "r");
  if (!fp) {
    printf("could not open set data %s\n", file);
  } 
  fgets(buff, 100000, fp);
  n = 0;
  while (!feof(fp)) {
    n++;
    fgets(buff, 100000, fp);
    //n++;
  }  
  fclose(fp);

  free(buff);
  return n;
}


//
//  read in an aa sequence
//
void readAASequences(char* filename, int* m, int* n, int*** data, char*** rownames, char** symbol_table, int* nbsymbols) {
  
  char* s;
  int   mym;
  int   myn;
  char* buff;
  FILE* fp;
  int   i;
  int   mynmax = 100000;
  int   l = 0;
  int   lookup_symbol[256];
  int   nbsymb = 0;



  buff  = (char*)malloc(100000 * sizeof(char));

  for (i=0; i<256; i++) {
    lookup_symbol[ i ] = -1;
  }

  (*symbol_table) = (char*)calloc(256, sizeof(char));

  fp    = fopen(filename, "r");
  if (!fp) {
    die("readFloatTable: please enter a valid filename\n");
  }

  //
  // get the first line and estimate the number of columns
  //
  
  fgets(buff, 100000, fp);
  chomp(buff);
  mym = 0;

  //
  //   allocate a big chunk of data
  //
  *rownames = (char**)malloc(mynmax * sizeof(char*));
  
    
  
  //
  //   allocate a big chunk of data
  //
  *data     = (int**)malloc(mynmax * sizeof(int*));

  myn = 0;
  while (!feof(fp)) {
    fgets(buff, mynmax, fp);
    chomp(buff);

    if (feof(fp))
      break; 

    //printf("buff=%s\n", buff);
    
    //
    // allocate a line
    //
    

    //
    // get the row name 
    //
    s = mystrtok(buff, '\t');
    (*rownames)[myn] = strdup(s);
    free(s);

    //
    //  get the aa sequence
    //
    s = mystrtok(0, '\t');
    l = strlen(s);
    (*data)[myn] = (int*)malloc(l * sizeof(int));

    for (i=0; i<l; i++) {
      
      if (lookup_symbol[ (int)s[i] ] == -1) {
	// create a new symbol
	
	(*symbol_table)[ nbsymb ] = s[i];
	
	lookup_symbol[ (int)s[i] ] = nbsymb ++;

	
      }

      (*data)[myn][i] = lookup_symbol[ (int)s[i] ];
      

    }

    free(s);
    
    myn ++;
  
   
  }



  //printf("ehat ?\n");

  *n = myn;
  *m = l;
  *nbsymbols = nbsymb;

  fclose(fp);

}


void readFloatTable(char* filename, int* m, int* n, float*** data, char*** rownames, char*** colnames, int logt, int header) {
  
  char* s;
  int   mym;
  int   myn;
  char* buff;
  FILE* fp;
  int   i;
  int   mynmax = 100000;

  buff  = (char*)malloc(mynmax * sizeof(char));
  fp    = fopen(filename, "r");
  if (!fp) {
    printf("readFloatTable: please enter a valid filename (%s invalid)\n", filename);
    exit(0);
  }

  //
  // get the first line and estimate the number of columns
  //
  
  //getline(&buff, &len, fp);
  fgets(buff, mynmax, fp);
  chomp(buff);
  mym = 0;

  s = mystrtok(buff, '\t'); free(s);
  while ((s = mystrtok(0, '\t'))) {
    mym ++; free(s);
  }
  
  
  //
  // allocate the right number of columns
  //

  if (header == 1) {
    *colnames = (char**)malloc(mym * sizeof(char*));
    s   = mystrtok(buff, '\t'); free(s);
    i   = 0;
    while ((s = mystrtok(0, '\t'))) {
      (*colnames)[i] = strdup(s);  free(s);
      i++;
    }
  } else {
    rewind(fp);
  }

  //
  //   allocate a big chunk of data
  //
  *rownames = (char**)malloc(mynmax * sizeof(char*));
  *data     = (float**)malloc(mynmax * sizeof(float*));

  myn = 0;
  while (!feof(fp)) {

    fgets(buff, mynmax, fp);

    if (feof(fp))
      break; 

    chomp(buff);

    //printf("buff=%s\n", buff);
    
    //
    // allocate a line
    //
    (*data)[myn] = (float*)malloc(mym * sizeof(float));
 
    //
    // get the row name 
    //
    s = mystrtok(buff, '\t');
    (*rownames)[myn] = strdup(s);
    free(s);


    // get the values 
    i = 0;
    while ((s = mystrtok(0, '\t'))) {
      //printf("%d:%s\n", i, s);
      
      if ((strcmp(s, "") != 0) && (strcmp(s, "NA") != 0)) {
	(*data)[myn][i] = atof(s); //printf("%s\t%3.2f\n", s, (*data)[myn][i]);
      } else 
	(*data)[myn][i] = NAN;
      
      if ( logt == 1 ) {
	(*data)[myn][i] = log( (*data)[myn][i] );
      }

       free(s);
       i ++;
    }
    
    
    //s = mystrtok(buff, '\t');
    //(*colnames)[i] = strdup(s);
    //i++;
    //}
  
    //printf("m=%d\n", mym);
  
    //*m = mym;
    myn ++;

    if (myn == mynmax) {
      die("readFloatTable: running out of memory, please recompile ..\n"); 
    }

    #ifdef GENEMAX
    if (myn == GENEMAX) {
      break;
    }
    #endif

  }



  //printf("ehat ?\n");

  *n = myn;
  *m = mym;

  fclose(fp);

}



void readIntTable(char* filename, int* m, int* n, int*** data, char*** rownames, char*** colnames) 
{
  
  char* s;
  int   mym;
  int   myn;
  char* buff;
  FILE* fp;
  int   i;
  int   mynmax = 100000;

  buff  = (char*)malloc(100000 * sizeof(char));
  fp    = fopen(filename, "r");
  if (!fp) {
    die("readIntTable: please enter a valid filename\n");
  }

  //
  // get the first line and estimate the number of columns
  //
  
  fgets(buff, mynmax, fp);
  chomp(buff);
  mym = 0;

  s = mystrtok(buff, '\t'); free(s);
  while ((s = mystrtok(0, '\t'))) {
    mym ++; free(s);
  }
  
  
  //
  // allocate the right number of columns
  //
  *colnames = (char**)malloc(mym * sizeof(char*));
  s         = mystrtok(buff, '\t'); free(s);
  i         = 0;
  while ((s = mystrtok(0, '\t'))) {
    (*colnames)[i] = strdup(s); free(s);
    i++;
  }

  //
  //   allocate a big chunk of data
  //
  *rownames = (char**)malloc(mynmax * sizeof(char*));
  *data     = (int**)malloc(mynmax * sizeof(int*));

  myn = 0;
  while (1) {
    fgets(buff, mynmax, fp);

    if (feof(fp)) {
      break;
    }

    //printf("buff=%s\n", buff);
    
    //
    // allocate a line
    //
    (*data)[myn] = (int*)malloc(mym * sizeof(int));

    //
    // get the row name 
    //
    s = mystrtok(buff, '\t');
    (*rownames)[myn] = strdup(s);
    free(s);


    // get the values 
    i = 0;
    while ((s = mystrtok(0, '\t'))) {
      //printf("%d:%s\n", i, s);
      ;if (strcmp(s, "inf") == 0) {
	(*data)[myn][i] = INT_MAX;
      } else {
	(*data)[myn][i] = atoi(s);
      }

       free(s);
       i ++;

       
    }
    
    
    //s = mystrtok(buff, '\t');
    //(*colnames)[i] = strdup(s);
    //i++;
    //}
  
    //printf("myn=%d\n", myn);
  
    //*m = mym;
    myn ++;

    if (myn == mynmax) {
      die("readFloatTable: running out of memory, please recompile ..\n"); 
    }

    #ifdef GENEMAX
    if (myn == GENEMAX) {
      break;
    }
    #endif
  }

  //printf("ehat ?\n");

  *n = myn;
  *m = mym;
  fclose(fp);
  free(buff);

}



//
//     CHAR* 
//
void readStringTable(char* filename, char**** data, int* n, int* m) 
{
  
  char* s;
  int   mym;
  int   myn;
  char* buff;
  FILE* fp;
  int   i;
  int   mynmax = 100000;

  buff  = (char*)malloc(100000 * sizeof(char));
  fp    = fopen(filename, "r");
  if (!fp) {
    die("readFloatTable: please enter a valid filename\n");
  }

  //
  // get the first line and estimate the number of columns
  //

  fgets(buff, mynmax, fp);
  chomp(buff);

  s = mystrtok(buff, '\t'); free(s);
  mym = 1;
  while ((s = mystrtok(0, '\t'))) {
    mym ++; free(s);
  }

  //
  // allocate the right number of columns
  //
  *data      = (char***)malloc(mynmax * sizeof(char**));
  (*data)[0] = (char**)malloc(mym * sizeof(char*));
  chomp(buff);
  s   = mystrtok(buff, '\t');
  (*data)[0][0] = strdup(s); free(s);
  i   = 1;
  while ((s = mystrtok(0, '\t'))) {
    (*data)[0][i] = strdup(s); free(s);
    i++;
  }


  myn = 1;
  while (1) {

    fgets(buff, mynmax, fp);

    
    
    if (feof(fp)) {
      break;
    }
    
    chomp(buff);


    (*data)[myn] = (char**)malloc(mym * sizeof(char*));

    s = mystrtok(buff, '\t');
    (*data)[myn][0] = strdup(s);  free(s);
    i = 1;
    while ((s = mystrtok(0, '\t'))) {

      (*data)[myn][i] = strdup(s);
      free(s);
      i ++;

    }
    

	  
    myn ++;

    if (myn == mynmax) {
      die("readFloatTable: running out of memory, please recompile ..\n"); 
    }
    
    #ifdef GENEMAX
    if (myn == GENEMAX) {
      break;
    }
    #endif
  }

  *n = myn;
  *m = mym;

  free(buff); fclose(fp);
}




//
//  read a array of floats
//
float* readFloatSet(char* setfile, int* n) 
{

  
  FILE*   fp;
  int     nb = 0;
  char*   buff;  
  float*  oset;
  float   f;

  // create a line buffer
  buff = (char*)malloc(10000 * sizeof(char));
   
  // open setfile
  fp = fopen(setfile, "r");
  if (!fp) {
    printf("could not open set data %s\n", setfile);
    exit(0);
  }

  // allocate a big chunk of memory
  oset = (float*)malloc(10000 * sizeof(float));
  if (!oset) {
    printf("could not allocate memory for oset\n");
    exit(0);
  }

  //
  // go thru setfile
  //
  nb = 0;
  while (!feof(fp)) {
    fscanf(fp, "%f\n", &f);
    oset[nb] = f;
    nb++;
  }

  *n = nb;
  
  fclose(fp);
  
  return oset;
}

char** readStringSet(char* setfile, int* n) 
{

   // ok now get the set of ORFS
   
  FILE*  fp;
  int    nb = 0;
  char*  buff;
  char** oset;
  
 
    
  nb   = nbLinesInFile(setfile);
  
  *n   = nb;
 

  buff = (char*)malloc(10000 * sizeof(char));
  fp = fopen(setfile, "r");
  if (!fp) {
    printf("could not open set data %s\n", setfile);
  }
 
  // alloc memory
  oset = (char**)malloc(nb * sizeof(char*));
  if (!oset) {
    printf("could not allocate memory for oset\n");
    exit(0);
  }

  nb = 0;
  while (!feof(fp)) {
    fscanf(fp, "%s\n", buff);
    oset[nb] = (char*)calloc(300, sizeof(char));
    strcat(oset[nb], buff);
    nb++;
  }
  
  fclose(fp);
  
  return oset;
}



//
//  new implementation of strtok
//
char* mystrtok(char* s, char delim) 
{
 
  int posc;
  int l;
  char* myretstr;
  
  static char* sp;
  static int   posi;
  
  

  if (s != NULL) {
    sp   = s;
    posi = 0;
  }

  l = strlen(sp);
  
  if (posi > l)
    return NULL;


  // move to the next delim
  posc = posi;
  while ((posc < l) && (sp[posc] != delim)) { 
    ///printf("posc++(becomes %d,c=%c)\n", posc+1, sp[posc+1]);
    posc++;
  }
  

  // alloc enough memory
  myretstr = (char*)calloc(posc - posi + 1, sizeof(char));

  // cpy 
  strncat(myretstr, sp + posi, posc - posi);
  myretstr[posc - posi] = '\0';

  // update posi for next call  (posc points to the current delimiter)
  posi = posc+1;
  
  
  //printf("posi updated to %d\n", posi);

  return myretstr;
  
}




// remove the last \n, wherever it is ..
void chomp(char* s) 
{
  
  int j;
  
  j = strlen(s);
  while (j >= 0) {
    if (s[j] == '\n') {
       s[j] = '\0';
       return;
    }
    j--;
  }
  
  
}




char* complement(char* s) 
{

  int l = strlen(s);
  char* c;
  char  d;
  int i;
  
  c = (char*)calloc(l+1, sizeof(char));

  for (i=l-1; i>=0; i--) {
    if (s[i] == 'A')
      d = 'T';
    else if (s[i] == 'T')
      d = 'A';
    else if (s[i] == 'G')
      d = 'C';
    else if (s[i] == 'C') 
      d = 'G';
    else if (s[i] == '[')
      d = ']';
    else if (s[i] == ']')
      d = '[';
    else if (s[i] == '.')
      d = '.';
    else
      d = 'N';

    strncat(c, &d, 1);
  }
	 
  return c;
  
}


char* mybasename(char* filename)
{
  int l = strlen(filename);
  int posi, posj;
  int i;

  //  get the position of the first . from the right
  posi = -1;
  for (i=l-1; i>=0; i--) {
    if (filename[i] == '.') {
      posi = i;
      break;
    }
  }

  //printf("posi = %d\n", posi);
    
  //  get the position of the first / from the right
  posj = -1;
  for (i=posi-1; i>=0; i--) {
    if (filename[i] == '/') {
      posj = i;
      break;
    }
  }
  
  //printf("posj = %d\n", posj);


  if (posj < -1)
    posj = 0;

  if (posi < -1)
    posi = 0;

  return substr(filename, posj+1, posi - posj - 1);
  
}







char* get_parameter(int argc, char** argv, char* param)
{
  int i = 0;
  while ((i < argc) && (strcmp(param, argv[i]))) 
    i++;
  if (i<argc)
    return (argv[i+1]);
  else
    return "";
}


int exist_parameter(int argc, char** argv, char* param) {
  if (strlen(get_parameter(argc, argv, param)) > 0) 
    return 1;
  else
    return 0;
}


void die(char* s) 
{
  printf("%s", s);
  exit(0);
}


//
//  returns the next sequence in the file fp
//
/* USAGE

 int   nextSequence_started1 = 0;
 int   nextSequence_ended1   = 0;
 char* nextSequence_currentLine1; 
 int   size1;
 char* name1;
 char* fastafile1;
 FILE* fp1;
 char* seq1;

 nextSequence_currentLine1 = (char*)malloc(50000 * sizeof(char));
 fp1 = fopen(fastafile1, "r");
 if (!fp1) {
  printf("cannot open %s ..\n", fastafile1);
  exit(0);
 }

 while ( (seq1 = nextSequence(fp1, &name1, &size1, &nextSequence_started1, &nextSequence_ended1, nextSequence_currentLine1)) ) {   
 

 }

*/


char* nextSequence(FILE* fp, char** name, int* size, int* nextSequence_started, int* nextSequence_ended, char* nextSequence_currentLine) 
{
  
  int len = 20000;
  char* seq;
  

  // not yet started, get the first line
  if (*nextSequence_started == 0) {
    fgets(nextSequence_currentLine, len, fp);  
    chomp(nextSequence_currentLine);      
    *nextSequence_started = 1;
 
    
 }
  
  // if file is finished, return 0
  if (*nextSequence_ended == 1) {
    return 0;
  }

  // if started, line should be filled with ">..."
  if (nextSequence_currentLine[0] == '>') {
    

    *name = strdup(nextSequence_currentLine + 1);
    
    // create a new line
    seq = (char*)calloc(100000, sizeof(char) ); 

    while (1) {
      
      fgets(nextSequence_currentLine, len, fp);    
      
      if (feof(fp)) {
	*nextSequence_ended  = 1;
	realloc( seq, strlen(seq) + 1);
	return seq;
      }

      chomp(nextSequence_currentLine);      
      if (strlen(nextSequence_currentLine) == 0) { 
	continue;
      }
      
      if (nextSequence_currentLine[0] == '>') {	
	realloc( seq, strlen(seq) + 1);
	return seq;	
      } else {
	strcat(seq, nextSequence_currentLine);
      }
      
    }

  } else return 0;
  
}


/*
 *  open a FASTA file
 *
 */
int seqI_open(seqI* s, char* file)
{
  
  s->fp = fopen(file, "r");
  if (s->fp == 0) {
    return 0;
  }  
  s->nextSequence_started = 0;
  s->nextSequence_ended   = 0;
  s->nextSequence_currentLine = (char*)malloc(50000 * sizeof(char));

  return 1;
}

void seqI_close(seqI* s) 
{
  fclose(s->fp);
}

/*
 *  each call to this function returns the 
 *             next sequence in the FASTA file
 *
 */
char* seqI_nextSequence(seqI* s, char** name, int* size)
{
  
  int len        = 20000;
  int seqlen     = 100000;
  int seqlen_inc = 100000;
  char *seq_new;
  char* seq;
  

  // not yet started, get the first line
  if (s->nextSequence_started == 0) {
    fgets(s->nextSequence_currentLine, len, s->fp);  
    chomp(s->nextSequence_currentLine);      
    s->nextSequence_started = 1;
  }
  
  // if file is finished, return 0
  if (s->nextSequence_ended == 1) {
    return 0;
  }
  
  // if started, line should be filled with ">..."
  if ((s->nextSequence_currentLine)[0] == '>') {
    
    // copy line, starting from pos 1
    *name = strdup(s->nextSequence_currentLine + 1);
    
    // create a new line
    seq = (char*)calloc(seqlen, sizeof(char) ); 

    while (1) {
      
      fgets(s->nextSequence_currentLine, len, s->fp);    
      
      if (feof(s->fp)) {
	s->nextSequence_ended  = 1;
	realloc( seq, strlen(seq) + 1);
	return seq;
      }

      chomp(s->nextSequence_currentLine);      
      if (strlen(s->nextSequence_currentLine) == 0) { 
	continue;
      }
      
      if ((s->nextSequence_currentLine)[0] == '>') {	
	realloc( seq, strlen(seq) + 1);
	return seq;	
      } else {
	
	// do we have enough memory reserved ?
	if (strlen(seq) + strlen(s->nextSequence_currentLine) > seqlen) {
	  seq_new = (char*)malloc( (seqlen + seqlen_inc) * sizeof(char));
	  memcpy(seq_new, seq, strlen(seq));
	  free(seq);
	  seq = seq_new;
	  seqlen += seqlen_inc;
	}	
	strcat(seq, s->nextSequence_currentLine);
	
      }
      
    }

  } else return 0;
  
}




char* uc(char* seq) 
{
  int   l;
  int   i;
  char* sequp;
  
  l     = strlen(seq);
  sequp = strdup(seq);
  for (i=0; i<l; i++) {
    sequp[i] = toupper(seq[i]);
  }
  
  return sequp;
}

char*** initialize_aacode() 
{
  int   i = 0;
  char*** myaacode;
  myaacode = (char***)malloc(64*sizeof(char**));
  for (i=0; i<64; i++) {
    myaacode[i] = (char**)malloc(2*sizeof(char*));
  }
  

  myaacode[0][0]="TCA"; myaacode[0][1]="S";
  myaacode[1][0]="TCG"; myaacode[1][1]="S";
  myaacode[2][0]="TCC"; myaacode[2][1]="S";
  myaacode[3][0]="TCT"; myaacode[3][1]="S";
  myaacode[4][0]="TTT"; myaacode[4][1]="F";
  myaacode[5][0]="TTC"; myaacode[5][1]="F";
  myaacode[6][0]="TTA"; myaacode[6][1]="L";
  myaacode[7][0]="TTG"; myaacode[7][1]="L";
  myaacode[8][0]="TAT"; myaacode[8][1]="Y";
  myaacode[9][0]="TAC"; myaacode[9][1]="Y";
  myaacode[10][0]="TAA"; myaacode[10][1]="*";
  myaacode[11][0]="TAG"; myaacode[11][1]="*";
  myaacode[12][0]="TGT"; myaacode[12][1]="C";
  myaacode[13][0]="TGC"; myaacode[13][1]="C";
  myaacode[14][0]="TGA"; myaacode[14][1]="*";
  myaacode[15][0]="TGG"; myaacode[15][1]="W";
  myaacode[16][0]="CTA"; myaacode[16][1]="L";
  myaacode[17][0]="CTG"; myaacode[17][1]="L";
  myaacode[18][0]="CTC"; myaacode[18][1]="L";
  myaacode[19][0]="CTT"; myaacode[19][1]="L";
  myaacode[20][0]="CCA"; myaacode[20][1]="P";
  myaacode[21][0]="CCG"; myaacode[21][1]="P";
  myaacode[22][0]="CCC"; myaacode[22][1]="P";
  myaacode[23][0]="CCT"; myaacode[23][1]="P";
  myaacode[24][0]="CAT"; myaacode[24][1]="H";
  myaacode[25][0]="CAC"; myaacode[25][1]="H";
  myaacode[26][0]="CAA"; myaacode[26][1]="Q";
  myaacode[27][0]="CAG"; myaacode[27][1]="Q";
  myaacode[28][0]="CGA"; myaacode[28][1]="R";
  myaacode[29][0]="CGG"; myaacode[29][1]="R";
  myaacode[30][0]="CGC"; myaacode[30][1]="R";
  myaacode[31][0]="CGT"; myaacode[31][1]="R";
  myaacode[32][0]="ATT"; myaacode[32][1]="I";
  myaacode[33][0]="ATC"; myaacode[33][1]="I";
  myaacode[34][0]="ATA"; myaacode[34][1]="I";
  myaacode[35][0]="ATG"; myaacode[35][1]="M";
  myaacode[36][0]="ACA"; myaacode[36][1]="T";
  myaacode[37][0]="ACG"; myaacode[37][1]="T";
  myaacode[38][0]="ACC"; myaacode[38][1]="T";
  myaacode[39][0]="ACT"; myaacode[39][1]="T";
  myaacode[40][0]="AAT"; myaacode[40][1]="N";
  myaacode[41][0]="AAC"; myaacode[41][1]="N";
  myaacode[42][0]="AAA"; myaacode[42][1]="K";
  myaacode[43][0]="AAG"; myaacode[43][1]="K";
  myaacode[44][0]="AGT"; myaacode[44][1]="S";
  myaacode[45][0]="AGC"; myaacode[45][1]="S";
  myaacode[46][0]="AGA"; myaacode[46][1]="R";
  myaacode[47][0]="AGG"; myaacode[47][1]="R";
  myaacode[48][0]="GTA"; myaacode[48][1]="V";
  myaacode[49][0]="GTG"; myaacode[49][1]="V";
  myaacode[50][0]="GTC"; myaacode[50][1]="V";
  myaacode[51][0]="GTT"; myaacode[51][1]="V";
  myaacode[52][0]="GCA"; myaacode[52][1]="A";
  myaacode[53][0]="GCG"; myaacode[53][1]="A";
  myaacode[54][0]="GCC"; myaacode[54][1]="A";
  myaacode[55][0]="GCT"; myaacode[55][1]="A";
  myaacode[56][0]="GAT"; myaacode[56][1]="D";
  myaacode[57][0]="GAC"; myaacode[57][1]="D";
  myaacode[58][0]="GAA"; myaacode[58][1]="E";
  myaacode[59][0]="GAG"; myaacode[59][1]="E";
  myaacode[60][0]="GGA"; myaacode[60][1]="G";
  myaacode[61][0]="GGG"; myaacode[61][1]="G";
  myaacode[62][0]="GGC"; myaacode[62][1]="G";
  myaacode[63][0]="GGT"; myaacode[63][1]="G";

  //printf("%s\n", myaacode[49][0]);
  
  return myaacode;
}

int indexOfKmer(char* kmer) 
{
  
  int idx = 0;
  int i;
  int l   = strlen(kmer);
  int *N;

  N = (int*)calloc(256, sizeof(int));
  
  N['A'] = 0;
  N['C'] = 1;
  N['T'] = 2;
  N['G'] = 3;


  for (i=0; i<l; i++) {
    idx += N[ (int)(kmer[i]) ] * powint(4,(l - i - 1)); 
  }

  free(N);

  return idx;

}

int powint(int x, int y) 
{
  
  int o = 1;
  int i = 0;

  for (i=0; i<y; i++) {
    o *= x; 
  }
  
  return o;
}


char codon2aa(char* codon) 
{

  if (strcmp(codon, "TCA") == 0) return 'S';
  if (strcmp(codon, "TCG") == 0) return 'S';
  if (strcmp(codon, "TCC") == 0) return 'S';
  if (strcmp(codon, "TCT") == 0) return 'S';
  if (strcmp(codon, "TTT") == 0) return 'F';
  if (strcmp(codon, "TTC") == 0) return 'F';
  if (strcmp(codon, "TTA") == 0) return 'L';
  if (strcmp(codon, "TTG") == 0) return 'L';
  if (strcmp(codon, "TAT") == 0) return 'Y';
  if (strcmp(codon, "TAC") == 0) return 'Y';
  if (strcmp(codon, "TAA") == 0) return '*';
  if (strcmp(codon, "TAG") == 0) return '*';
  if (strcmp(codon, "TGT") == 0) return 'C';
  if (strcmp(codon, "TGC") == 0) return 'C';
  if (strcmp(codon, "TGA") == 0) return '*';
  if (strcmp(codon, "TGG") == 0) return 'W';
  if (strcmp(codon, "CTA") == 0) return 'L';
  if (strcmp(codon, "CTG") == 0) return 'L';
  if (strcmp(codon, "CTC") == 0) return 'L';
  if (strcmp(codon, "CTT") == 0) return 'L';
  if (strcmp(codon, "CCA") == 0) return 'P';
  if (strcmp(codon, "CCG") == 0) return 'P';
  if (strcmp(codon, "CCC") == 0) return 'P';
  if (strcmp(codon, "CCT") == 0) return 'P';
  if (strcmp(codon, "CAT") == 0) return 'H';
  if (strcmp(codon, "CAC") == 0) return 'H';
  if (strcmp(codon, "CAA") == 0) return 'Q';
  if (strcmp(codon, "CAG") == 0) return 'Q';
  if (strcmp(codon, "CGA") == 0) return 'R';
  if (strcmp(codon, "CGG") == 0) return 'R';
  if (strcmp(codon, "CGC") == 0) return 'R';
  if (strcmp(codon, "CGT") == 0) return 'R';
  if (strcmp(codon, "ATT") == 0) return 'I';
  if (strcmp(codon, "ATC") == 0) return 'I';
  if (strcmp(codon, "ATA") == 0) return 'I';
  if (strcmp(codon, "ATG") == 0) return 'M';
  if (strcmp(codon, "ACA") == 0) return 'T';
  if (strcmp(codon, "ACG") == 0) return 'T';
  if (strcmp(codon, "ACC") == 0) return 'T';
  if (strcmp(codon, "ACT") == 0) return 'T';
  if (strcmp(codon, "AAT") == 0) return 'N';
  if (strcmp(codon, "AAC") == 0) return 'N';
  if (strcmp(codon, "AAA") == 0) return 'K';
  if (strcmp(codon, "AAG") == 0) return 'K';
  if (strcmp(codon, "AGT") == 0) return 'S';
  if (strcmp(codon, "AGC") == 0) return 'S';
  if (strcmp(codon, "AGA") == 0) return 'R';
  if (strcmp(codon, "AGG") == 0) return 'R';
  if (strcmp(codon, "GTA") == 0) return 'V';
  if (strcmp(codon, "GTG") == 0) return 'V';
  if (strcmp(codon, "GTC") == 0) return 'V';
  if (strcmp(codon, "GTT") == 0) return 'V';
  if (strcmp(codon, "GCA") == 0) return 'A';
  if (strcmp(codon, "GCG") == 0) return 'A';
  if (strcmp(codon, "GCC") == 0) return 'A';
  if (strcmp(codon, "GCT") == 0) return 'A';
  if (strcmp(codon, "GAT") == 0) return 'D';
  if (strcmp(codon, "GAC") == 0) return 'D';
  if (strcmp(codon, "GAA") == 0) return 'E';
  if (strcmp(codon, "GAG") == 0) return 'E';
  if (strcmp(codon, "GGA") == 0) return 'G';
  if (strcmp(codon, "GGG") == 0) return 'G';
  if (strcmp(codon, "GGC") == 0) return 'G';
  if (strcmp(codon, "GGT") == 0) return 'G';
  
  return -1;

}

char* translate(char* seq, char*** aacode) 
{
  
  char* seqaa;
  int   l, i;
  int*  codon2idx;
  char* codon;
  int   rabe;

  codon2idx = (int*)malloc(64 * sizeof(int));
  for (i=0; i<64; i++) {
    codon2idx[ indexOfKmer(aacode[i][0]) ] = i;
  }
  
  
  l = strlen(seq);

  seqaa = (char*)calloc((1 + l / 3), sizeof(char));

  rabe = l % 3;

  for (i=0; i<l-rabe-2; i+=3) {
    codon = substr(seq, i, 3);
    strcat(seqaa, aacode[ codon2idx[ indexOfKmer(codon) ] ][1]);
    free(codon);
  }
  
    
  free(codon2idx);
  
  return seqaa;

}


//
//  substring
//
char *substr(char *string, int start, int length) 
{
  char *substring;
  
  substring = (char*)calloc(length+1, sizeof(char));

  memcpy(substring, string+start, length*sizeof(char));
  substring[length] = '\0';

  return substring;
}



void showIntVector(int* v, int m) 
{
  int i;
  
  for (i=0; i<m; i++) {
    printf("%d\t", v[i]);
  }
  printf("\n");
}


void showFloatVector(float* v, int m) 
{
  int i;
  
  for (i=0; i<m; i++) {
    printf("%3.2f\t", v[i]);
  }
  printf("\n");
}



short* ftos_vector(float* v, int n) 
{
  int    i;
  short* f = (short*)malloc(n * sizeof(short));

  for (i=0; i<n; i++) {
    f[ i ] = (short)(v[i]);
  }

  return f;
  
}


float* itof_vector(int* v, int n) 
{
  int    i;
  float* f = (float*)malloc(n * sizeof(float));

  for (i=0; i<n; i++) {
    f[ i ] = (float)v[i];
  }

  return f;
  
}


//
//  get the column for a matrix
//
float* f_matrix_column(float** data, int j, int n) 
{
  float* v;
  int    i;

  v = (float*) malloc ( n * sizeof( float ));
  for (i=0; i<n; i++) {
    v[i] = data[i][j];
  }
  
  return v;
}



//
//  get the column for a matrix
//
short* s_matrix_column(short** data, int j, int n) 
{
  short* v;
  int    i;

  v = (short*) malloc ( n * sizeof( short ));
  for (i=0; i<n; i++) {
    v[i] = data[i][j];
  }
  
  return v;
}



//
//  get the column for a matrix
//
float* stof_matrix_column(short** data, int j, int n) 
{
  float* v;
  int    i;

  v = (float*) malloc ( n * sizeof( float ));
  for (i=0; i<n; i++) {
    v[i] = (float)(data[i][j]);
  }
  
  return v;
}



//
//  get the column for a matrix
//
float* itof_matrix_column(int** data, int j, int n) 
{
  float* v;
  int    i;

  v = (float*) malloc ( n * sizeof( float ));
  for (i=0; i<n; i++) {
    v[i] = (float)(data[i][j]);
  }
  
  return v;
}


//
//  get the column for a matrix
//
int* stoi_matrix_column(short** data, int j, int n) 
{
  int* v;
  int    i;

  v = (int*) malloc ( n * sizeof( int ));
  for (i=0; i<n; i++) {
    v[i] = (int)(data[i][j]);
  }
  
  return v;
}


int* i_matrix_column(int** data, int j, int n) 
{
  int* v;
  int    i;

  v = (int*) malloc ( n * sizeof( int ));
  for (i=0; i<n; i++) {
    v[i] = (int)(data[i][j]);
  }
  
  return v;
}



int* c_matrix_column(char** data, int j, int n) 
{
  int* v;
  int    i;

  v = (int*) malloc ( n * sizeof( int ));
  for (i=0; i<n; i++) {
    v[i] = (int)(data[i][j]);
  }
  
  return v;
}
