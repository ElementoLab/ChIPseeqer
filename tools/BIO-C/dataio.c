#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <sys/stat.h>

#include "dataio.h"
#include "statistics.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif


int verbose;

#ifndef NAN
static const double NAN = (0.0 / 0.0);
#endif

// read entire file into string
char * readFile(char * fname){
	FILE *fp;
	long len;
	char *buf;
	fp=fopen(fname,"rb");
	if (fp == NULL)
			return (NULL);
	fseek(fp,0,SEEK_END); //go to end
	len=ftell(fp); //get position at end (length)
	rewind(fp); //go to beg.

	buf=(char *)malloc(len+1); //malloc buffer
	fread(buf, 1, len, fp); //read into buffer
	buf[len]= '\0';
	fclose(fp);
	return (buf);
}



void LineIclose(LineI* li)
{
  fclose(li->fp);
}

int LineIopen(LineI* li, char* file)
{
  int cnt        = 0;

  li->bufflen = 1000000;

  li->fp         = fopen(file, "r");
  li->eof       = 0;

  if (li->fp == 0) {
    printf("LineIopen: cannot open file %s\n", file);
    exit(0);
  }

  // alloc 
  li->buff   = (char*)calloc(li->bufflen, sizeof(char)); 
  if (!(li->buff)) {
    die("Cannot allocate li->buff\n");
  }

  // read first chunk 
  cnt = fread(li->buff, 1, li->bufflen, li->fp); 
  //printf("Read %d chars from %s.\n", cnt, file);
  if (cnt == 0) {
    free(li->buff);
    li->eof = 1;
    return 0;
  }

  // store chunck address
  li->buffptr = li->buff; // this is for free()
  
  return 1;
}


char* nextLine(LineI* li)
{
  int cnt = 0;
  char* fullline = 0;
  char* curLine = 0;
  char* firstline = 0;

  if (li->eof == 1)
    return 0;

  // get the first line
  curLine = strdup(strsep(&(li->buff), "\n"));

  //printf("buff=%d, curLine=%s\n", li->buff, curLine);
 

  // if end of buff was reached, need to read one more chunk, so that we can add end of line
  if (li->buff == 0) {
    
    // free current 
    //free(li->buffptr);
    
    // read another chunk
    li->buff = li->buffptr;
    cnt = fread(li->buff, 1, li->bufflen, li->fp); 
    //printf("read %d chars\n", cnt);
    //getchar();
    if (cnt != 0) { // not end of file
      
      // get the first line
      firstline = strsep(&(li->buff), "\n");
      
      // alloc new fullline
      //printf("Alloc full line of size %d+%d\n", strlen(curLine), strlen(firstline)); 
      //getchar();
      fullline = (char*)calloc(strlen(curLine)+strlen(firstline)+1, sizeof(char));
      
      // concatenate what is left of prev line
      strcat(fullline, curLine);
      free(curLine);
      //printf("fulline=%s\n", fullline);
      // concatenate first line
      strcat(fullline, firstline);
      //printf("fulline=%s\n", fullline);

    } else {

      // set eof
      li->eof = 1;

      // free ptr
      free(li->buffptr);

      return 0;
      //getchar();
    }

  }
  
  if (fullline != 0) 
    return fullline;
  else 
    return curLine;
  
 
}




// A = 00
// C = 01
// G = 10
// T = 11



void readmiRNASeeds(char* seedfile, char*** seeds, char*** mirnames, int* numseeds)
{
  FILE* f;
  int   mynmax = 10000;
  char** a;
  int    m;
  char* buff   = (char*)calloc(mynmax, sizeof(char));
  
  *numseeds    = nbLinesInFile(seedfile);  
  *seeds       = (char**)calloc((*numseeds+1), sizeof(char*));
  *mirnames    = (char**)calloc((*numseeds+1), sizeof(char*));

   // read snvfile
  int idx = 0;
  f = fopen(seedfile, "r"); 
  if (f == 0) 
    die("Cannot open seed file.\n");
  buff = (char*)calloc(mynmax, sizeof(char));
  while (fgets(buff, mynmax, f) != NULL) {
    chomp(buff);
    split_line_delim_sep(buff, "\t", &a, &m);
    (*seeds)[idx]    = strdup(a[0]);
    (*mirnames)[idx] = strdup(a[2]);
    free(a);
    idx++;
  }

}



int file_exists (char * fileName)
{
   struct stat buf;
   int i = stat ( fileName, &buf );
   /* File found */
   if ( i == 0 )
     {
       return 1;
     }
   return 0;
   
}


unsigned char* encodeDNAseq(char* s) 
{
  
  unsigned char* p;
  int i;
  int l = strlen(s);
  int o = l / 4 + 1;
  p = (unsigned char*)calloc(o, sizeof(unsigned char));
  if (p == 0) {
    die("Cannot allocate binarized DNA seq.\n");
  }
  
  for (i=0; i<l; i++) {
    int  o = i / 4;  // determine which char to modify, knowing that 1 char = 
    int  n = 2 * (i % 4);  // which offset to modify
    // 1

    //printf("o=%d, n=%d\n", o, n);
    
    if (s[i] == 'A') {
      // don't do anything (00)
    } else if (s[i] == 'C') {
      p[o] = p[o] | (1 << (n+1));  // (01)
      
      //if ((( p[o] & (1 << (n)) ) == 0) && ((p[o] & (1 << (n+1))) == 1)) 
      //printf("we should have a C here. p1=%d p2=%d\n", (int)( p[o] & (1 << (n)) ), (int)( p[o] & (1 << (n+1)) ) );
      //else
      //printf("we don't have a C\n");

    } else if (s[i] == 'G') {
      p[o] = p[o] | (1 << n);    // (10)
      //printf("we should have a C here. p1=%d p2=%d\n", (int)( p[o] & (1 << (n)) ), (int)( p[o] & (1 << (n+1)) ) );
    } else if (s[i] == 'T') {
      p[o] = p[o] | (1 << n);    // (11)
      p[o] = p[o] | (1 << (n+1));  // 
      //printf("we should have a C here. p1=%d p2=%d\n", (int)( p[o] & (1 << (n)) ), (int)( p[o] & (1 << (n+1)) ) );
    }
    //printf("p[%d] = %d\n", o, (int)(p[o]));
  }
  
  return p;
}

  
char* decodeDNAseq(unsigned char* p, int l) 
{
  
  char* c;
  int i;
  c = (char*)calloc(l+1, sizeof(char));
  if (c == 0) {
    die("Cannot allocate binarized DNA seq.\n");
  }
  
  for (i=0; i<l; i++) {
    int  o = i / 4;        // determine which char to modify, knowing that 1 char = 
    int  n = 2 * (i % 4);  // which offset to modify    
    if ((      (int)( p[o] & (1 << (n)) ) == 0) && ( (int) (p[o] & (1 << (n+1))) == 0))
      c[i] = 'A';
    else if (( (int)( p[o] & (1 << (n)) ) == 0) && ( (int) (p[o] & (1 << (n+1))) >  0))
      c[i] = 'C';
    else if (( (int)( p[o] & (1 << (n)) ) > 0 ) && ( (int) (p[o] & (1 << (n+1))) == 0))
      c[i] = 'G';
    else if (( (int)( p[o] & (1 << (n)) ) > 0 ) && ( (int) (p[o] & (1 << (n+1))) >  0))
      c[i] = 'T';
  }

  c[l] = '\0';
 
  return c;
}
    


// DEBUG:
// 	safer memory management 
void create_safe_binarized_array(int n, unsigned char ** p) {

  // unsigned char* p;
  int o = n / 8 + 1;
  
  *p = (unsigned char**)calloc(o, sizeof(unsigned char*));
  if (*p == 0) {
    die("Cannot allocate binarized array.\n");
  }
  
  // return p;
}


unsigned char* create_binarized_array(int n) {

  unsigned char* p;
  int o = n / 8 + 1;
  
  p = (unsigned char*)calloc(o, sizeof(unsigned char));
  if (p == 0) {
    die("Cannot allocate binarized array.\n");
  }
  
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


void split_line(char* l, char*** p) {
  
  int mym = 0;
  (*p) = (char**)malloc(200000*sizeof(char*));
  
  char* s = mystrtok(l, '\t'); 
  (*p)[mym++] = s;  
  while ((s = mystrtok(0, '\t'))) {
    (*p)[mym++] = s;
  }

  
}


void split_line_delim(char* l, char* delim, char*** p, int* m) {
  
  int mym = 0;
  (*p) = (char**)malloc(200000*sizeof(char*));
  if (*p == 0) {
    die("Cannot alloc mem to store line elements in split_line_delim\n");
  }
  char* s = strtok(l, delim); 
  (*p)[mym++] = s;  
  while ((s = strtok(0, delim))) {
    (*p)[mym++] = s;
  }

  *m = mym;
  
}


void split_line_delim_sep(char* l, char* delim, char*** p, int* m) 
{
  
  int mym = 0;
  (*p) = (char**)malloc(200000*sizeof(char*));

  char* s = strsep(&l, delim);

  (*p)[mym++] = s;  
  while ((s = strsep(&l, delim))) {
    (*p)[mym++] = s;
  }

  *m = mym;
  
}



int getNumPositiveEntries(int* v, int n)
{
  
  int i;
  int cnt = 0;

  for (i=0; i<n; i++) {
    if (v[i] > 0) 
      cnt ++;
      
  }
  
  return cnt;
  
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


int nbSequencesInFastaFile(char* file) 
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
    if (buff[0] == '>') {
      n++;
    }
    fgets(buff, 100000, fp);
    //n++;
  }  
  fclose(fp);

  free(buff);
  return n;
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
  //fgets(buff, 100000, fp);
  n = 0;
  while (fgets(buff, 100000, fp) != 0) {
    n++;
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
    die("readAASequences: please enter a valid filename\n");
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


//
//  read in a collection of sites from AlignACE
//

void readAlignACEMotif(char* filename, int* m, int* n, int*** data) 

{
  
  //char* s;
  int   mym;
  int   myn;
  char* buff;
  FILE* fp;
  int   i;
  int   mynmax = 500;
  int   l = 0;

  int   cnt_stars = 0;

  buff  = (char*)malloc(100000 * sizeof(char));

  fp    = fopen(filename, "r");
  if (!fp) {
    die("readAlignACEmotif: please enter a valid filename\n");
  }


  mym = 0;

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

    l = strlen(buff);

    if (buff[0] == 'M')
      continue;

    // if cnt
    cnt_stars = 0;
    for (i=0; i<l; i++) {
      if (buff[i] == '*') {
	cnt_stars++;
	break;
      }
    }

    if (cnt_stars > 0) 
      continue;
    
    (*data)[myn] = (int*)malloc(l * sizeof(int));

    for (i=0; i<l; i++) {
      (*data)[myn][i] = data_nucl( buff[i] );
    }
    
    myn ++;
  
   
  }


  *n = myn;
  *m = l;

  
  fclose(fp);
}


int data_nucl(char c)
{

  if (c == 'A')
    return 0;
  else if (c == 'C')
    return 1;
  else if (c == 'G')
    return 2;
  else if (c == 'T')
    return 3;
  //  else if (c == 'N')
  //  return 4;
  else { 
    //fprintf(stderr, "Warning: %c\n", c); fflush(stdout);
    return -1;
  }


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
// read in a FIRE profile (integers)
//
void readIntProfile(char* filename, int* n, int** data, char*** rownames) 
{
  
  char* s;
  //int   mym;
  int   myn;
  char* buff;
  FILE* fp;
  //int   i;
  int   mynmax = 100000;

  buff  = (char*)malloc(100000 * sizeof(char));
  fp    = fopen(filename, "r");
  if (!fp) {
    die("readIntProfile: please enter a valid filename\n");
  }

  //
  // get the first line and estimate the number of columns
  //  
  fgets(buff, mynmax, fp);
  
  
  //
  //   allocate a big chunk of data
  //
  *rownames = (char**)malloc(mynmax * sizeof(char*));
  *data     = (int*)malloc(mynmax * sizeof(int));

  myn = 0;
  while (1) {

    fgets(buff, mynmax, fp);

    if (feof(fp)) {
      break;
    }

    //
    // get the row name 
    //
    s = mystrtok(buff, '\t');
    (*rownames)[myn] = strdup(s);
    free(s);

    // get the value 
    s = mystrtok(0, '\t');
    if (strcmp(s, "inf") == 0) {
      (*data)[myn] = INT_MAX;
    } else {
      (*data)[myn] = atoi(s);
    }    
    free(s);

    myn ++;

    if (myn == mynmax) {
      die("readIntProfile: running out of memory, please recompile ..\n"); 
    }

    #ifdef GENEMAX
    if (myn == GENEMAX) {
      break;
    }
    #endif
  }

  *n = myn;

  fclose(fp);
  free(buff);

}




//
// read in a FIRE profile (float)
//
void readFloatProfile(char* filename, int* n, float** data, char*** rownames) 
{
  
  char* s;
  //int   mym;
  int   myn;
  char* buff;
  FILE* fp;
  //int   i;
  int   mynmax = 100000;

  buff  = (char*)malloc(100000 * sizeof(char));
  fp    = fopen(filename, "r");
  if (!fp) {
    die("readFloatProfile: please enter a valid filename\n");
  }

  //
  // get the first line and estimate the number of columns
  //  
  fgets(buff, mynmax, fp);
  
  
  //
  //   allocate a big chunk of data
  //
  *rownames = (char**)malloc(mynmax * sizeof(char*));
  *data     = (float*)malloc(mynmax * sizeof(float));

  myn = 0;
  while (1) {

    fgets(buff, mynmax, fp);

    if (feof(fp)) {
      break;
    }

    //
    // get the row name 
    //
    s = mystrtok(buff, '\t');
    (*rownames)[myn] = strdup(s);
    free(s);

    // get the value 
    s = mystrtok(0, '\t');
    if (strcmp(s, "inf") == 0) {
      (*data)[myn] = (float)INT_MAX;
    } else {
      (*data)[myn] = atof(s);
    }    
    free(s);

    myn ++;

    if (myn == mynmax) {
      die("readFloatProfile: running out of memory, please recompile ..\n"); 
    }

    #ifdef GENEMAX
    if (myn == GENEMAX) {
      break;
    }
    #endif
  }

  *n = myn;

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
  int   nb = -1;
  //char* buffc = 0;

  buff  = (char*)malloc(100000 * sizeof(char));
  fp    = fopen(filename, "r");
  if (!fp) {
    die("readFloatTable: please enter a valid filename\n");
  }

  // get the first line and estimate the number of columns
  nb   = nbLinesInFile(filename);

  // get line 1
  fgets(buff, mynmax, fp);
  chomp(buff);
  //buffc = strdup(buff); // copy line 1 
  s = mystrtok(buff, '\t'); free(s);
  mym = 1;
  while ((s = mystrtok(0, '\t'))) {
    mym ++; free(s);
  }

  //
  // allocate the right number of columns
  //
  *data      = (char***)malloc(nb * sizeof(char**));
  rewind(fp);
  

  /*
  (*data)[0] = (char**)malloc(mym * sizeof(char*));
  chomp(buff);
  s   = mystrtok(buff, '\t');
  (*data)[0][0] = strdup(s); free(s);
  i   = 1;
  while ((s = mystrtok(0, '\t'))) {
    (*data)[0][i] = strdup(s); free(s);
    i++;
  }
  */

  myn = 0;
  while (fgets(buff, mynmax, fp) != 0) {
    
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


char char_complement(char s) 
{

  if (s == 'A')
    return 'T';
  else if (s == 'T')
    return 'A';
  else if (s == 'G')
    return 'C';
  else if (s == 'C') 
    return 'G';
  else 
    return 'N';
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



int exist_flag(int argc, char** argv, char* param)
{
  int i = 0;
  while ((i < argc) && (strcmp(param, argv[i])))
    i++;
  if (i<argc)
    return 1;
  else
    return 0;
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
  exit(1);
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

/* USAGE:

  char* name;
  int   size;
  seqI  si;
  char* myseq;

  // working for the human genome (fast)
  seqI_set_max_seqlen(&si, 3500000000UL);
  seqI_set_seqlen_inc(&si, 300000000UL);
  seqI_set_fread_chunk(&si, 100000000UL);

  seqI_open_fast(&si, argv[1]);
  
  while ( myseq = seqI_nextSequence_fast(&si, &name, &size) ) {
   
    printf("%d\n", strlen(myseq));
    free(myseq);

  }

*/

int seqI_open_fast_old(seqI* s, char* file)
{
  
  int   cnt       = 0;
  long int   total_cnt = 0;
  char* buff;
  char* seq;
  int   bufflen   = s->fread_chunksize;

  //printf("opening %s\n", file); fflush(stdout);

  s->fp = fopen(file, "r");
  if (s->fp == 0) {
    return 0;
  }  

  seq    = (char*)calloc(s->max_seqlen, sizeof(char));
  if (seq == 0) {
    printf("cannot allocate seq with n=%ld\n", s->max_seqlen);
    exit(0);
  }
  buff   = (char*)calloc(bufflen, sizeof(char)); 

  while (!feof(s->fp)) {
    cnt = fread(buff, 1, bufflen, s->fp);    
    printf("read %d char.\n", cnt);
    total_cnt += cnt;
    strcat(seq, buff);    
  }
  
  seq[total_cnt] = '\0';

  printf("len = %ld\n", strlen(seq));
    
  s->seq = seq;
  s->pos = 0;
  
  s->verbose    = 0;
  return 1;
}





void seqI_set_max_seqlen(seqI* s, unsigned long int n)
{
  
  //printf("setting max_seqlen to %ld\n", n);
  
  s->max_seqlen = n;

}

void seqI_set_fread_chunk(seqI* s, unsigned long int n)
{
  
  //printf("setting max_seqlen to %ld\n", n);
  
  s->fread_chunksize = n;

}

// remember to free s->chunk
int seqI_open_fast(seqI* s, char* file)
{
  
  int   cnt       = 0;
  //long int   total_cnt = 0;
  char* buff;
  //char* seq;
  int   bufflen   = s->fread_chunksize;

  //printf("opening %s\n", file); fflush(stdout);

  s->fp = fopen(file, "r");
  if (s->fp == 0) {
    return 0;
  }  

  // allocate buffer
  buff   = (char*)calloc(bufflen, sizeof(char)); 
  if (!buff) {
    die("Cannot allocate buff in seqI_open_fast\n");
  }

  // read chunk 
  cnt = fread(buff, 1, bufflen, s->fp);    
  s->chunk        = buff;

  // set len
  s->chunk_len = cnt;

  // advance to first >
  s->pos_in_chunk = 0;
  while ((s->pos_in_chunk < s->chunk_len) && (s->chunk[ s->pos_in_chunk ] != '>'))
    s->pos_in_chunk ++;

  if (s->verbose == 1)
    printf("read %d char.\n", cnt);

  s->enter_new_seq = 1;
  s->eof           = 0;
  s->freadeof      = 0;
  return 1;
}



// returns a sequence read from a fasta file
char* seqI_nextSequence_fast(seqI* s, char** name, int* size)
{


  int cnt    = 0;
  char* buff = 0;
  //char* seq  = 0;
  int   bufflen   = s->fread_chunksize;
  char  c;

  if (s->eof == 1)
    return 0;

  while (1) {  // loop over the chunks

    //
    // concatenate until end of chunk or other >
    //
    while (s->pos_in_chunk < s->chunk_len) {

      //if (s->pos_in_chunk % 1000 == 0) 
      //printf("Pos in chunk = %ld\n", s->pos_in_chunk);
      
      // we get a character      
      c = s->chunk[ s->pos_in_chunk ];

      // ignore this stupid chracter
      if (c == '\r') {
	s->pos_in_chunk ++;
	continue;
      }

      // what to do with it (depends on the context)
      
      if (s->chunk[ s->pos_in_chunk ] == '>') {
	
	if (s->verbose == 1)
	  printf("Found > ! at pos = %ld\n", s->pos_in_chunk);

	if (s->enter_new_seq == 0) {
	  
	  // for next time
	  s->enter_new_seq = 1;

	  // deliberately do not increase counter  s->pos_in_chunk

	  // return name and seq
	  if (s->verbose == 1) {
	    printf("Returning seq\n"); fflush(stdout);
	    getchar();
	  }

	  *name = s->name;
	  return s->seq;

	} else if (s->enter_new_seq == 1) {
	
	  // allocate a new sequence
	  if (s->verbose == 1) {
	    printf("creating new seq\n"); getchar();
	  }
	  s->seq = (char*)calloc(s->seqlen_inc, sizeof(char));
	  if (s->seq == 0) {
	    printf("Problem allocating seq, n=%d\n", s->seqlen_inc);
	    exit(0);
	  }
	  
	  // create name
	  s->name = (char*)calloc(100, sizeof(char));
	  
	  // init positions
	  s->pos_in_name = 0;
	  s->pos_in_seq  = 0;
	  
	  // we are now in the name, not in the seq
	  s->inname = 1;
	  s->inseq  = 0;

	  // move
	  s->pos_in_chunk ++;
	  
	  // not entering new seq anymore
	  s->enter_new_seq = 0;
	  
	} // end if enter new sequence


      } else if (s->inname == 1) {

	// if we are in the middle of the name
	
	
	// if not a return carriage, continue
	if (c != '\n') {
	  
	  // add to name
	  s->name[ s->pos_in_name ] = c;
	  s->pos_in_name ++;
	  
	} else {
	  
	  // fiish name
	  s->name[ s->pos_in_name ] = '\0';
	  s->pos_in_name ++;
	  
	  // get out
	  s->inname = 0;
	  
	  // enter sequence
	  s->inseq  = 1;
	  
	}

	// move
	s->pos_in_chunk ++;
		
	
      } else if (s->inseq == 1) {
	
	if ((c != ' ') && (c != '\n')) {
	  s->seq[ s->pos_in_seq ] = c;
	  s->pos_in_seq ++;
	}
	
	// move
	s->pos_in_chunk ++;
		
      }

    } // end while in chunk
   
    
    // we have reached the end of the chunk
    //  but not the end of the filq
    if (s->freadeof != 1) {

      // free up current chunk
      free(s->chunk);
      
      // allocate buffer for new chunk
      buff   = (char*)calloc(bufflen, sizeof(char)); 
      if (!buff) {
	die("Cannot allocate buff in seqI_open_fast\n");
      }

      // read chunk 
      cnt = fread(buff, 1, bufflen, s->fp); 

      if (cnt < bufflen) {
	s->freadeof = 1;
      }
   
      s->chunk        = buff;

      // set len
      s->chunk_len = cnt;
      s->pos_in_chunk = 0;
      
    } else {

      //printf("end of file\n");

      s->eof = 1;
      
      // end is reached, return the seq and name
      *name = s->name;
      return s->seq;

    }
    
  } // end while (1)

}



char* seqI_nextSequence_fast_old(seqI* s, char** name, int* size)
{

  unsigned long int l = strlen(s->seq);
  int i = 0, j = 0;
  int inseq = 0;
  char* seq = 0;


  j = 0;
  while (s->pos < l) {

    if (s->seq[ s->pos ] == '>') {
  
      
    
      // two cases, if not yet started, enter seq
      if (inseq == 0) {
	
	//printf("allocating %d\n", s->seqlen_inc);
	
	seq = (char*)calloc(s->seqlen_inc, sizeof(char));
	if (seq == 0) {
	  printf("Problem allocating seq, n=%d\n", s->seqlen_inc);
	  exit(0);
	}

	*name = (char*)calloc(100, sizeof(char));
	i     = 0;
	s->pos++;
	while (s->seq[ s->pos ] != '\n') {
	  (*name)[i] = s->seq[ s->pos ];
	  s->pos++;
	  i ++;
	}

	(*name)[i] = '\0';
	
	//printf("name = %s\n", *name);
	fflush(stdout);

	inseq = 1;

      } else {
	
	// means we are reaching the begining of a new seq
	// stay at the same position in file
	
	return seq;
      }

    } else {

      // add to seq
      if ((s->seq[ s->pos ] != '\n') && (s->seq[ s->pos ] != '\r')) {	
	seq[j] = s->seq[ s->pos ];
	j ++;
      }
      
      s->pos ++;
    }
    
    
  }

  return seq;
  
}






/*

Usage of new module:

  seqI     si;
  char*    seq;
  char*    name;
  int      size;

  if (seqI_open(&si, fastafile) == 0) {
    die("Error opening FASTA file\n");
  while (seq = seqI_nextSequence(&si, &name, &size)) {

  }
  seqI_close(&si);

*/


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
  s->seqlen_inc = 100000;
  s->verbose    = 0;
  return 1;
}

void seqI_close(seqI* s) 
{
  fclose(s->fp);
}

void seqI_set_seqlen_inc(seqI* s, int n)
{
  s->seqlen_inc = n;

}

void seqI_set_verbose(seqI* s, int v) {
  s->verbose = 1;
}


/*
 *  each call to this function returns the 
 *             next sequence in the FASTA file
 *
 */
char* seqI_nextSequence(seqI* s, char** name, int* size)
{
  
  int len        = 20000;
  int seqlen     = s->seqlen_inc;
  int seqlen_inc = s->seqlen_inc;


  char* seq;
  int   cnt = 0;

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
  
  //
  // if started, line should be filled with ">..."
  //
  if ((s->nextSequence_currentLine)[0] == '>') {
    
    // copy line, starting from pos 1
    *name = strdup(s->nextSequence_currentLine + 1);
    
    // create a new line
    seq = (char*)calloc(seqlen, sizeof(char) ); 
    if (! seq ) {
      die("dataio.c: Cannot allocate memory for seq in seqI_nextSequence\n");
    }
    while (1) {
      
      fgets(s->nextSequence_currentLine, len, s->fp);    

      if (s->verbose == 1) {
	cnt ++;
	printf("read %d lines.      \r", cnt);
      }
      
      if (feof(s->fp)) {
	s->nextSequence_ended  = 1;
	realloc( seq, strlen(seq) + 1);
	return seq;
      }

      chomp(s->nextSequence_currentLine);      

      //printf("%s\n", s->nextSequence_currentLine);

      if (strlen(s->nextSequence_currentLine) == 0) { 
	continue;
      }
      
      if ((s->nextSequence_currentLine)[0] == '>') {	
	//realloc( seq, strlen(seq) + 1);
	return seq;	
      } else {
	
	// do we have enough memory reserved ?
	if (strlen(seq) + strlen(s->nextSequence_currentLine) > seqlen) {
	  
	  /*
	  seq_new = (char*)malloc( (seqlen + seqlen_inc) * sizeof(char));
	  memcpy(seq_new, seq, strlen(seq));
	  free(seq);
	  seq = seq_new;
	  */

	  seq = (char*)realloc(seq, (seqlen+seqlen_inc)*sizeof(char));
	  if (seq == 0) {
	    printf("realloc failed\n");  fflush(stdout);
	    exit(0);
	  }
	  seqlen += seqlen_inc;
	}	

	strcat(seq, s->nextSequence_currentLine);
	
      }
      
    }

  } else return 0;
  
}


void uctransform(char* seq) 
{
  int   l;
  int   i;
  l     = strlen(seq);
  for (i=0; i<l; i++) {
    seq[i] = toupper(seq[i]);
  }

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



void lctransform(char* seq) 
{
  int   l;
  int   i;
  l     = strlen(seq);
  for (i=0; i<l; i++) {
    seq[i] = tolower(seq[i]);
  }

}


char* lc(char* seq) 
{
  int   l;
  int   i;
  char* sequp;
  
  l     = strlen(seq);
  sequp = strdup(seq);
  for (i=0; i<l; i++) {
    sequp[i] = tolower(seq[i]);
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

//585-2108

void showTwoIntVectors(int* v1, int* v2, int m) 
{
  int i;
  
  for (i=0; i<m; i++) {
    printf("%d\t%d\n", v1[i], v2[i]);
  }
  //printf("\n");
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


int* ctoi_vector(char* v, int n) 
{
  int    i;
  int* f = (int*)malloc(n * sizeof(int));

  for (i=0; i<n; i++) {
    f[ i ] = (int)(v[i]);
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
