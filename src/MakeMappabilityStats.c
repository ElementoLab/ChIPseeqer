#include <stdio.h>
#include <stdlib.h>
#include <search.h>
#include <string.h>
#include <math.h>


#include "dataio.h"
#include "statistics.h"
#include "information.h"
#include "hashtable.h"
#include "sequences.h"


#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)


int main(int argc, char** argv)
{

  // hash
  ENTRY   e;
  ENTRY*  ep;
  int     hashret;
  struct my_hsearch_data* hash_chr;

  // input
  FILE* fl;
  FILE* fn;
  FILE* fk;
  char* buff;
  int   mynmax = 10000;
  char** p;
  int    m;
  
  char**  chr;

  // chrlen
  int*  chrlen;

  // map
  unsigned char** map;
  
  // general
  int    i, j;
  
  int    k = 30;

  int*   nummap;

  if (argc == 1) {
    die("Args: chrlen chrN kmercount (the chr in chrlen must be in the same order as the chr in kmercount)\n");
  }

  // to index chromosomes
  hash_chr  = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(1000, hash_chr);
  if (hashret == 0) {
    printf("Could not create hash table ...\n");
    exit(0);
  } 

  
  buff  = (char*)malloc(mynmax * sizeof(char));
  
  // alloc map
  map    = (unsigned char**)calloc(50, sizeof(unsigned char*));
  nummap = (int*)calloc(50, sizeof(int));

  // load chr lengths
  chrlen = (int*)calloc(50, sizeof(int));
  chr    = (char**)calloc(50, sizeof(char*));

  fl = fopen(argv[1], "r");  

  int numchr = 0;
  while (!feof(fl)) {

    fgets(buff, mynmax, fl);
    if (feof(fl))
      break; 
    chomp(buff);
    split_line_delim(buff, "\t", &p, &m);

    e.key   = strdup(p[0]);
    e.data  = (char*)numchr;
    //printf("Entering %d\n", );
    hashret = my_hsearch_r(e, ENTER, &ep, hash_chr);
    if (hashret == 0) {
      printf("Could not enter entry into hash table ...\n");
      exit(0);
    }     
    
    // len
    chrlen[numchr] = atoi(p[1]);

    chr   [numchr] = strdup(p[0]);

    // alloc binary array for chr
    //map[numchr]    = create_binarized_array(chrlen[numchr]);

    numchr++;
    if (numchr == 50) {
      die("max # chromosomes reached.\n");
    }

  }
  fprintf(stderr, "Done reading chrlen\n");
  
  // load Ns
  fn = fopen(argv[2], "r");
  j = 0;  
  while (!feof(fn)) {

    fgets(buff, mynmax, fn);
    if (feof(fn))
      break; 
    chomp(buff);
    split_line_delim(buff, "\t", &p, &m);

    e.key   = p[0];
    my_hsearch_r(e, FIND, &ep, hash_chr);

    if (ep) {

      int id = (int)(ep->data);
      int st = atoi(p[1]);
      int en = atoi(p[2]);
      int le = chrlen[id];
      for (i=max(0,st-k); i<=min(le-1,en+k); i++) {
	nummap[id]++;
	//set_entry_in_binarized_array(map[id], i);
      }

    } else {
      die("cannot find chr in Ns\n");
      
    }

    j++;	    
    //if (j % 100000 == 0) {
     // fprintf(stderr, "read %d lines              \r", j);
    //}
   free(p);

  }
  fprintf(stderr, "Done reading where Ns are\n");

  // load k-mer counts ... the big file.
  fk = fopen(argv[3], "r");  
  j = 0;
  int id;
  int po;
  int co;
  while (!feof(fk)) {

    fgets(buff, mynmax, fk);
    if (feof(fk))
      break; 
    chomp(buff);
    split_line_delim(buff, "\t", &p, &m);

    // that's assuming that the chr id is the same as in wg.fa
    id = atoi(p[0]);
    po = atoi(p[1]);
    co = atoi(p[2]);
    
    if (po > chrlen[id]+10) {
      printf("Problem: pos in readmap = %d exceeds chr len = %d\n", po, chrlen[id]);
      exit(1);
    }
    if (co > 1) 
      nummap[id] ++;
    
    //for (i=po; i<po+k; i++) {
    //set_entry_in_binarized_array(map[id], i);
    //}

    
    free(p);
    j++;
    //if (j % 1000000 == 0) {
      //fprintf(stderr, "read %d lines              \r", j);
    //}
    
  }

  //long nummappable = 0;
  //long total       = 0;
  for (i=0; i<numchr; i++) {
    //printf("%s\tNon-mappable\t%d\tTotal\t%d", chr[i], nummap[i], chrlen[i]);
    printf( "%s\t%d\t%4.3f\n", chr[i], chrlen[i], 1.0 - nummap[i]/(double)chrlen[i] );
    //for (j=0; j<chrlen[i]; j++) {
    //if (get_entry_in_binarized_array(map[i], j) == 0) {
    //nummappable ++;
    //}
    // total++;
    //}	
  }

  //printf("num bases mappable = %ld out of %ld\n", nummappable, total);
  
  return 0;
}
