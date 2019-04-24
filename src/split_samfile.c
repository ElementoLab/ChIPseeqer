//
// implement buffered splitting of SAM files
//
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistics.h"
#include "dataio.h"
#include <math.h>
#include <search.h>
#include <ctype.h>
#include <errno.h>
#include "hashtable.h"
#include "sequences.h"

#define MAXNUMENTRIES 100000
#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

typedef struct _Line {
  char* line;
  struct _Line* next;
  int num;
} Line;
/* linked list
 
   1. array of Line pointers ?
   Line** lines
   lines[i] = 0;
 
   2. insert new pointer
   tmpline  = lines[i];
   lines[i] = new;
   lines[i]->next = tmpline; 
*/

int main(int argc, char** argv) {
	
  FILE** a_outfiles = 0;
	
  int   maxnumfiles = 100000;
	
  int   bufflen = 10000000;
  char** p;
  int   m;
  int   i =0;
  int    hashret;
  ENTRY  e;
  ENTRY* ep;
  int idxfile = -1;
  int cntfile = 0;
  FILE* thisfp;
  Line* lineptr;
  Line* curlineptr;
  Line* newline;
  Line* tmpline;
  int   numwritten = 0;
  int   numlinespergene = 0;
  char*	outdir = 0;
  struct my_hsearch_data* hash_genes;
  int verbose = 0;
  
  if (argc < 2) {
    printf("Usage: split_samfile FILE1 FILE2 ... [ -outdir DIR] \n");
    exit(0);
  }
	
	
  if (exist_parameter(argc, argv, "-outdir")) {
    outdir  = get_parameter(argc, argv, "-outdir");		
    mkdir(outdir, S_IRWXU | S_IRWXO | S_IRWXG);		
  } else {
    die("Please specify -outdir\n");
  }

  if (exist_parameter(argc, argv, "-verbose")) {
    verbose = atoi(get_parameter(argc, argv, "-verbose"));		
    printf("# verbose mode\n");
  }
	
  char* outfile = (char*)calloc(10000, sizeof(char));
	
  char** chrnames = (char**)calloc(maxnumfiles, sizeof(char*));
  a_outfiles = (FILE**)malloc(maxnumfiles * sizeof(FILE));
	
  // create an array of line pointers
  Line** a_lines    = (Line**)calloc(maxnumfiles, sizeof(Line*));
  int*   a_numlines = (int*)  calloc(maxnumfiles, sizeof(int));
	
  /* build a hash table */
  hash_genes = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(maxnumfiles, hash_genes);
  if (hashret == 0) {
    printf("Could not create hash table ...\n");
    exit(0);
  }
  cntfile = 0;
	
  int cntlines = 0;
  int cntsam   = 1;
  char* lco;
  // iterate thru sam files
  while (cntsam <= argc) { 

    if (strcmp(argv[cntsam], "-outdir") == 0) {
      // we have reached outdir, stop
      break;
    }

    LineI li;
    char* l;
    
    li.bufflen = bufflen;
    if (!LineIopen(&li, argv[cntsam])) {
      printf("Could not open %s\n", argv[cntsam]);
      exit(1);
    }
    
    printf("# Opened %s\n", argv[cntsam]);

    while ((l = nextLine(&li))) {

      if (l[0] == '@') {
	free(l);
	continue;
      } 

      lco = strdup(l);  // make a copy
      split_line_delim(lco, "\t", &p, &m); // split it 
      
      if (strcmp(p[5], "*") == 0) {
	free(p);
	free(l);
	free(lco);
	continue;
      } 
      
      // inc line num
      cntlines ++;

      // chr name 
      char* chrname = p[2];
      
      // we don't need p anymore
      

      /* query chr ... do we have it already */
      e.key = chrname;        
      my_hsearch_r(e, FIND, &ep, hash_genes) ;
      thisfp = 0;
      
      if (ep) {
	/* found existing file */
	idxfile = (int)(ep->data);
	
      } else {
	
	/* enter key/value pair into hash */
	e.key   = strdup( chrname );
	e.data  = (char*)cntfile;
	hashret = my_hsearch_r(e, ENTER, &ep, hash_genes);
	if (hashret == 0) {
	  printf("Could not enter entry into hash table ...\n");
	  exit(0);
	}
	idxfile = cntfile;
	
	// add to chr db
	chrnames[cntfile] = strdup(chrname);					
	cntfile++;
					
	if (cntfile == maxnumfiles) {
	  die("max number of sequences reached. You need to increase the hard limit and recompile.\n");
	}
      } // else
      
    
      // add line to linked list at idxfile
      newline            = (Line*)malloc(sizeof(Line));  // get new line
      tmpline            = a_lines[ idxfile ];           // then add it
      newline->next      = tmpline;
      newline->num       = cntlines;
      newline->line      = l;
      a_lines[ idxfile ] = newline;
      a_numlines[ idxfile ] ++;                          // increase num
      
      if ( (cntlines % 100000) == 0) {                   // if 100000 inc reached
	
	if (verbose == 1) {
	  printf("# %d lines reached, write\n", cntlines);
	  fflush(stdout);
	}

	numlinespergene = 0;
	// iterate thru existing chr files
	for (i=0; i<cntfile; i++) {

	  // if we have lines to write
	  if (a_lines[i] != 0) {
	    
	    // make file name
	    if (outdir != 0) 
		sprintf(outfile, "%s/reads.%s", outdir, chrnames[i]);
	      else
		sprintf(outfile, "reads.%s", chrnames[i]);
	    
	    // make hadler
	    thisfp = fopen(outfile, "a");
	    if (thisfp == 0) {
	      printf("Cannot open file %s for appending\n", outfile);
	      exit(1);
	    }
	    
	    if (verbose == 1) {
	      printf("# about to write %d lines to %s\n", a_numlines[i], outfile);
	    }

	    // iterate thru lines
	    lineptr   = a_lines[i];
	    curlineptr = 0;
	    numlinespergene = 0;
	    while (lineptr != 0) {
	      fprintf(thisfp, "%s\n", lineptr->line);
	      numlinespergene++;
	      numwritten ++;
	      if (numwritten % 100000 == 0) {
		printf("Wrote %d lines.          \r", numwritten);
		fflush(stdout);
	      }
	      curlineptr = lineptr;
	      lineptr    = lineptr->next;
	      
	      // destroy current, we don't need it anymore
	      free(curlineptr->line);
	      free(curlineptr);
	    }
	    a_lines[i] = 0; // set init pointer to 0
	    fclose(thisfp); // close fp
	
	    // quick sanity check
	    if (numlinespergene != a_numlines[i]) {
	      printf("Problem(1): %d lines stored, %d lines written\n", a_numlines[i], numlinespergene);
	      exit(1);
	    }
							
	  } // if (a_lines[i] != 0
						
	} // for (i<cntfile
					
	for (i=0; i<cntfile; i++) // redundant
	  a_numlines[i] = 0;
	
      } // if 100000
      
      free(lco);  // free line copy 
      free(p);    // free p

      //cntlines++; // inc lines

    } // end line iterator

    LineIclose(&li); // close iterator
  
    // write one more time; we still have lines in buffer
    numlinespergene = 0;
    for (i=0; i<cntfile; i++) {
      if (a_lines[i] != 0) {
	
	if (outdir != 0) 
	  sprintf(outfile, "%s/reads.%s", outdir, chrnames[i]);
	else
	  sprintf(outfile, "reads.%s", chrnames[i]);
	
	thisfp = fopen(outfile, "a");
	if (thisfp == 0) {
	  printf("Cannot open file %s for appending\n", outfile);
	  exit(1);
	}
	lineptr   = a_lines[i];
	curlineptr = 0;
	numlinespergene = 0;
	while (lineptr != 0) {
	  fprintf(thisfp, "%s\n", lineptr->line);
	  numlinespergene++;
	  numwritten ++;
	  if (numwritten % 100000 == 0) {
	    printf("Wrote %d lines.            \r", numwritten);
	    fflush(stdout);
	    
	  }
	  curlineptr = lineptr;
	  lineptr = lineptr->next;
	  
	  // destroy current
	  free(curlineptr->line);
	  free(curlineptr);
	}
	a_lines[i] = 0;
	//printf("\n");
	fclose(thisfp);
	
	if (numlinespergene != a_numlines[i]) {
	  printf("Problem(2): %d lines stored, %d lines written\n", a_numlines[i], numlinespergene);
	  exit(1);
	}
	a_numlines[i] = 0;
      } // if (a_lines[i] != 0
      
    }
    
    // we have reached the end of the processing of 1 sam file
    printf("# Finished processing %s\n", argv[cntsam]);
    cntsam++;
    
  }
	
  return 0;
	
}
