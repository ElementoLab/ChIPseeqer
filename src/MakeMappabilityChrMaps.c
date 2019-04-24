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

	FILE* fk		= NULL;
	char* buff		= NULL;
	int   mynmax	= 10000;
	char** p		= NULL;
	int    m		= 0;
  
  // general
  int    i = 0;
  int    k = 30;

  //int*   nummap;
  
  float*	chrmap		= 0;
  int		numchroms	= 0;
  char*		chrdata		= 0;
  char**	chrnames	= 0;
  int*		chrlens		= 0;
  char*		readmapfile = 0;
  char*		outdir		= 0;

  if (argc == 1) {
    die("Args: -chrdata FILE -readmap FILE -outdir DIR \n");
  }
  
  if (exist_parameter(argc, argv, "-chrdata")) {
    chrdata  = get_parameter(argc, argv, "-chrdata");
    readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);
    //genome   = 0;
  }
  
  if (exist_parameter(argc, argv, "-readmapfile")) {
    readmapfile = get_parameter(argc, argv, "-readmapfile");
  }
  
  if (exist_parameter(argc, argv, "-outdir")) {
    outdir = get_parameter(argc, argv, "-outdir");
  }


  int chrnum = -1;
  unsigned char* mapmap = 0;
  int kmercnt;
  int chrpos;
  int kmersize;
  int prevchrnum = chrnum;
  char chromfile[1000];
  FILE* fo ;

  // load k-mer counts ... the big file.
  fk = fopen(readmapfile, "r");  
  if (!fk) {
    die("Cannot read readmapfile\n");
  }
  
  buff = (char*)calloc(1000, sizeof(char));

  while (fgets(buff, mynmax, fk) != 0) { 
    chomp(buff);

    split_line_delim(buff, "\t", &p, &m);

    // that's assuming that the chr id is the same as in wg.fa
    chrnum   = atoi(p[0]);
    chrpos   = atoi(p[1]);

    if (chrpos > chrlens[chrnum]+10) {
      fprintf(stderr, "Problem: pos in readmap = %d exceeds chr len = %d\n", chrpos, chrlens[chrnum]);
      //exit(1);
    }

    kmercnt  = atoi(p[2]);
    kmersize = strlen(p[3]);

    if (chrnum != prevchrnum) {

      if (prevchrnum != -1){ 
      
	// output      
	sprintf(chromfile, "%s/map.%s", outdir, chrnames[prevchrnum]);
	fo = fopen(chromfile, "w");
	if (fo == 0) {
	  printf("Can't open %s\n", chromfile);
	  exit(1);
	}
	
	for (i=0; i<chrlens[prevchrnum]; i++) {
	  int score = 0;
	  for (k=max(0,i-kmersize+1); k<=i; k++) {
	    score += get_entry_in_binarized_array(mapmap, k);
	  }
	  fprintf(fo, "%d\n", score);
	  //printf("%s\t%d\n", chrnames[prevchrnum], score);
	}
	fclose(fo);

	// free
	free(mapmap);
	
      }
	
      // alloc new one
      mapmap = create_binarized_array(chrlens[ chrnum ]);
    

    }
    
    if (kmercnt == 1) { // mappable
      
      for (k=chrpos; k<chrpos+kmersize; k++) {
	set_entry_in_binarized_array(mapmap, k);
      }

    }

    prevchrnum = chrnum;
    
    free(p);
    
  }

  //
  // output      
  //
  sprintf(chromfile, "%s/map.%s", outdir, chrnames[prevchrnum]);
  fo = fopen(chromfile, "w");
  if (fo == 0) {
    printf("Can't open %s\n", chromfile);
    exit(1);
  }
  
  for (i=0; i<chrlens[prevchrnum]; i++) {
    int score = 0;
    for (k=max(0,i-kmersize+1); k<=i; k++) {
      score += get_entry_in_binarized_array(mapmap, k);
    }
    fprintf(fo, "%d\n", score);
    //printf("%s\t%d\n", chrnames[prevchrnum], score);
  }
  fclose(fo);
  
  // free
  free(mapmap);
  
  return 0;
}
