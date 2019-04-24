#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <search.h>

#include "dataio.h"

#include "statistics.h"
#include "hashtable.h"
#include "sequences.h"



#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)
char* this_complement(char* s); 

char uc_char(char c);

int* N;
int* C;

int main(int argc, char ** argv)
{

  seqI     si;
  char*    fastafile = 0;

  int      size;
  char*    name;
  int      i, l;

  int*     matches_pos    = 0;
  char*    matches_ori    = 0;
  int      np=0;
  int      j = 0;
  char*    seq1;
  int      verbose = 0;
  
  float**  fwm;
  float**  wm;
  char**   sites;
  int      w;
  int      nsites;
  int*     nsites_j;

  int*     stars;
  float    t = -100000;

  int**    intwm;

  //float   human_bkg_raw[4] = {0.23, 0.27, 0.23, 0.27};
  //float   human_bkg[4]     = {-1.47, -1.31, -1.47, -1.31};

  float bkg[4];

  char*   asite;


  char*   motiffile  = 0;
  char*   jasparfile = 0;
  char*   bulykfile  = 0;
  float   g          = 0.5;
  float   c = 2.0;
  float*  matches_sco;
  int     p = 0;
  char*   mn = 0;
  
  // geneset hash
  ENTRY   e;
  ENTRY*  ep;
  int     hashret;

  char*   gs = 0;
  char**  genenames;
  int     numgenes = 0;
  struct my_hsearch_data* hash_genes = 0;

  int     rd = 0;
  int     pos;
  int     to = 0;
  char*   wmscorefile = 0;
  int     format = 0;
  char*   formattext = 0;
  float   minscore, maxscore;
  int     allowmatchoverlap = 1;
  int     maxnumntdiffwithcons = -1;
  char*   consensus = 0;
  int     k = 0;
  int     maxnummatches = 100000;
  char*   intfile = 0;
  int***   intervals = 0;
  int*     numints = 0;
  int     cidx = -1;
  char*   chr = 0;
  int     header = 1;

  if (argc == 1) {
    printf("Usage: MyScanACE -i,j,jb motifile -z fastafile -g gcbackround -c numstddevbelowavg -p [0/1]/2 -gs FILE -wmscorefile FILE -allowmatchoverlap INT -output gmod|counts|cs|wig\n");
    exit(0);
  }

  default_set_seed(122);
  
  if (exist_parameter(argc, argv, "-i"))
    motiffile  = get_parameter(argc, argv, "-i");
  else if (exist_parameter(argc, argv, "-j")) {
    jasparfile = get_parameter(argc, argv, "-j");
  } else if (exist_parameter(argc, argv, "-jb")) {
    bulykfile  = get_parameter(argc, argv, "-jb");
  } else
    die("Please set -i or -j (motif file).\n");

  if (exist_parameter(argc, argv, "-z"))
    fastafile = get_parameter(argc, argv, "-z");
  else 
    die("Please set -z (fasta file).\n");

  if (exist_parameter(argc, argv, "-g"))
    g = atof(get_parameter(argc, argv, "-g"));

  if (exist_parameter(argc, argv, "-header"))
    header = atoi(get_parameter(argc, argv, "-header"));

  if (exist_parameter(argc, argv, "-chr"))
    chr = get_parameter(argc, argv, "-chr");

  if (exist_parameter(argc, argv, "-mn"))
    mn = get_parameter(argc, argv, "-mn");
  
  if (exist_parameter(argc, argv, "-gs"))
    gs = get_parameter(argc, argv, "-gs");

  if (exist_parameter(argc, argv, "-wmscorefile"))
    wmscorefile = get_parameter(argc, argv, "-wmscorefile");

  if (exist_parameter(argc, argv, "-maxnummatches"))
   maxnummatches = atoi(get_parameter(argc, argv, "-maxnummatches"));

  if (exist_parameter(argc, argv, "-maxnumntdiffwithcons"))
    maxnumntdiffwithcons = atoi(get_parameter(argc, argv, "-maxnumntdiffwithcons"));

  if (exist_parameter(argc, argv, "-allowmatchoverlap"))
    allowmatchoverlap = atoi(get_parameter(argc, argv, "-allowmatchoverlap"));

  if (exist_parameter(argc, argv, "-amo"))
    allowmatchoverlap = atoi(get_parameter(argc, argv, "-amo"));
 
  if (exist_parameter(argc, argv, "-maxnummatches"))
    maxnummatches = atoi(get_parameter(argc, argv, "-maxnummatches"));

  if (exist_parameter(argc, argv, "-intervals")) {
    intfile = get_parameter(argc, argv, "-intervals");
  }

  if (exist_parameter(argc, argv, "-t"))
    t = atof(get_parameter(argc, argv, "-t"));

  if (exist_parameter(argc, argv, "-c"))
    c = atof(get_parameter(argc, argv, "-c"));

  if (exist_parameter(argc, argv, "-p"))
    p = atoi(get_parameter(argc, argv, "-p"));
 
  if (exist_parameter(argc, argv, "-to"))
    to = atoi(get_parameter(argc, argv, "-to"));
  
  if (exist_parameter(argc, argv, "-rd"))
    rd = atoi(get_parameter(argc, argv, "-rd"));
  
  if (exist_parameter(argc, argv, "-v"))
    verbose = atoi(get_parameter(argc, argv, "-v"));

  if (exist_parameter(argc, argv, "-output")) {
    formattext = get_parameter(argc, argv, "-output");
    if (strcmp(formattext, "gmod") == 0)
      format = 1;
    else if (strcmp(formattext, "wig") == 0)
      format = 2;
    else if (strcmp(formattext, "counts") == 0) 
      format = 3;	
    else if (strcmp(formattext, "density") == 0) 
      format = 4;	
    else if (strcmp(formattext, "cs") == 0)
      format = 10;
    if (verbose == 1) 
      fprintf(stderr, "# formt = %d\n", format);
  }

  if ((format == 3) && (p == 0)) 
    p = 1;
    
  
  bkg[0] = (1.0 - g) / 2.0;
  bkg[1] =        g  / 2.0;
  bkg[2] = (1.0 - g) / 2.0;
  bkg[3] =        g  / 2.0;

  initialize_nt(); 

  if (motiffile != 0) {
    //
    // load up AlignACE motif						
    //
    
    // read in count matrix
    readACEintWM(motiffile, &intwm, &w, &sites, &nsites, &stars);
    
    // convert to log2 bkg-corrected
    ACEintWMtologWM(intwm, w, nsites, bkg, &wm);

    if (verbose == 1)
      printWM(wm, w);

   
  } else if (jasparfile != 0) {
    

    

    // read in count matrix
    read_JASPAR_WM(jasparfile, &intwm, &w, &nsites_j);
    
    if (verbose == 1) {
      printf("w=%d\n", w);
      for (i=0; i<w; i++) {
	//printf("n[col=%d]=%d\n", i, nsites_j[i]);
      }
    }
    
    stars = (int*)malloc(w * sizeof(int));
    for (i=0; i<w; i++) {
      stars[i] = 1;
    }
 
    // convert to log2 bkg-corrected
    ACEintWMtologWM_J(intwm, w, nsites_j, bkg, &wm);
 
   
    if (verbose == 1) {
      printIntegerWM(intwm, w);
      printWM(wm, w);

    }
	  
  } else if (bulykfile != 0) {

    if ( !file_exists(bulykfile) )
      die("Motif file does not exist.\n");
    
    // read in MB matrix
    read_BULYK_WM(bulykfile, &fwm, &w);
    
    // convert to log
    // ACEintWMtologWM(intwm, w, nsites, bkg, &wm);
    WMtologWMbkg(fwm, w, bkg, &wm); 
    //printWM(wm, w);

    stars = (int*)malloc(w * sizeof(int));
    for (i=0; i<w; i++) {
      stars[i] = 1;
    }

    if (verbose == 1)
      printWM(wm, w);


  } else {
    die("Please input a motif using -i, -j or -jb.\n");
  }

  

  if (!exist_parameter(argc, argv, "-t") && (p != 2)) {

    //if (jasparfile != 0)  {
    //  die("Cannot calculate a threshold this way with Jaspar matrices.\n");
    //}

    if (motiffile != 0) {
      // get a threshold (2 means that log2 bkg-corrected matrix is assumed)
      t = getWMThreshold(wm, w, stars, 0, 2, sites, nsites, c, &minscore, &maxscore);

      if (wmscorefile != 0)
	writeWMScoresToFile(wm, w, stars, 0, 2, sites, nsites, wmscorefile); 


    } else if (bulykfile != 0) {

      // Jaspar or Bulyk
      // 1. generate       
      generateManyRandomSequences(fwm, w, 10000, &sites);
      t = getWMThreshold(wm, w, stars, 0, 2, sites, 10000, c, &minscore, &maxscore);
      
      if (wmscorefile != 0)
	writeWMScoresToFile(wm, w, stars, 0, 2, sites, 10000, wmscorefile); 

      //exit(0);

    } else if (jasparfile != 0) {
      
      // need to generate a fwm from intwm
      intWMtofloatWM_J(intwm, w, nsites_j, &fwm); 

      generateManyRandomSequences(fwm, w, 100000, &sites);
      t = getWMThreshold(wm, w, stars, 0, 2, sites, 100000, c, &minscore, &maxscore);
      
      if (wmscorefile != 0)
	writeWMScoresToFile(wm, w, stars, 0, 2, sites, 100000, wmscorefile); 

      //die("No Jaspar support for this yet.\n");
    }

    if (to == 1) {
      //printf("Threshold = %f\n", t);
      //exit(0);
    }
  }

  if (maxnumntdiffwithcons >= 0) {
    
    // get consensus
    getWMConsensus(fwm, w, &consensus);
    //printf("consensus=%s\n", consensus);
    
  }

  if (gs != 0) {    
    genenames = readStringSet(gs, &numgenes);

    // to index genes
    hash_genes  = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
    hashret = my_hcreate_r(100000, hash_genes);
    if (hashret == 0) {
      printf("Could not create hash table ...\n");
      exit(0);
    } 
    
    // enter genes
    for (i=0; i<numgenes; i++) {
      //printf("Entering %s\n", genenames[i]);
      e.key   = strdup(genenames[i]);
      e.data  = (char*)1;
      hashret = my_hsearch_r(e, ENTER, &ep, hash_genes);
      if (hashret == 0) {
	printf("Could not enter entry into hash table ...\n");
	exit(0);
      }     
    }
    
  }

  HASH hc;
  int  numchr = 0;
  

  if (intfile != 0) {

    HASH_init(&hc, 100);
    numints = (int*)calloc(100, sizeof(int));

    char** p;
    int  m;
    FILE* fi = fopen(intfile, "r");
    if (fi == 0) {
      die("Cannot open intervals\n");
    }

    char* buff = (char*)calloc(10000, sizeof(char));
    cidx = -1;
    while (fgets(buff, 10000, fi) != 0) {
      chomp(buff);
      split_line_delim(buff, "\t", &p, &m);
      if (!HASH_find(&hc, p[0], &cidx)) {
	HASH_insert(&hc, p[0], numchr);
	cidx = numchr;
	numchr++;
      }
      numints[cidx] ++;
      free(p);      
    }
  
    // alloc memory
    intervals = (int***)calloc(numchr, sizeof(int**));
    for (i=0; i<numchr; i++) {
      intervals[i] = (int**)calloc( numints[i], sizeof(int*));
      for (j=0; j<numints[i]; j++) {
	intervals[i][j] = (int*)calloc(2, sizeof(int));
      }
    }

    // start over
    rewind(fi);
    
    for (i=0; i<numchr; i++) {
      if (verbose == 1)
	printf("Read %d intervals for chr %d\n", numints[i], i);
      numints[i] = 0;
    }

    while (fgets(buff, 10000, fi) != 0) {
      chomp(buff);
      split_line_delim(buff, "\t", &p, &m);
      if (!HASH_find(&hc, p[0], &cidx)) {
	die("Mkes no sense\n");
      }      
      intervals[ cidx ][ numints[cidx] ][0] = atoi(p[1]);
      intervals[ cidx ][ numints[cidx] ][1] = atoi(p[2]);
      numints[cidx] ++;
      free(p);      
    }

    
    
    fclose(fi);

  }
  

  //						
  // get first sequence (ref)
  //
  if (verbose == 1)
    printf("Reading sequences.\n");

  if ((format == 2) && (mn != 0)) {
    printf("track name=\"%s\" description=\"%s\" visibility=dense autoScale=off color=0,0,0\n", mn, mn);
  }
  
  int numtgtgenes = 0;
  int numseqs     = 0;
  int effseqlen   = 0;

  if ((p == 1) && (format != 10) && (header == 1)) {
    printf("GENE\tMOTIF\n");
  }

  // working for the human genome (fast)
  seqI_set_max_seqlen(&si, 500000000UL);
  seqI_set_seqlen_inc(&si, 300000000UL);
  seqI_set_fread_chunk(&si, 100000000UL);
  if (seqI_open_fast(&si, fastafile) == 0) {
    die("Error opening file\n");
  }

  while ((seq1     = seqI_nextSequence_fast(&si, &name, &size))) {

    if ((chr != 0) && (strcmp(chr, name) != 0)) {
      free(seq1);
      free(name);
      continue;
    }

    if (intfile != 0) {
      if (!HASH_find(&hc, name, &cidx)) {
	free(seq1);
	free(name);
	continue;
      }
    }

    if (gs != 0) {      
      e.key   = name;
      my_hsearch_r(e, FIND, &ep, hash_genes);
      if (!ep) {
	free(seq1);
	continue;
      }      
    }

    l        = strlen(seq1);  
    if (verbose == 1) {
      printf("Seq Size = %d\n", l);
      fflush(stdout);
    }
    
    for (j=0; j<l; j++) {
      seq1[j] = toupper(seq1[j]);
    }

    // calculate effective length if needed
    if (format == 4) {
      effseqlen = 0;
      for (j=0; j<l; j++) {
	if (seq1[j] != 'N') 
	effseqlen ++;
      }
    }

    int runfind = 1;
    if (p < 2) {

      np = 0;
      if (intfile == 0) {
	findAllWeightMatrixMatches(wm, w, stars, 0, t, seq1, 0, &matches_pos, &np, &matches_ori, &matches_sco, maxnummatches);
      } else {
	cidx = -1;
	if (HASH_find(&hc, name, &cidx))       
	  findAllWeightMatrixMatchesInIntervals(wm, w, stars, 0, t, seq1, 0, &matches_pos, &np, &matches_ori, &matches_sco, maxnummatches, intervals[cidx], numints[cidx]);
	else 
	  runfind = 0;
      }

      if (verbose == 1)
	printf("Found %d sites\n", np);

      if (p == 0) {

	if (format == 1)
	  printf("reference = %s\n", name);
	
	int prevmatchpos = -1000;
	for (i=0; i<np; i++) {

	  // if overlaps between matches are not allowed and pres match overlaps with prev, continue
	  if ((allowmatchoverlap == 0) && (matches_pos[i] < prevmatchpos+w)) {
	    
	    continue;
	  }


	  asite = substr(seq1, matches_pos[i], w);      

	  // if min edit distance with consensus is specified ... 
	  if (maxnumntdiffwithcons >= 0) {
	    int nummm = basicEditDistance(asite, consensus);
	    if (nummm > maxnumntdiffwithcons) {
	      free(asite);
	      continue;
	    }
	  }

	  prevmatchpos = matches_pos[i];

	  if (format == 0) {
	    printf("%s\t", name);
	    if (mn != 0) {
	      printf("%s\t", mn);
	    }
	  }

	  pos = matches_pos[i];
	  if (rd == 1)
	    pos =  - ( l - (pos + 1));
	  
	  if (format == 0)
	    printf("%d\t%d\t%3.2f\t%s\n", pos, matches_ori[i], matches_sco[i], asite);
	  else if (format == 1)
	    printf("MotifMatch\t%s\t%d-%d\tscore=%d\n", (mn!=0?mn:bulykfile), pos, pos+w, (int)(0.5+matches_sco[i]));
	  else if (format == 2)
	    printf("%s\t%d\t%d\n", name, pos, pos+w);

	  free(asite);	
	}
      } else if (p == 1) {

	if (format == 0) {

	  printf("%s\t%d\n", name, (np>0?1:0) );
	  if (np > 0)
	    numtgtgenes ++;

        } else if (format == 3) {
	  //counts
	  printf("%s\t%d\n", name, np);
	   
	} else if (format == 4) {
	  // density
	  printf("%s\t%3.2f\n", name, (np==0?0:1000.0*np/(double)effseqlen));

	} else if (format == 10) {
	  // special format for ChIPseeqer
	  if (np > 0) {
	    char** p;
	    int    m;	    
	    split_line_delim_sep(name, "/", &p, &m);
	    printf("%s", p[0]);
	    for (k=1; k<m; k++) {
	      printf("\t%s", p[k]);
	    }


	    // add position of best match ?
	    int posbest = -1;
	    int scobest = -100.0;
	    for (i=0; i<np; i++) {
	      if (matches_sco[i] > scobest) {
		posbest = matches_pos[i];
		scobest = matches_sco[i];
	      }
	    }
	    
	    int actpos = posbest+w/2;
	    printf("\t%d", actpos);

	    printf("\t%3.1f", 100*actpos/(float)l);

	    printf("\n");

	    free(p);
	    numtgtgenes ++;
	  }
	  

	}
      }
    
    } else if (p == 2) {
      
      
      int posidx = -1;
      float maxscore = findMaxWMScore(seq1, wm, 0, w, stars, bkg, 0, &posidx); 
      
      asite = substr(seq1, posidx, w);      

      pos = posidx;
      if (rd == 1)
	pos =  - ( l - (pos + 1));

      printf("%s\t%4.3f\t%d\t%s\n", name, maxscore, pos, asite );
                   
    }
  
    numseqs ++;
    
    if (runfind == 1) {
      free(matches_pos);
      free(matches_sco);
      free(matches_ori);
    }
    free(seq1);
    
  }
    
  seqI_close(&si);

  if (to == 1) {
    printf("threshold=%f, #target genes = %d\n", t, numtgtgenes);
  }

  if ((p == 1) && (format != 10)) {
    printf("# %d peaks, %d with motif match (%3.1f %%)\n", numseqs, numtgtgenes, 100*numtgtgenes/(float)numseqs);
  }
  
  return 0;
}


char uc_char(char c) {
  if (c == 'a') {
    return 'A';
  } else if (c == 'c') {
    return 'C';
  } else if (c == 't') {
    return 'T';
  } else if (c == 'g') {
    return 'G';
  } else {
    return c;
  }
}




char* this_complement(char* s) 
{

  int   l = strlen(s);
  char* c;
  char  d;
  int   i;


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
    else if (s[i] == '-')
      d = '-';
    else
      d = 'N';

    strncat(c, &d, 1);
  }
	 
  return c;
  
}
