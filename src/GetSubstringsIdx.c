#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <search.h>

#include "dataio.h"
#include "prefix.h"
#include "sequences.h"
#include "statistics.h"
#include "hashtable.h"

#include "lib/third_party/samtools/sam.h"
#include "lib/third_party/samtools/faidx.h"

#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)

typedef struct _GenInt {
	int i;
	int j;
	char* c;
	char* id;
	char*  dbsnp;
	char* str;
} GenInt;

char* my_fai_fetch(faidx_t* idx, char* c, int i, int j);

int main(int argc, char ** argv)
{
	
	char*    fastafile	= NULL;  
	int      i			= 0;
	char*    intervals	= NULL;
	char*    buff		= NULL;
	FILE*    f			= NULL;
	int      mynmax		= 100000;
	char**   p			= NULL;
	int      m			= 0;
	
	int      hasid     = 0;
	int      randomize = 0;
	int      minlen    = -1;
	int      forcelen  = -1;
	int      extend    = -1;
	int      iswig     = 0;
	char*    chr             = 0;
	int      seed            = 2345;
	int      useall          = 0;
	char*    h               = 0;
	int      header          = 0;  
	int      cidx            = 0;
	int	  rand_close	  = 0;
	int	  distance	  = 10000;
	int	  show_clust	  = 0;
	int	  print_both	  = 0;
	int      show_peakdesc   = 0;
	int		 show_seq_only	= 0;
	HASH     clh;
	int      verbose         = 0;
	char*    genome          = "hg18"; 
	int      output          = 0;
	char*    output_txt      = 0;
	char*    chrdata         = 0;
	int      numchroms       = 0;
	int*     chrlens         = 0;
	char**   chrnames        = 0;
	float*   chrmaps         = 0;
	int      chrlen          = 0;
	int      idx             = 0;
	int      actextleft      = -1;
	int      actextright     = -1;
	int      newen           = -1;
	int      newst           = -1;
	int      blankpeaks      = -1;
	
	if (argc == 1) {
		die("Args: -peakfile FILE -fastafile FILE [ -extend INT -forcelen INT -minlen INT ] \n");
	}
	
	if (exist_parameter(argc, argv, "-fastafile"))
		fastafile  = get_parameter(argc, argv, "-fastafile");
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals  = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-peakfile"))
		intervals  = get_parameter(argc, argv, "-peakfile");
	
	if (exist_parameter(argc, argv, "-hasid"))
		hasid  = atoi(get_parameter(argc, argv, "-hasid"));
	
	if (exist_parameter(argc, argv, "-header"))
		header  = atoi(get_parameter(argc, argv, "-header"));
	
	if (exist_parameter(argc, argv, "-seed"))
		seed  = atoi(get_parameter(argc, argv, "-seed"));
	
	if (exist_parameter(argc, argv, "-useall"))
		useall  = atoi(get_parameter(argc, argv, "-useall"));
	
	if (exist_parameter(argc, argv, "-randomize"))
		randomize  = atoi(get_parameter(argc, argv, "-randomize"));
	
	if (exist_parameter(argc, argv, "-minlen"))
		minlen = atoi(get_parameter(argc, argv, "-minlen"));
	
	if (exist_parameter(argc, argv, "-forcelen"))
		forcelen = atoi(get_parameter(argc, argv, "-forcelen"));
	
	if (exist_parameter(argc, argv, "-extend"))
		extend = atoi(get_parameter(argc, argv, "-extend"));
	
	if (exist_parameter(argc, argv, "-blankpeaks"))
		blankpeaks = atoi(get_parameter(argc, argv, "-blankpeaks")); 
	
	if (exist_parameter(argc, argv, "-chr"))
		chr = get_parameter(argc, argv, "-chr");
	
	if (exist_parameter(argc, argv, "-iswig"))
		iswig = atoi(get_parameter(argc, argv, "-iswig"));
	
	if (exist_parameter(argc, argv, "-distance"))
		distance = atoi(get_parameter(argc, argv, "-distance"));
	
	if (exist_parameter(argc, argv, "-rand_close"))
		rand_close = atoi(get_parameter(argc, argv, "-rand_close"));
  	
	if (exist_parameter(argc, argv, "-show_clust"))
		show_clust = atoi(get_parameter(argc, argv, "-show_clust"));
	
	if (exist_parameter(argc, argv, "-print_both"))
		print_both = atoi(get_parameter(argc, argv, "-print_both"));
	
	if (exist_parameter(argc, argv, "-show_peakdesc"))
		show_peakdesc = atoi(get_parameter(argc, argv, "-show_peakdesc"));
	
	if (exist_parameter(argc, argv, "-genome"))
		genome = get_parameter(argc, argv, "-genome");
	
	if (exist_parameter(argc, argv, "-verbose"))
		verbose = atoi(get_parameter(argc, argv, "-verbose"));
	
	if (exist_parameter(argc, argv, "-output")) {
		output_txt = get_parameter(argc, argv, "-output");
		if (strcmp(output_txt, "tab") == 0) {
			output = 1;
		}   
	}
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmaps, &numchroms);
	} else {
		get_hg18_chroms(&chrnames, &numchroms, &chrlens);
	}
	
	if (exist_parameter(argc, argv, "-show_seq_only"))
		show_seq_only = atoi(get_parameter(argc, argv, "-show_seq_only"));

	
	// RNG init
	if (seed == 0) {
		default_set_seed(time(NULL));
	}
	else {
		default_set_seed(seed);
	}
	
	
	// store chromsomal indices  
	HASH_init(&clh, numchroms);
	for (i=0; i<numchroms; i++) {
		HASH_insert(&clh, chrnames[i], i);
	}
	
	char* id  = (char*)calloc(10000, sizeof(char));
	
	
	faidx_t* fidx = 0;
	
	// load index
	fidx = fai_load(fastafile);
	if (fidx == 0) {
		die("Cannot load fasta index\n");
	}
	
	char* peakid = (char*)calloc(10000, sizeof(char));
	buff         = (char*)malloc(mynmax * sizeof(char));
	
	//
	// load up intervals
	//
	f = fopen(intervals, "r");
	if (f == 0)
		die("cannot open interval file.\n");    
	if ((iswig == 1) || (header == 1)) {
		fgets(buff, mynmax, f);
		h = strdup(buff);
		chomp(h);    
		if (output == 1) {
			//printf("%s\tContext\n", h);
		}
	}
	while (fgets(buff, mynmax, f) != 0) {
		chomp(buff);
		if (buff[0] == '#')
			continue;
		char* line = strdup(buff);
		int linelen = strlen(line);
		for (i=0; i<linelen; i++) {
			if (line[i] == '\t')
				line[i] = '/';
		}
		split_line_delim(buff, "\t", &p, &m);
		
		idx  = 0;    
		
		if (id != NULL)
			free(id);
			
		if (hasid == 1) {
			idx = 1;
			id = strdup(p[0]);
		} else {
			id = 0;      
		}
		
		char* c     =      p[idx  ];
		int   i     = atoi(p[idx+1]);
		int   j     = atoi(p[idx+2]);
		
		if (j < i) {
			printf("Problem with interval: j=%d < i=%d\n", j, i);
			exit(-1);
		}
		
		// if chr specified, keep or skip
		if ((chr != 0) && (strcmp(chr, c) != 0)) {
			continue;
		}
		
		// look up
		if (!HASH_find(&clh, c, &cidx)) {
			fprintf(stdout, "%s not found\n", c);
			free(p);
			if (hasid == 1) {
				free(id);
				id = NULL;
			}
			free(line);
			continue;
		}
		
		// get chrom length
		chrlen = chrlens[cidx];    
		
		// for dm9
		if (strcmp(genome, "dm9") == 0) {
			c = c + 3;
		}
		
		// next make/adjust peak coordinate 
		int st    = i;
		int en    = j;
		int le    = abs(en - st + 1);
		
		// extend to match min length
		if ((forcelen != -1) && (le < forcelen)) {
			int sideext = (int)( 0.5 + ( forcelen - le ) / 2 );
			st = max(0,st - sideext);
			en = min(en + sideext,chrlen-1);
		}    
		
		int orile = abs(en - st + 1);
		
		// if extension required
		if (extend != -1) { 
			
			newst = max(0, st - extend);     
			newen = min(en + extend,chrlen-1);
			
			actextleft  = abs(st - newst);
			actextright = abs(newen - en);
			
			if (verbose == 1)
				printf("%d %d\n", actextleft, actextright);
			
			if (verbose == 1)
				printf("was %d %d, is %d %d\n", st, en, newst, newen);
			
			st = newst;
			en = newen;
			
		}
		
		// if minlen
		if ((minlen >= 0) && (le < minlen)) {
			free(line);
			free(p);
			if (hasid == 1) {
				free(id);
				id = NULL;
			}	
			continue;
		}
		
		// update length
		le = abs(en - st + 1);
		
		
		// generate randomly located peak
		// if rand_close == 1, random peak is selected from local distance .... 
		if (randomize == 1) {	 
			
			while (1) {	
				// draw an interval starting point at random from chromosome
				st = rand_chr_interval2(c, le, chrlen);	
				
				// random scheme 1
				// draw local intervals !!!!!! this is wayyyyyyy sub-optimal; it's much easier to select an interval starting at a distance d from peak
				if (rand_close) {
					if(((st > j) && (st < j + distance)) || ((st < i) && (st > i-distance))){
						en = st + le;
						
						
						//char* ss =  uc(substr(seq, st, abs(en - st + 1)));
						char* ss =  my_fai_fetch(fidx, c, st, en);
						if (ss == 0) {
							printf("Problem: could not extract %s:%d-%d\n", c, st, en);
							exit(-1);
						}
						uctransform(ss);
						
						int cntN = 0;
						for (j=0; j<le; j++) {
							if (ss[j] == 'N') 
								cntN ++;	       
						}
						
						free(ss);
						
						if ( (cntN / (float)le) < 0.1) {
							break;
						} 
					}
				} // end rand_close
				else {
					en = st + le;
					
					//char* ss =  uc(substr(seq, st, abs(en - st + 1)));
					char* ss =  my_fai_fetch(fidx, c, st, en);
					if (ss == 0) {
						printf("Problem: could not extract %s:%d-%d\n", c, st, en);
						exit(-1);
					}
					uctransform(ss);
					
					int cntN = 0;
					for (j=0; j<le; j++) {
						if (ss[j] == 'N') 
							cntN ++;	       
					}
					free(ss);
					
					if ( (cntN / (float)le) < 0.1) {
						break;
					} 
				}				
			} // end while(1) loop
			
			sprintf(peakid, "%s-%d-%d-random", c, st, en);	
			
		} // end if (randomize == 1)
		else if (randomize == 2) {
			
			// not really random; pull out sequences from flanking region
			int ra = default_rand();
			if (ra < 0.5) {  // left side
				st = max(0, i - le);
				en = i;
			} else {
				st = j;
				en = min(chrlen-1, j + le);
			}
			
			sprintf(peakid, "%s\t%d\t%d", c, st, en);	
		}
		// randomize == 0
		else if (randomize == 0) {
			if (id != 0) {
				if(show_clust == 1) {
					sprintf(peakid, "%s-%d-%d-%s", c, st, en, id);
				}
			} else {
				
				if (show_peakdesc == 0) 
					sprintf(peakid, "%s-%d-%d", c, st, en);
				else 
					sprintf(peakid, "%s", line);
			}
			
		}
		// randomize == 3. When it comes from CGI mode, we want to add the "random" notation
		else if (randomize == 3) {
			if (id != 0) {
				if(show_clust == 1) {
					sprintf(peakid, "%s-%d-%d-%s-random", c, st, en, id);
				}
			} else {
				
				if (show_peakdesc == 0) 
					sprintf(peakid, "%s-%d-%d-random", c, st, en);
				else 
					sprintf(peakid, "%s", line);
			}
		}
		
		else {
			die("Randomize mode unknown\n");
		}
		
		//int   len = en - st + 1;
		char* ss =  my_fai_fetch(fidx, c, st, en);
		if (ss == 0) {
			printf("Problem: could not extract %s:%d-%d\n", c, st, en);
			exit(-1);
		}
		uctransform(ss);
		
		if (blankpeaks == 1) {
			for (i=actextleft; i<min(chrlen-1,actextleft+orile); i++) {
				ss[i] = 'N';
			}
		}
		if(show_seq_only)
			printf("%s\n", ss);
		else
			printf(">%s\n%s\n\n", peakid, ss);
		
		free(ss);
		
		free(p);
		free(line);
		
	} // while fgets    
	
	/* cleanup */	
	if (id != NULL)
		free(id);
	free(peakid);
	free(buff);
	fai_destroy(fidx);
	HASH_destroy(&clh);
	for (i=0; i<numchroms; i++)
		free(chrnames[i]);
	free(chrnames);
	free(chrlens);
	fclose(f);

	return 0;
}



char* my_fai_fetch(faidx_t* idx, char* c, int i, int j)
{
	int len;
	char id[1000];
	sprintf(id, "%s:%d-%d", c, i, j);
	return fai_fetch(idx, id, &len);
}



