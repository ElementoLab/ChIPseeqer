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
	int     i			= 0;
	int		k			= 0;
	char*   intervals	= NULL;
	char*   buff		= NULL;
	FILE*   f			= NULL;
	int     mynmax		= 100000;
	char**  p			= NULL;
	int     m			= 0;
	int     hasid		= 0;
	int     randomize	= 0;
	int     minlen		= -1;
	int     forcelen	= -1;
	int     extend		= -1;
	int     iswig		= 0;
	char*   chr			= 0;
	int     seed		= 2345;
	char*   h			= 0;
	int     header		= 0;  
	int     cidx		= 0;
	int		distance	= 10000;
	HASH    clh;
	int     verbose         = 0;
	char*   genome          = "hg18"; 
	char*   chrdata         = 0;
	int     numchroms       = 0;
	int*    chrlens         = 0;
	char**  chrnames        = 0;
	float*  chrmaps         = 0;
	int     chrlen          = 0;
	int     idx             = 0;
	int     actextleft      = -1;
	int     actextright     = -1;
	int     newen           = -1;
	int		newst           = -1;
	int		length			= 200;
	int		rand_en			= -1;
	int		rand_st			= -1;
	
	if (argc == 1) {
		die("Args: -peakfile FILE -fastafile FILE [ -extend INT -forcelen INT -minlen INT ] \n");
	}
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals  = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-hasid"))
		hasid  = atoi(get_parameter(argc, argv, "-hasid"));
	
	if (exist_parameter(argc, argv, "-header"))
		header  = atoi(get_parameter(argc, argv, "-header"));
	
	if (exist_parameter(argc, argv, "-seed"))
		seed  = atoi(get_parameter(argc, argv, "-seed"));
	
	if (exist_parameter(argc, argv, "-randomize"))
		randomize  = atoi(get_parameter(argc, argv, "-randomize"));
	
	if (exist_parameter(argc, argv, "-minlen"))
		minlen = atoi(get_parameter(argc, argv, "-minlen"));
	
	if (exist_parameter(argc, argv, "-forcelen"))
		forcelen = atoi(get_parameter(argc, argv, "-forcelen"));
	
	if (exist_parameter(argc, argv, "-extend"))
		extend = atoi(get_parameter(argc, argv, "-extend"));
	
	if (exist_parameter(argc, argv, "-iswig"))
		iswig = atoi(get_parameter(argc, argv, "-iswig"));
	
	if (exist_parameter(argc, argv, "-distance"))
		distance = atoi(get_parameter(argc, argv, "-distance"));
	
	if (exist_parameter(argc, argv, "-genome"))
		genome = get_parameter(argc, argv, "-genome");
	
	if (exist_parameter(argc, argv, "-verbose"))
		verbose = atoi(get_parameter(argc, argv, "-verbose"));
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmaps, &numchroms);
	} else {
		get_hg18_chroms(&chrnames, &numchroms, &chrlens);
	}
	
	if (exist_parameter(argc, argv, "-length"))
		length = atoi(get_parameter(argc, argv, "-length"));
	
	// RNG init
	if (seed == 0) {
		int t = (unsigned)(getpid() + time(NULL));
		default_set_seed(t);
		//printf("ran1_seed:%d\t", t);
	}
	else {
		default_set_seed(seed);
	}
		
	// store chromsomal indices  
	HASH_init(&clh, numchroms);
	for (i=0; i<numchroms; i++) {
		HASH_insert(&clh, chrnames[i], i);
	}
	
	char* id		= (char*)calloc(10000, sizeof(char));
	char* peakid	= (char*)calloc(10000, sizeof(char));
	buff			= (char*)malloc(mynmax * sizeof(char));
	
	//
	// load up intervals file
	//
	f = fopen(intervals, "r");
	if (f == 0)
		die("cannot open interval file.\n");    
	if ((iswig == 1) || (header == 1)) {
		fgets(buff, mynmax, f);
		h = strdup(buff);
		chomp(h);    
	}
	
	//for each line
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
			//if (hasid == 1)
			//	free(id);
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
			//if (hasid == 1)
			//	free(id);      
			continue;
		}
		
		// update length
		length = abs(en - st + 1);
		
		// generate randomly located peak
		if (randomize == 1) {	 
						
			//rand_st = rand_chr_interval2(c, length, chrlen);	

			rand_st = (int)(0.5 + default_rand() * chrlen - length);
			
			if((rand_st + length) > chrlen) {
				rand_en = chrlen;
			}
			else {
				rand_en = rand_st + length;
			}
						
			sprintf(peakid, "%s\t%d\t%d", c, rand_st, rand_en);	
			
		} // end if (randomize == 1)
		else if (randomize == 2) {
			
			// not really random; pull out sequences from flanking regions
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
		else {
			die("Randomize mode unknown\n");
		}
		
		printf("%s\n", peakid);
		
		free(p);
		free(line);
		
	} // while fgets    
	
	/* cleanup */
	fclose(f);
	free(buff);
	free(peakid);
	free(h);
	for (k=0; k<numchroms; k++)
		free(chrnames[k]);
	free(chrnames);
	free(chrlens);
	free(chrmaps);
	if(id != NULL)
		free(id);
	HASH_destroy(&clh);
	
	return 0;
}

