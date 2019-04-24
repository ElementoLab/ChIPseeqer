#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <search.h>

#include "dataio.h"
#include "prefix.h"
#include "sequences.h"
#include "statistics.h"
#include "hashtable.h"

#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)


typedef struct _GenInt {
	int i;
	int j;
	char* c;
	char* id;
	char* str;
} GenInt;


int main(int argc, char ** argv)
{
	//int      k = 0; 
	char*    fastafile = 0; 
	
	seqI     si;
	char*    seq;
	int      size;
	char*    name;
	
	int      i, l, j;
	//char*    kmer;
	//char*    cmer;
	//int      nb1, nb2;
	//char**   kmers;
	//char**   newkmers;
	//int      nbkmers		= 20000000; // estimated
	//int      nbkmers_inc	= 100000;
	//int*     counts;
	//int      index_kmers	= 0;
	int      cnt_seq		= 0;
	
	char*    intervals		= 0;
	
	// store intervals
	GenInt* a_int;
	char*   buff;
	int     numint			= 0;
	FILE*   f;
	//int		maxnumint		= 10000000;
	int		mynmax			= 100000;
	char**	p;
	int		m;
	int		hasid			= 0;
	int		randomize		= 0;
	int		minlen			= -1;
	int		forcelen		= -1;
	int		extend			= -1;
	int		iswig			= 0;
	int		rand_close		= 0;
	int		distance		= 10000;
	int		show_clust		= 0;
	int		print_both		= 0;
	char*   id              = 0;
	int     show_peakdesc   = 0;
	char*   randintfile     = 0;
	HASH    hc;
	int     hasrandid       = 0;
	int     chridx          = -1;
	int     chrnum          = 0;
	int     verbose         = 0;
	int     genome          = 0; //0 for human, 1 otherwise
	int		defaultSeed		= 0;
	
	if (argc == 1) {
		die("Args: -fastafile FILE -intervals FILE [ -extend INT -forcelen INT -minlen INT ] \n");
	}
	
	
	if (exist_parameter(argc, argv, "-fastafile"))
		fastafile  = get_parameter(argc, argv, "-fastafile");
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals  = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-peakfile"))
		intervals  = get_parameter(argc, argv, "-peakfile");
	
	if (exist_parameter(argc, argv, "-hasid"))
		hasid  = atoi(get_parameter(argc, argv, "-hasid"));
	
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
	
	if (exist_parameter(argc, argv, "-rand_close"))
		rand_close = atoi(get_parameter(argc, argv, "-rand_close"));
	
	// intervals for random search
	if (exist_parameter(argc, argv, "-randintfile"))
		randintfile = get_parameter(argc, argv, "-randintfile");
	if (exist_parameter(argc, argv, "-hasrandid"))
		hasrandid  = atoi(get_parameter(argc, argv, "-hasrandid"));
	
	
	if (exist_parameter(argc, argv, "-show_clust"))
		show_clust = atoi(get_parameter(argc, argv, "-show_clust"));
	
	if (exist_parameter(argc, argv, "-print_both"))
		print_both = atoi(get_parameter(argc, argv, "-print_both"));
	
	if (exist_parameter(argc, argv, "-show_peakdesc"))
		show_peakdesc = atoi(get_parameter(argc, argv, "-show_peakdesc"));
	
	if (exist_parameter(argc, argv, "-genome"))
		genome = atoi(get_parameter(argc, argv, "-genome"));

	if (exist_parameter(argc, argv, "-defaultSeed"))
		defaultSeed = atoi(get_parameter(argc, argv, "-defaultSeed"));
	
	if (defaultSeed == 0) {
		default_set_seed(time(NULL));
	}
	else {
		default_set_seed(1234);
	}
	
	// read intervals
	int maxnumint = nbLinesInFile(intervals);
	
	a_int = (GenInt*)malloc(maxnumint * sizeof(GenInt));
	if (a_int == 0) {
		die("Problem allocating a_int.\n");
	}
	
	
	buff  = (char*)malloc(mynmax * sizeof(char));
	
	// read first set of intervals
	numint = 0;
	f = fopen(intervals, "r");
	if (f == 0)
		die("cannot open interval file.\n");
	
	// hash to store chr num
	HASH_init(&hc, 100);
	
	if (iswig == 1)
		fgets(buff, mynmax, f);
	
	id = (char*)calloc(1000, sizeof(char));
	
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
		
		
		
		int idx = 0;
		
		if (hasid == 1) {
			idx = 1;
			a_int[numint].id = strdup(p[0]);
		} else {
			a_int[numint].id = 0;      
		}
		a_int[numint].c = strdup(p[idx]);
		
		//HASH_enter(&hc, p[idx], chrnum);
		
		
		if (!HASH_find(&hc, p[idx], &chridx)) {
			HASH_enter(&hc, p[idx], chrnum);
			//printf("Enter new chr id = %d\n", chrnum);
			chrnum++;
		} else {
			//printf("%s already found\n", p[idx]);
		}
		
		a_int[numint].i = atoi(p[idx+1]);
		a_int[numint].j = atoi(p[idx+2]);
		a_int[numint].str = line;
		if (a_int[numint].j < a_int[numint].i) {
			printf("Problem reading intervals ... j (%d) < i (%d) \n", a_int[numint].j, a_int[numint].i);
			die("Problem reading intervals ... j < i \n");
		}
		
		free(p);
		numint++;
		
		//if (numint == maxnumint) {
		//	die("Max number of intervals reached.\n");
		//}
		
	}
	fclose(f);
	
	
	if (verbose == 1)
		printf("Found %d peaks\n", numint);
	
	// load intervals for random peaks
	if (randintfile != 0) {
		
		// read first set of intervals
		//int numrandint = 0;
		f = fopen(randintfile, "r");
		if (f == 0)
			die("cannot open random interval file.\n");
		
		while (fgets(buff, mynmax, f) != 0) {
			chomp(buff);
			if (buff[0] == '#')
				continue;
			split_line_delim(buff, "\t", &p, &m);
			
			int idx = 0;
			if (hasrandid == 1) {
				idx = 1;
				a_int[numint].id = strdup(p[0]);
			} else {
				a_int[numint].id = 0;      
			}
			a_int[numint].c = strdup(p[idx]);
			a_int[numint].i = atoi(p[idx+1]);
			a_int[numint].j = atoi(p[idx+2]);
			//a_int[numint].str = line;
			if (a_int[numint].j < a_int[numint].i) {
				die("Problem reading intervals ... j < i\n");
			}
			
			free(p);
			numint++;
			
			if (numint == maxnumint) {
				die("Max number of intervals reached.\n");
			}
			
		}
		
	} // if (randintfile)
	
	
	
	// working for the human genome (fast)
	si.verbose = 0;
	seqI_set_max_seqlen(&si, 500000000UL);
	seqI_set_seqlen_inc(&si, 300000000UL);
	seqI_set_fread_chunk(&si, 100000000UL);
	
	if (seqI_open_fast(&si, fastafile) == 0) {
		die("Error opening file\n");
	}
	
	
	cnt_seq = 1;
	
	char* modname = (char*)calloc(1000, sizeof(char));
	
	while ((seq = seqI_nextSequence_fast(&si, &name, &size))) {
		//printf("SEQ: %s; l=%d, %d\n", name, (int)(strlen(seq)), cnt_seq);

		int m = strlen(name);

		for (i=0; i<m; i++) {
			if (name[i] == ' ') {
 				name[i] = '\0';
				break;
			}
		}
		
		l = strlen(seq);
		
		if (genome == 1) {
			sprintf(modname, "chr%s", name);
			name = modname;
		}
		
		
		for (i=0; i<numint; i++) {
			if (strcmp(a_int[i].c, name) == 0) {
				
				//id        = 0;
				int st    = a_int[i].i;
				int en	  = a_int[i].j;
				int le    = abs(en - st + 1);
				
				if ((forcelen != -1) && (le < forcelen)) {
					int sideext = (int)( 0.5 + ( forcelen - le ) / 2 );
					st = max(0,st - sideext);
					en = min(en + sideext,l-1);
				}
				
				if (extend != -1) {
					st = max(0,st - extend);
					en = min(en + extend,l-1);
				}
				
				if ((minlen >= 0) && (le < minlen)) {
					continue;
				}
				
				// update length
				le = abs(en - st + 1);
				
				
				// generate randomly located peak
				// if rand_close == 1, random peak is selected from local distance .... 
				if (randomize == 1) {	 
					
					while (1) {
						
						st = rand_chr_interval2(a_int[i].c, le, (int)(strlen(seq)));
							
						if (rand_close) {
							if(((st > a_int[i].j) && (st < a_int[i].j + distance)) || ((st < a_int[i].i) && (st > a_int[i].i-distance))){
								en = st + le;
								char* ss =  uc(substr(seq, st, abs(en - st + 1)));
								
								int cntN = 0;
								for (j=0; j<le; j++) {
									if (ss[j] == 'N') 
										cntN ++;	       
								}
								
								if ( (cntN / (float)le) < 0.1) {
									break;
								} else {
									free(ss);
								}
							}
						}
						else {
							en = st + le;
							char* ss =  uc(substr(seq, st, abs(en - st + 1)));
							
							int cntN = 0;
							for (j=0; j<le; j++) {
								if (ss[j] == 'N') 
									cntN ++;	       
							}
							
							if ( (cntN / (float)le) < 0.1) {
								break;
							} else {
								free(ss);
							}
						}
						
						
					}
					
					//id = (char*)calloc(1000, sizeof(char));
					
					sprintf(id, "%s-%d-%d-random", a_int[i].c, st, en);	
					
					
					
				} // end if (randomize == 1)
				else if (randomize == 2) {
					
					// not really random; pull out sequences from flanking region
					int ra = default_rand();
					if (ra < 0.5) {  // left side
						st = max(0, a_int[i].i - le);
						en = a_int[i].i;
					} else {
						st = a_int[i].j;
						en = min(l-1, a_int[i].j + le);
					}
					sprintf(id, "%s-%d-%d-random", a_int[i].c, st, en);	
					
				}
				// randomize == 0
				else if (randomize == 0) {
					if (a_int[i].id != 0) {
						if(show_clust == 1) {
							// id = (char*)calloc(1000, sizeof(char));
							sprintf(id, "%s-%d-%d-%s", a_int[i].c, st, en, a_int[i].id);
						}
						else {
							id  = a_int[i].id;
						}
					} else {
						//id = (char*)calloc(1000, sizeof(char));
						if (show_peakdesc == 0) 
							sprintf(id, "%s-%d-%d", a_int[i].c, st, en);
						else 
							sprintf(id, "%s", a_int[i].str);
					}
					
				} // end else (randomize == 0)
				
				else {
					die("Randomize mode unknown\n");
				}
				
				
				// matters only when we want to print intervals and randintervals
				if (print_both == 1) {
					printf("%s-%d-%d\t%s\n", a_int[i].c, a_int[i].i, a_int[i].j, id);
				}
				else {
					
					char* ss = substr(seq, st, abs(en - st + 1));
					char* su = uc(ss);
					printf(">%s\n%s\n\n", id, su);
					free(su);
					free(ss);
					
				}
				//if (a_int[i].id != 0)
				//	free(id);
				
			} // if right chr
		} // loop over intervals
		
		cnt_seq++;
		free(seq);
		
		//free(name);
		
	} // while (seq == ...
	
	seqI_close(&si);
	
	return 0;
}