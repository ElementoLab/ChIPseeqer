#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <search.h>
#include "dataio.h"
#include "statistics.h"
#include "hashtable.h"
#include "sequences.h"

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

typedef struct _GenInt {
	int i;
	int j;
	char* c;
	char* id;
	float score;
	int del;
} GenInt;

int CmpGenIntStart(const void* _a, const void* _b); 

int main(int argc, char** argv) {
	
	// general
	long i;
	long k;
	
	// reading file
	unsigned short* chip_counts;
	char* chipfile			= NULL;
	long  chrlen			= 0; 
	int   fraglen			= 0;
	int   readlen			= -1; 
	int   numreads_chip		= 0;
	int   verbose			= 0;
	char* chrname			= NULL;
	int    format			= 0;
	char*  formattext		= NULL;
	int     ws				= 10;
	long    j				= 0;
	int     from            = -1;
	int     to              = -1;
	char*   desc			= NULL;
	char*   intervals		= NULL;
	FILE*   f1				= NULL;
	GenInt* a_int			= NULL;
	int     maxnumint		= 10000000;
	int     numint			= 0;
	int     hasid			= 1;
	
	// reading file
	char*	buff			= NULL;
	int		mynmax			= 100000;
	char**	p				= NULL;
	int		m				= 0;
	char*	tmpid			= NULL;
	int		idx				= 0;
	int		numchroms		= 0;
	int*	chrlens			= NULL;
	char**	chrnames		= NULL;
	int		c				= 0;
	char*	readdir			= NULL;
	char	readfile[1000];
	char*	chr				= NULL;
	int		bigwig			= 0;
	int		normalize		= 1;
	int		totalnumreads	= 0;
	int		numreads		= 0;
	double	normfactor		= 0;
	int		normto			= 10000000;
	char*	chrdata			= NULL;
	int		uniquereads		= 1;
	int		inc				= 10;
	int		numreadsclonal_chip = 0;
	if (argc < 2) {
		die("Usage: MakeGenomicReadDensityTrack -readdir DIR -ws INT -genome STR [ -normalize INT -bigwig INT -intervals FILE]\n");
	}
	
	if (exist_parameter(argc, argv, "-readdir"))
		readdir = get_parameter(argc, argv, "-readdir");
	
	if (exist_parameter(argc, argv, "-ws"))
		ws = atoi(get_parameter(argc, argv, "-ws"));
	
	if (exist_parameter(argc, argv, "-inc"))
		inc = atoi(get_parameter(argc, argv, "-inc"));
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-chipfile"))
		chipfile  = get_parameter(argc, argv, "-chipfile");
	
	if (exist_parameter(argc, argv, "-chrname"))
		chrname  = get_parameter(argc, argv, "-chrname");
	
	if (exist_parameter(argc, argv, "-fraglen"))
		fraglen  = atoi(get_parameter(argc, argv, "-fraglen"));
	
	if (exist_parameter(argc, argv, "-uniquereads"))
		uniquereads  = atoi(get_parameter(argc, argv, "-uniquereads"));
	
	if (exist_parameter(argc, argv, "-chr"))
		chr  = get_parameter(argc, argv, "-chr");
	//else 
	//chrlen  = hg18_chrlen(chrname);
	
	if (exist_parameter(argc, argv, "-verbose"))
		verbose = atoi(get_parameter(argc, argv, "-verbose"));
	
	if (exist_parameter(argc, argv, "-desc"))
		desc = get_parameter(argc, argv, "-desc");
	
	if (exist_parameter(argc, argv, "-bigwig"))
		bigwig = atoi(get_parameter(argc, argv, "-bigwig"));
	
	if (exist_parameter(argc, argv, "-format")) {
		formattext = get_parameter(argc, argv, "-format");
		if (strcmp(formattext, "mit") == 0)
			format = 1;
		else if (strcmp(formattext, "bed") == 0)
			format = 2;
		else if (strcmp(formattext, "sam") == 0)
			format = 3;    
	}
	
	if (exist_parameter(argc, argv, "-genome")) {
		
	} else if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, 0, &numchroms);
		//genome   = 0;
	} else {
		get_hg18_chroms(&chrnames, &numchroms, &chrlens);
	}
	
	if (exist_parameter(argc, argv, "-normalize"))
		normalize = atoi(get_parameter(argc, argv, "-normalize"));
	
	if (exist_parameter(argc, argv, "-from"))
		from = atoi(get_parameter(argc, argv, "-from"));
	
	if (exist_parameter(argc, argv, "-to"))
		to = atoi(get_parameter(argc, argv, "-to"));
	
	if (exist_parameter(argc, argv, "-chr"))
		chr = get_parameter(argc, argv, "-chr");
	
	if (intervals != 0) {
		
		//
		// start reading intervals
		//
		
		a_int = (GenInt*)malloc(maxnumint * sizeof(GenInt));
		if (a_int == 0) {
			die("Problem allocating a_int.\n");
		}
		
		buff  = (char*)malloc(mynmax * sizeof(char));
		
		// read first set of intervals
		numint = 0;
		f1 = fopen(intervals, "r");
		if (!f1) {
			die("Cannot open f1\n");
		}
		while (!feof(f1)) {
			fgets(buff, mynmax, f1);
			if (feof(f1))
				break; 
			chomp(buff);
			split_line_delim(buff, "\t", &p, &m);
			
			idx   = 0;
			tmpid = 0;
			if (hasid == 1) {
				idx = 1;
				tmpid = strdup(p[0]);
			} else {
				idx = 0;
				tmpid = 0;
			}
			
			if (strcmp(p[idx], chrname) == 0) {
				a_int[numint].id = tmpid;      
				a_int[numint].c  = strdup(p[idx]);
				a_int[numint].i  = atoi(p[idx+1]);
				a_int[numint].j  = atoi(p[idx+2]);
				a_int[numint].score = 0; 
				a_int[numint].del = 0;
				numint++;
				
				if (numint == maxnumint) {
					die("Max number of intervals reached.\n");
				}
			}      
			
			free(p);
		}
		
		//
		// end read intervals
		//
		
		// sort intervals
		qsort((void*)a_int, numint, sizeof(GenInt), CmpGenIntStart);
		
		
		// merge overlapping intervals
		for (i=0; i<numint-1; i++) {
			
			// int i      
			if (a_int[i].del == 1)
				continue;
			
			for (j=i+1; j<numint; j++) {
				
				// int j      
				if (a_int[j].del == 1)
					continue;
				
				if (sequencesOverlap(a_int[i].i, a_int[i].j, a_int[j].i, a_int[j].j)) {
					a_int[i].j   = max(a_int[i].j, a_int[j].j);
					a_int[j].del = 1;
				}				
			}
		}
	}
	
	
	// calc total number of reads for RPKM-style 
	if (normalize == 1) {
		totalnumreads = 0;
		for (c=0; c<numchroms; c++) {
			chrname = chrnames[c];
			sprintf(readfile, "%s/reads.%s", readdir, chrname); 
			numreads = CountReads(readfile, format);
			if (verbose == 1) {
				fprintf(stdout, "Found %d aligned reads in chr %s\n", numreads, chrname);
			}
			totalnumreads += numreads;
			
		}
		if (verbose == 1) {
			fprintf(stdout, "Found a total of %d aligned reads\n", totalnumreads);
		}
		
		normfactor = normto / (double)totalnumreads;
		if (verbose == 1) {
			fprintf(stdout, "Normalization factor = %4.3f\n", normfactor);
		}
		
	}
	
	
	if ((bigwig == 0) || (desc != 0)) {
		printf("track type=wiggle_0 name=\"%s\" description=\"%s\" visibility=\"full\" maxHeightPixels=\"64:64:11\" smoothingWindow=\"10\" viewLimits=\"0:100\" autoScale=\"off\"\n", desc, desc);
	}
	
	for (c=0; c<numchroms; c++) {
		
		chrname = chrnames[c];
		chrlen  = chrlens[c];
		
		sprintf(readfile, "%s/reads.%s", readdir, chrname); 
		
		// alloc counter vectors
		chip_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
		if (chip_counts == 0) {
			die("Problem allocating chip_counts.\n");
		}
		
		getCountFromReadFile(readfile,  format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &chip_counts, &numreads_chip, &numreadsclonal_chip);
		
		//getCountFromReadFile(inputfile, format, chrname, chrlen, readlen, fraglen, &input_counts, &numreads_input);
		
		if (verbose == 1)
			fprintf(stdout, "Number of reads from chr %s, len=%ld  = %d\n", chrname, chrlen, numreads_chip);
		
		
		printf("variableStep chrom=%s span=%d\n", chrname, inc);
		
		
		double avg     = 0;
		int    mid     = 0;
		int    halfinc = inc/2;
		int    halfws  = ws/2;
		
		if (intervals == 0) {
			
			for (i=0; i<chrlen-inc; i+=inc) {    
				avg = 0;
				
				mid = i + halfinc;
				for (j=max(0,mid-halfws); j<min(chrlen,mid+halfws); j++) 
					avg += chip_counts[j];
				avg /= ws;
				
				if (normalize == 1)
					avg = avg * normfactor;
				
				if (avg > 0) {
					printf("%ld\t%3.1f\n", i+1, avg);
				}
			}
		} else {
			
			for (i=0; i<numint; i++) {
				
				// skip if wrong chr
				if (a_int[i].del == 1)
					continue;
				
				// go thru interval, ws nt at a time
				for (j=a_int[i].i; j<=a_int[i].j; j+=inc) {
					mid = j + halfinc;
					avg = 0;
					for (k=max(0,mid-halfws); k<min(chrlen,mid+halfws); k++) 
						avg += chip_counts[k];
					avg /= ws;
					if (avg >= 1) {
						printf("%ld\t%3.1f\n", j+1, avg);
					}
				}
				
			} // end for     
			
		}
		free(chip_counts);
	}
	
	/* cleanup */
	for (c=0; c<numchroms; c++)
		free(chrnames[c]);
	free(chrnames);
	free(chrlens);
	
	return 0;
}


int CmpGenIntStart(const void* _a, const void* _b) 
{
	
	const GenInt* a = (const GenInt*) _a;
	const GenInt* b = (const GenInt*) _b;
	
	if (a->i > b->i) 
		return 1;
	else if (a->i == b->i) {
		return 0;
	} else {
		return -1;
	}
}

