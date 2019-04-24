#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dataio.h"
#include "statistics.h"
#include "sequences.h";

#include "protos.h"


#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

typedef struct _GenInt {
	int i;
	int j;
	char* c;
	char* id;
	float score;
} GenInt;
 




int main(int argc, char** argv) {
	
	// general
	long i;
	
	//FILE* fi;
	//FILE* fc;
	
	// reading file
	char* buff;
	int   mynmax = 100000;
	char** p;
	int    m;
	
	unsigned short*	chip_counts;
	unsigned short*	input_counts;
	
	char*				chipfile		= 0; //"SOLEXA/s_1_090126_hg18.eland.txt.chr17.fa";
	char*				inputfile		= 0; //"SOLEXA/s_8_090126_hg18.eland.txt.chr17.fa";  
	//char*				chr				= 0; //"chr17.fa";
	long				chrlen			= 0; //78774742; //247249719;
	
	int				fraglen         = 250;
	int				readlen         = 36;
	int				numreads_chip   = 0;
	int				numreads_input  = 0;
	
	int				t				= 5;
	//int				t_prev			= 0;
	//int				t_start			= 0;
	//int				t_end			= 0;
	
	//int				t_sumnumreads_chip  = 0;
	//int				t_sumnumreads_input = 0;
	
	//double			t_sum_lpr		= 0.0;
	//double			t_sum_lpi		= 0.0;
	//double			t_sum_lpc		= 0.0;
	
	double			fold_t			= 2.0;
	
	int				verbose			= 0;
	char*				chrname			= 0;
	double			mappability		= 1.0;
	int				fast			= 1;
	int				format			= 0;
	char*				formattext		= 0;
	// params
	// -chipfile
	// -inputfile
	// -fraglen (200?)
	// -chr
	// -chrlen
	// -wpois (10^6 ?) 
	// -wsmooth (1kb?)
	
	FILE*   f1;
	GenInt* a_int;
	int     maxnumint = 10000000;
	int     numint;
	int     hasid		= 0;
	char*   intervals = 0;
	int     ws		= 100;
	//float   eps		= 0.01;
	int     j;
	int	  k;
	char*                           chrdata         = 0;
	char**                          chrnames        = 0;
	int*                            chrlens         = 0;
	int numchroms = 24;

	if (argc < 2) {
		die("Usage: ChIPfinder -intervals FILE -chipfile FILE -inputfile FILE -t INT -fraglen INT -chrlen INT\n");
	}
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-hasid"))
		hasid = atoi(get_parameter(argc, argv, "-hasid"));
	
	if (exist_parameter(argc, argv, "-ws"))
		ws = atoi(get_parameter(argc, argv, "-ws"));
	
	if (exist_parameter(argc, argv, "-chipfile"))
		chipfile  = get_parameter(argc, argv, "-chipfile");
	
	if (exist_parameter(argc, argv, "-inputfile"))
		inputfile  = get_parameter(argc, argv, "-inputfile");
	
	if (exist_parameter(argc, argv, "-chrname"))
		chrname  = get_parameter(argc, argv, "-chrname");
	
	if (exist_parameter(argc, argv, "-t"))
		t  = atoi(get_parameter(argc, argv, "-t"));
	
	if (exist_parameter(argc, argv, "-fold_t"))
		fold_t  = atof(get_parameter(argc, argv, "-fold_t"));
	
	if (exist_parameter(argc, argv, "-fraglen"))
		fraglen  = atoi(get_parameter(argc, argv, "-fraglen"));
	
	if (exist_parameter(argc, argv, "-readlen"))
		readlen  = atoi(get_parameter(argc, argv, "-readlen"));
	
	if (exist_parameter(argc, argv, "-chrlen"))
		chrlen  = atoi(get_parameter(argc, argv, "-chrlen"));
	else {
		chrlen = hg18_chrlen(chrname);
		if (chrlen < 0) {
			printf("Problem ... chr %s unrecognized.\n", chrname);
		}
	}
	
	if (exist_parameter(argc, argv, "-mappability"))
		mappability = atof(get_parameter(argc, argv, "-mappability"));
	
	if (exist_parameter(argc, argv, "-verbose"))
		verbose = atoi(get_parameter(argc, argv, "-verbose"));
	
	if (exist_parameter(argc, argv, "-fast"))
		fast = atoi(get_parameter(argc, argv, "-fast"));
	
	if (exist_parameter(argc, argv, "-format")) {
		formattext = get_parameter(argc, argv, "-format");
		if (strcmp(formattext, "mit") == 0)
			format = 1;
		else if (strcmp(formattext, "bed") == 0)
			format = 2;
		else if (strcmp(formattext, "sam") == 0)
		  format = 3;
	}
	
	
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
		
		int idx = 0;
		if (hasid == 1) {
			idx = 1;
			a_int[numint].id = strdup(p[0]);
		} else {
			a_int[numint].id = 0;      
		}
		a_int[numint].c = strdup(p[idx]);
		a_int[numint].i = atoi(p[idx+1]);
		a_int[numint].j = atoi(p[idx+2]);
		a_int[numint].score = 0; 
		
		free(p);
		numint++;
		
		if (numint == maxnumint) {
			die("Max number of intervals reached.\n");
		}
		
	}
	//
	// end read intervals
	//
	
	// alloc counte vectors
	chip_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
	if (chip_counts == 0) {
		
		die("Problem allocating chip_counts.\n");
	}
	input_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
	if (input_counts == 0) {
		die("Problem allocating chip_counts.\n");
	}
	
	getCountFromReadFile(chipfile, format, chrname, chrlen, readlen, fraglen, &chip_counts, &numreads_chip);
	
	if (inputfile != 0)
		getCountFromReadFile(inputfile, format, chrname, chrlen, readlen, fraglen, &input_counts, &numreads_input);
	
	if (verbose == 1)
		printf("Number of reads = %d (chip) and %d (input)\n", numreads_chip, numreads_input);
	
	double scale_numreads = numreads_chip / (double)numreads_input;
	
	// now go thru all intervals
	float* tmp_rat;
	for (i=0; i<numint; i++) {
		
		if ( strcmp(a_int[i].c, chrname) != 0 )
			continue;
		
		// int length
		int li = a_int[i].j - a_int[i].i + 1;
		
		// allocate temporary interval
		tmp_rat   = (float*)calloc(li, sizeof(float));
		//tmp_rat_s = (float*)calloc(li, sizeof(float));
		
		// go thru interval
		for (j=a_int[i].i,k=0; j<=a_int[i].j; j++,k++) {
			//tmp_rat[k] = log( chip_counts[j] + eps) / ( (input_counts[j] + eps) * scale_numreads);
			
			if (inputfile != 0)
				tmp_rat[k] = max(0, chip_counts[j] - input_counts[j] * scale_numreads) ;
			else 
				tmp_rat[k] = chip_counts[j] ;
			
			//if (tmp_rat[k] < 0) 
			//	tmp_rat[k] = 0;
			
			//tmp_rat[k] = log( chip_counts[j] - input_counts[j] + eps)  scale_numreads);
			
			// printf(" %f (%d / %d*%f)\n", tmp_rat[k],  chip_counts[j],  input_counts[j], scale_numreads);
		}
		
		// smooth intervals
		float maxpeak = -100000.0;
		
		for (j=0; j<li-ws; j++) {
			float sum = 0.0;
			for (k=j; k<j+ws; k++) {
				sum += tmp_rat[k]; 
			}
			sum /= ws;
			
			if (sum > maxpeak) {
				maxpeak = sum;
			}
		}
		
		
		printf("%s\t", a_int[i].id);
		printf("%4.3f\n", maxpeak);
		
	}
	
	return 0;
	
}

