#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "dataio.h"
#include "statistics.h"
#include "sequences.h";

typedef struct _GenInt {
	int i;			//start position
	int j;			//end position
	char* c;			//chromosome
	char* id;			//id
	float score;		//score
	int strand;
} GenInt;


int main(int argc, char** argv) {
	
	// general
	long i, j;
	
	// files to compare
	FILE* f1;
	FILE* f2;			
	
	// reading file
	char* buff;
	int   mynmax = 100000;
	char** p;
	int    m;
	
	char* intervals1 = 0;
	char* intervals2 = 0;
	
	int   verbose = 0;
	
	int   maxnumint = 10000000;
	
	// store intervals
	GenInt* a_int1;
	GenInt* a_int2;
	int     numint1 = 0;
	int     numint2 = 0;
	
	int     ext1   = 0;
	int     ext2   = 0;
	int     hasid1 = 0;
	int     hasid2 = 0;
	int     rand2  = 0;
	
	int     showscores	= 0;
	int     showprofile = 0;
	int     show_ov_int = 0;
	char*   featname	= 0;
	
	int*    a_ov_int = 0;
	
	int     showpos  = 0;
	int     iswig    = 0;
	int     minx;
	int     maxx;
	int     showunion	= 0;
	int     hasheader1	= 0;
	int     hasheader2	= 0;
	
	int		distance_min	= 5000; //minimum distance between intervals
	long	distance_estim	= 0;	//holds the computed distance
	
	// params
	// -chipfile
	// -inputfile
	// -fraglen (200?)
	// -chr
	// -chrlen
	// -wpois (10^6 ?) 
	// -wsmooth (1kb?)
	
	int showstrand = 0;
	
	if (argc < 2) {
		die("Usage: FindIsolatedPeaks -intervals1 FILE -ext1 INT -hasid1 INT -intervals2 FILE -ext2 INT -hasid2 INT -showprofile INT -distance_min INT \n");
	}
	
	if (exist_parameter(argc, argv, "-intervals1"))
		intervals1 = get_parameter(argc, argv, "-intervals1");
	
	if (exist_parameter(argc, argv, "-ext1"))
		ext1 = atoi(get_parameter(argc, argv, "-ext1"));
	
	if (exist_parameter(argc, argv, "-hasid1"))
		hasid1 = atoi(get_parameter(argc, argv, "-hasid1"));
	
	if (exist_parameter(argc, argv, "-hasheader1"))
		hasheader1 = atoi(get_parameter(argc, argv, "-hasheader1"));
	
	if (exist_parameter(argc, argv, "-ext2"))
		ext2 = atoi(get_parameter(argc, argv, "-ext2"));
	
	if (exist_parameter(argc, argv, "-hasid2"))
		hasid2 = atoi(get_parameter(argc, argv, "-hasid2"));
	
	if (exist_parameter(argc, argv, "-hasheader2"))
		hasheader2 = atoi(get_parameter(argc, argv, "-hasheader2"));
	
	if (exist_parameter(argc, argv, "-rand2"))
		rand2 = atoi(get_parameter(argc, argv, "-rand2"));
	
	if (exist_parameter(argc, argv, "-intervals2"))
		intervals2 = get_parameter(argc, argv, "-intervals2");
	
	if (exist_parameter(argc, argv, "-verbose"))
		verbose = atoi(get_parameter(argc, argv, "-verbose"));
	
	if (exist_parameter(argc, argv, "-showscores"))
		showscores = atoi(get_parameter(argc, argv, "-showscores"));
	
	if (exist_parameter(argc, argv, "-showprofile"))
		showprofile = atoi(get_parameter(argc, argv, "-showprofile"));
	
	if (exist_parameter(argc, argv, "-showpos"))
		showpos = atoi(get_parameter(argc, argv, "-showpos"));
	
	if (exist_parameter(argc, argv, "-showstrand"))
		showstrand = atoi(get_parameter(argc, argv, "-showstrand"));
	
	if (exist_parameter(argc, argv, "-show_ov_int"))
		show_ov_int = atoi(get_parameter(argc, argv, "-show_ov_int"));
	
	if (exist_parameter(argc, argv, "-featname"))
		featname = get_parameter(argc, argv, "-featname");
	
	if (exist_parameter(argc, argv, "-iswig"))
		iswig = atoi(get_parameter(argc, argv, "-iswig"));
	
	if (exist_parameter(argc, argv, "-showunion"))
		showunion = atoi(get_parameter(argc, argv, "-showunion"));
	
	if (exist_parameter(argc, argv, "-distance_min"))
		distance_min = atoi(get_parameter(argc, argv, "-distance_min"));
	
	default_set_seed(time(NULL));
	
	// alloc counter vectors
	a_int1 = (GenInt*)malloc(maxnumint * sizeof(GenInt));
	if (a_int1 == 0) {
		die("Problem allocating a_int1.\n");
	}
	a_int2 = (GenInt*)malloc(maxnumint * sizeof(GenInt));
	if (a_int2 == 0) {
		die("Problem allocating a_int2.\n");
	}
	
	buff  = (char*)malloc(mynmax * sizeof(char));  
	
	// read first set of intervals
	numint1 = 0;
	f1 = fopen(intervals1, "r");
	if (!f1) {
		die("Cannot open f1\n");
	}
	
	if ((iswig == 1) || (hasheader1 == 1))
		fgets(buff, mynmax, f1);
	
	while (!feof(f1)) {
		fgets(buff, mynmax, f1);
		if (feof(f1))
			break; 
		chomp(buff);
		split_line_delim(buff, "\t", &p, &m);
		
		int idx = 0;
		if (hasid1 == 1) {
			idx = 1;
			a_int1[numint1].id = strdup(p[0]);
		} else {
			a_int1[numint1].id = 0;      
		}
		a_int1[numint1].c = strdup(p[idx]);
		a_int1[numint1].i = atoi(p[idx+1]);
		a_int1[numint1].j = atoi(p[idx+2]);
		a_int1[numint1].strand = atoi(p[idx+3]);
		if (showscores == 1) {
			a_int1[numint1].score = atoi(p[idx+4]);
		}
		
		free(p);
		numint1++;
		
		if (numint1 == maxnumint) {
			die("Max number of intervals reached.\n");
		}
		
	}
	
	printf("\n");
	
	// read second set of intervals
	numint2 = 0;
	f2 = fopen(intervals2, "r");
	//fgets(buff, mynmax, f2); // skip first line
	
	if ((iswig == 1) || (hasheader2 == 1))
		fgets(buff, mynmax, f2);
	
	
	while (!feof(f2)) {
		fgets(buff, mynmax, f2);
		if (feof(f2))
			break; 
		chomp(buff);
		split_line_delim(buff, "\t", &p, &m);
		
		
		int idx = 0;
		if (hasid2 == 1) {
			idx = 1;
			a_int2[numint2].id = strdup(p[0]);
		} else {
			a_int2[numint2].id = 0;      
		}
		a_int2[numint2].c = strdup(p[idx]);
		
		int p1 = atoi(p[idx+1]);
		int p2 = atoi(p[idx+2]);
		
		if (rand2 == 1) {
			int li = p2 - p1 + 1;
			p1 = rand_chr_interval(p[idx], li);
			p2 = p1 + li;
		}
		
		a_int2[numint2].i = p1;
		a_int2[numint2].j = p2;
		free(p);
		numint2++;
		if (numint2 == maxnumint) {
			die("Max number of intervals reached.\n");
		}
	}
	
	if (verbose == 1) {
		printf("Read both files.\n");
	}
	
	if ((showprofile >= 1) && (featname != 0)) {
		printf("FEATURE\t%s\nFEATTYPES\t%d\n", featname, (showprofile==1?0:1));
	}
	
	if (showprofile == 1) 
		printf("GENE\tBOUND\n");
	
	a_ov_int = (int*)calloc(1000, sizeof(int));
	
	int numinter = 0;
	
	for (i=0; i<numint1; i++) {
		
		if (hasid1 == 1)
			printf("%s\t", a_int1[i].id);
		if (showprofile == 0) {
			printf("%s\t%d\t%d\t", a_int1[i].c, a_int1[i].i, a_int1[i].j);
		}
		
		if (showstrand == 1)
			printf("%d\t", a_int1[i].strand);
		if (showscores == 1) 
			printf("%f\t", a_int1[i].score);
		
		int cnt = 0;
		
		minx = a_int1[i].i;
		maxx = a_int1[i].j;
		
		for (j=0; j<numint2; j++) {
			
			distance_estim = sequencesDistance(a_int1[i].i-ext1, a_int1[i].j+ext1, a_int2[j].i-ext2, a_int2[j].j+ext2);
			
			if ((strcmp(a_int1[i].c, a_int2[j].c) == 0) && (distance_estim > distance_min)) {
				if (a_int2[j].i < minx) {
					minx = a_int2[j].i;
				}
				if (a_int2[j].j > maxx)
					maxx = a_int2[j].j;
				
				a_ov_int[cnt] = j;
				cnt++;
			}
		}
		
		int mycnt = cnt;
		if ((showprofile == 1) && (mycnt >= 1))
			mycnt = 1;
		
		printf("%d\n", mycnt);
		
		if (cnt > 0) 
			numinter ++;
	}
	
	return 0;
}
