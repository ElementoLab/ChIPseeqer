//
// this program returns the average or max read count for a list of interval (or motifs)
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <search.h>
#include "hashtable.h"
#include "dataio.h"
#include "statistics.h"
#include "sequences.h"


#define NUMCHROMS 25

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

typedef struct _GenInt {
	int   i;
	int   j;
	int   ci;
	int   cj;
	char* chr;
	char* tsname;
	char* genename;
	int   fr; // frame
	int   ef; // exon frame
	int   num; // exon number
	char  utrtype;
	
	int   mut[4];
	int   N;
	int   k;
	float score;
	char* dbsnp;
	float entropy;
	char* str;
} GenInt;

void readFileSumUpColumnPerChrom(char* file, int col, int colchr, int** counts, int** countgenes);
void readRefGenePromoters(char* file, GenInt** a_proms, int* numproms, int up, int dn);


int main(int argc, char** argv) {
	
	// general
	long i;
	long j;
	
	unsigned short*	chip_counts;
	
	char*			chipfile		= 0; //"SOLEXA/s_1_090126_hg18.eland.txt.chr17.fa";
	char*			inputfile		= 0; //"SOLEXA/s_8_090126_hg18.eland.txt.chr17.fa";  
	long			chrlen			= 0; //78774742; //247249719
	int				fraglen         = -1;
	int				readlen         = 0;
	int				numreads_chip   = 0;
	int				t				= 5;
	double			fold_t			= 2.0;
	
	int				verbose			= 0;
	char*			chrname			= 0;
	double			mappability		= 1.0;
	int				fast			= 1;
	int				format			= 0;
	char*			formattext		= 0;
	
	int				hasid			= 0;
	char*			intervals		= 0;
	int				ws				= 100;
	char*			chrdata         = 0;
	char**			chrnames        = 0;
	int*			chrlens         = 0;
	int				numchroms		= 24;
	char*			chipdir			= 0;
	
	char	readfile[1000];
	int		normalize		= 1;
	double	normfactor		= 0;
	int		numreads		= 0;
	int		totalnumreads	= 0;
	int		normto			= 10000000;
	int		c;
	int		genome			= 0;
	char*	inputdir		= 0;
	int		ext				= 0;
	char*	output_text		= 0;
	int		output			= 0;
	float	tmpval			= 0.0;
	int		motifmatches	= 0;
	int		randomize		= 0;
	int		seed			= 1000;
	float*	chrmap			= 0;
	
	char*	fastadir		= 0;
	int		checkN			= 0;
	int		uniquereads		= 1;
	int		refgene			= 0;
	int		up				= 2000;
	int		dn				= 2000;
	long*   dist			= 0;
	char*	mapdir			= 0;
	unsigned char* mapmap;
	char	mapfile[1000];
	int		mapt			= -1;
	int		numreadsclonal_chip = 0;
	
	if (argc < 2) {
		die("Usage: ChIPreadCountDistribution [ -refgene FILE -up INT -dn INT | -intervals FILE ] -chipdir DIR -inputdir FILE -t INT -fraglen INT -chrlen INT -output max|avg -chrdata FILE\n");
	}
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-motifmatches")) {
		intervals = get_parameter(argc, argv, "-motifmatches");
		motifmatches = 1;
	}
	
	// refgene
	if (exist_parameter(argc, argv, "-refgene")) {
		intervals = get_parameter(argc, argv, "-refgene");
		refgene = 1;
	}
	if (exist_parameter(argc, argv, "-dn"))
		dn = atoi(get_parameter(argc, argv, "-dn"));  
	if (exist_parameter(argc, argv, "-up"))
		up = atoi(get_parameter(argc, argv, "-up"));
	
	if (exist_parameter(argc, argv, "-hasid"))
		hasid = atoi(get_parameter(argc, argv, "-hasid"));
	
	if (exist_parameter(argc, argv, "-ws"))
		ws = atoi(get_parameter(argc, argv, "-ws"));
	
	if (exist_parameter(argc, argv, "-chipfile"))
		chipfile  = get_parameter(argc, argv, "-chipfile");
	
	if (exist_parameter(argc, argv, "-chipdir"))
		chipdir  = get_parameter(argc, argv, "-chipdir");
	
	if (exist_parameter(argc, argv, "-inputdir"))
		inputdir  = get_parameter(argc, argv, "-inputdir");
	
	if (exist_parameter(argc, argv, "-mapdir"))
		mapdir  = get_parameter(argc, argv, "-mapdir");  
	
	if (exist_parameter(argc, argv, "-mapt"))
		mapt  = atoi(get_parameter(argc, argv, "-mapt"));  
	
	if (exist_parameter(argc, argv, "-fastadir"))
		fastadir  = get_parameter(argc, argv, "-fastadir");
	
	if (exist_parameter(argc, argv, "-checkN"))
		checkN  = atoi(get_parameter(argc, argv, "-checkN"));
	
	if (exist_parameter(argc, argv, "-inputfile"))
		inputfile  = get_parameter(argc, argv, "-inputfile");
	
	if (exist_parameter(argc, argv, "-chrname"))
		chrname  = get_parameter(argc, argv, "-chrname");
	
	if (exist_parameter(argc, argv, "-uniquereads"))
		uniquereads  = atoi(get_parameter(argc, argv, "-uniquereads"));
	
	if (exist_parameter(argc, argv, "-t"))
		t  = atoi(get_parameter(argc, argv, "-t"));
	
	if (exist_parameter(argc, argv, "-normalize"))
		normalize  = atoi(get_parameter(argc, argv, "-normalize"));
	
	if (exist_parameter(argc, argv, "-normto"))
		normto  = atoi(get_parameter(argc, argv, "-normto"));  
	
	if (exist_parameter(argc, argv, "-randomize"))
		randomize  = atoi(get_parameter(argc, argv, "-randomize"));
	
	if (exist_parameter(argc, argv, "-fold_t"))
		fold_t  = atof(get_parameter(argc, argv, "-fold_t"));
	
	if (exist_parameter(argc, argv, "-fraglen"))
		fraglen  = atoi(get_parameter(argc, argv, "-fraglen"));
	
	if (exist_parameter(argc, argv, "-readlen"))
		readlen  = atoi(get_parameter(argc, argv, "-readlen"));
	
	if (exist_parameter(argc, argv, "-normto"))
		normto  = atoi(get_parameter(argc, argv, "-normto"));
	
	if (exist_parameter(argc, argv, "-mappability"))
		mappability = atof(get_parameter(argc, argv, "-mappability"));
	
	if (exist_parameter(argc, argv, "-ext"))
		ext = atoi(get_parameter(argc, argv, "-ext"));
	
	if (exist_parameter(argc, argv, "-output")) {
		output_text = get_parameter(argc, argv, "-output");
		if (strcmp(output_text, "avg") == 0)
			output = 0;
		else if (strcmp(output_text, "max") == 0)
			output = 1;	  
	}
	
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
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata  = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);
		genome   = 0;
	}
	
	if (exist_parameter(argc, argv, "-seed")) {
		seed = atoi(get_parameter(argc, argv, "-seed"));  
	}
	default_set_seed(seed);
    
	
	
	//
	// end read intervals
	//
	
	// calc total number of reads for RPKM-style 
	if (normalize == 1) {
		totalnumreads = 0;
		for (c=0; c<numchroms; c++) {
			chrname = chrnames[c];
			sprintf(readfile, "%s/reads.%s", chipdir, chrname); 
			
			if (!file_exists(readfile))
				continue;
			
			numreads = CountReads(readfile, format);
			if (verbose == 1) {
				fprintf(stderr, "Found %d aligned reads in chr %s\n", numreads, chrname);
			}
			totalnumreads += numreads;
			
		}
		if (verbose == 1) {
			fprintf(stderr, "Found a total of %d aligned reads\n", totalnumreads);
		}
		
		normfactor = normto / (double)totalnumreads;
		if (verbose == 1) {
			fprintf(stderr, "Normalization factor = %4.3f\n", normfactor);
		}
		
	}
	
	Histogram h;
	Histogram_create(&h, 10, 10);
	
	
	dist = (long*)calloc(100, sizeof(long));
	long numpos = 0;
	
	for (c=0; c<numchroms; c++) {
		
		if (chrdata == 0)
			chrname = 0; //(char*)(chroms[c]);
		else // means that chrdata was read from file
			chrname = chrnames[c];
		
		if (genome == 0) {
			chrlen      = chrlens[c];
		} else {
			die("not supportted yet\n");
		}
		
		// alloc counter vectors
		chip_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
		if (chip_counts == 0) {	    
			printf("Problem allocating chip_counts (chrlen=%ld).\n", chrlen);
			exit(1);
		}
		
		sprintf(readfile, "%s/reads.%s", chipdir, chrname); 
		
		if (!file_exists(readfile))
			continue;
		
		numreadsclonal_chip = 0;
		numreads_chip       = 0;
		getCountFromReadFile(readfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &chip_counts, &numreads_chip, &numreadsclonal_chip);
		
		if (mapdir != 0) {
			sprintf(mapfile, "%s/map.%s", mapdir, chrname);
			mapmap = create_binarized_array(chrlen);
			int mapb = readMapData(mapfile, &mapmap, mapt);
			printf("mapb = %d, out of %ld\n", mapb, chrlen);
		}
		
		int normrc = 0;
		int mapb2 = 0;
		
		for (i=0; i<chrlen; i++) {
			tmpval = (float)(chip_counts[i]);
			if (normalize == 1)
				tmpval *= normfactor;
			normrc = (int)(0.5 + tmpval);
			if (normrc < 100) {
				dist[ normrc ] ++;
			}       
			
			if ((mapdir == 0) || ((mapdir != 0) && (get_entry_in_binarized_array(mapmap, i) > 0))) {	
				mapb2 ++;
				Histogram_add(&h, tmpval); //(float)normrc);
				//Histogram_add(&h, (float)normrc);
			}
		}
		//printf("mapb2 = %d, out of %ld\n", mapb2, chrlen);
		numpos += chrlen;
		
		if (mapdir != 0)
			free(mapmap);
		free(chip_counts);
		
	}  // loop over chr
	
	//Histogram_print(&h);
	
	//printf("q=0.25 = %f\n", Histogram_quantile(&h, 0.25));
	//printf("q=0.50 = %f\n", Histogram_quantile(&h, 0.50));
	//printf("q=0.75 = %f\n", Histogram_quantile(&h, 0.75));
	
	printf("COVERAGE\tNUCLEOTIDES\t%%NUCLEOTIDES\tNUCLEOTIDES(>=X)\n");
	
	for (i=0; i<=10; i++) {
		
		long hcnm2 = 0;

		for (j=i; j<100; j++) {
			hcnm2 += dist[j];
		}

		printf("%ldX\t%ld\t%3.2f%%\t%ld\n", i, dist[i], 100*dist[i]/(double)numpos, hcnm2);
	}
	long hcnm = 0;
	for (i=11; i<100; i++) {
		hcnm += dist[i];
	}
	printf(">10X\t%ld\t%3.2f%%\n", hcnm, 100*hcnm/(double)numpos);

	return 0;	
}



