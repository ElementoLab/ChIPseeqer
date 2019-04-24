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

int main(int argc, char** argv) {
	
	unsigned short*	chip_counts;
	char			readfile[1000];
	char*			chipdir			= 0;
	long			chrlen			= 0;
	int				fraglen         = -1;
	int				readlen         = 0;
	int				numreads_chip   = 0;
	int				c				= 0;
	int				verbose			= 0;
	char*			chrname			= 0;
	char**          chrnames        = 0;
	int				format			= 0;
	char*			formattext		= 0;
	char*           chrdata         = 0;
	int*            chrlens         = 0;
	int				numchroms		= 24;
	int				totalnumreads	= 0;
	int				genome			= 0;
	float*			chrmap			= 0;
	int				uniquereads		= 1;
	
	int				totalnumreadsclonal	= 0;
	int				numreadsclonal_chip = 0;
	char*          singlereadfile = 0;

	if (argc < 2) {
	  die("Usage: ChIPseeqerGetNumClonalReads -chipdir DIR -chrdata FILE -format STR\n");
	}
	
	if (exist_parameter(argc, argv, "-chipdir"))
	  chipdir  = get_parameter(argc, argv, "-chipdir");

	if (exist_parameter(argc, argv, "-singlereadfile"))
	  singlereadfile  = get_parameter(argc, argv, "-singlereadfile");

	if (exist_parameter(argc, argv, "-verbose"))
	  verbose = atoi(get_parameter(argc, argv, "-verbose"));

	if (exist_parameter(argc, argv, "-format")) {
	  formattext = get_parameter(argc, argv, "-format");
	  format     = formatToId(formattext);				
	}
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata  = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);
		genome   = 0;
	} else {
		die("Please specify -chrdata\n");
	}
	
	//for each chromosome
	for (c=0; c<numchroms; c++) {
		
		if (chrdata == 0)
			chrname = 0;
		else // means that chrdata was read from file
			chrname = chrnames[c];
		
		if (genome == 0) {
			chrlen      = chrlens[c];
		} else {
			die("not supportted yet\n");
		}
		
		// allocate chip_count vector
		chip_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
		if (chip_counts == 0) {	    
			printf("Problem allocating chip_counts (chrlen=%ld).\n", chrlen);
			exit(1);
		}
	
		if (singlereadfile != 0) 
		  sprintf(readfile, "%s", singlereadfile);
		else 
		  sprintf(readfile, "%s/reads.%s", chipdir, chrname); 
		
		if (!file_exists(readfile)) {
			free(chip_counts);
			continue;
		}
		
		printf("# processing %s\n", readfile);

		numreadsclonal_chip = 0;
		numreads_chip       = 0;
		getCountFromReadFile(readfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &chip_counts, &numreads_chip, &numreadsclonal_chip);
		
		if (verbose == 1) {
			printf("chr %s, numclonal = %d, total num reads = %d\n", chrname, numreadsclonal_chip, numreadsclonal_chip + numreads_chip);
		}
		
		totalnumreads       += numreads_chip;
		totalnumreadsclonal += numreadsclonal_chip;
		
		free(chip_counts);
		
	}
	
	printf("clonal fraction that will be removed = %4.1f%% (%d/%d)\n", 100.0 * totalnumreadsclonal / ((double)totalnumreads + totalnumreadsclonal), totalnumreadsclonal, totalnumreads+totalnumreadsclonal);
	
	/* cleanup */
	for (c=0; c<numchroms; c++)
		free(chrnames[c]);
	free(chrnames);
	free(chrlens);
	free(chrmap);
	
	return 0;
}
