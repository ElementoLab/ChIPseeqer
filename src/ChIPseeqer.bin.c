#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dataio.h"
#include "statistics.h"
#include "sequences.h"
#include "protos.h"


#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

#define mynmax 100000


//struct defining a Genomic Interval
typedef struct _GenInt {
	int i;					//start position
	int j;					//end position
	int l;					//length
	char* c;				//chromosome name
	char* id;				//interval id
	float score;			//interval score
	double avg_chip;		//average height of ChIP peaks
	double avg_input;		//average height of input DNA peaks
	double fold;			//defines how much higher ChIP peaks should be, compared to input DNA peaks
	double avg_input_scaled;//average height of scaled input DNA peaks
	double avglpc;			//average p-value for each nucleotide in a ChIP peak region
	double avglpi;			//average p-value for each nucleotide in an input DNA peak region
	double avglpr;			//average p-value for each nucleotide in the normalized peak region
} GenInt;

int main(int argc, char** argv) {
	
	int					c				= 0;	//variable for the chromosomes loop
	long				i				= 0.0;	//variable for nucleotides loop
	int					k				= 0;	//variable for intervals loop
	unsigned short*		chip_counts		= NULL;	//vector that holds the number of ChIP peaks per position	
	unsigned short*		input_counts	= NULL; //vector that holds the number of input DNA peaks per position
	char*				chipfile		= 0;	//name of the file that contains the ChIP reads
	char*				inputfile		= 0;	//name of the file that contains the input DNA reads
	
	//2 new variables
	char*				chipdir         = 0;	//name of the directory that contains the ChIP files
	char*				inputdir        = 0;	//name of the directory that contains the input DNA files
	
	//char*				chr             = 0;	//chromosome name
	long				chrlen          = 0;	//chromosome length
	int					fraglen         = 170;	//fragment length
	int					readlen         = 36;	//read length
	int					numreads_chip   = 0;	//number of ChIP peaks per position
	int					numreads_input  = 0;	//number of input DNA peaks per position
	
	int					t               = 15;	//threshold
	int					t_prev          = 0;	//previous position in the peak
	int					t_start         = 0;	//start position of the peak
	int					t_end           = 0;	//end position of the peak
	
	int					t_sumnumreads_chip  = 0;	//sum of ChIP reads per position
	int					t_sumnumreads_input = 0;	//sum of input DNA reads per position
	
	double				t_sum_lpr       = 0.0;	//sum of normalized reads p-values per position
	double				t_sum_lpi       = 0.0;	//sum of input DNA reads p-values per position
	double				t_sum_lpc       = 0.0;	//sum of ChIP reads p-values per position
	
	double				fold_t          = 2.0;	//fold threshold
	
	int					verbose         = 0;		//verbose mode (value=1) displays program running details
	char*				chrname         = 0;		//chromosome name
	double				mappability     = 1.0;		//chromosome mappability
	int					format          = 0;		//format type of the data (default value 0: eland)
	char*				formattext      = "eland";	//format name of the data
	int					genome			= 1;		//version of genome (default value 1: hg18)
	char*				genome_text     = 0;		//version name of genome
	int					countreads		= 0;		//(value=1) creates a file that contains the number of reads per position
	
	GenInt*				a_int;						//variable of GenInt type
	int					maxnumint       = 1000000;	//maximum number of allowed intervals
	int					numint          = 0;		//number of intervals
	int					len				= 0;		//interval length
	double				avg_chip		= 0.0;		//average height of ChIP peaks
	double				avg_input		= 0.0;		//average height of input DNA peaks
	double				fold			= 0.0;		//defines how much higher ChIP peaks should be, compared to input DNA peaks
	int					mindist         = 100;		//minimum distance between peaks
	int					minlen          = 100;		//minimum peak width
	char*		        outfile			= 0;
	FILE*				countReadsFile	= 0;
	char*               chrdata         = 0;
	char**              chrnames        = 0;
	int*                chrlens         = 0;
	float*              chrmap          = 0;
	int                 ext				= 0;
	int					defmaxheight	= 0;
	//more new variables
	char*				chroms[] = {"chr1",		//array that holds the chromosomes names
		"chr10",
		"chr11",
		"chr12",
		"chr13",
		"chr14",
		"chr15",
		"chr16",
		"chr17",
		"chr18",
		"chr19",
		"chr2",
		"chr20",
		"chr21",
		"chr22",
		"chr3",
		"chr4",
		"chr5",
		"chr6",
		"chr7",
		"chr8",
		"chr9",
		"chrX",
		"chrY"};
	
	int	numchroms			= 24;//sizeof(chroms)/sizeof(*chroms);		//number of chromosomes
	const char*	nameformat	= "%s/reads.%s";
	int uniquereads			= 1;
	int j					= 0;
	int posmaxheight		= -1;
	int maxheight			= 0;
	int numpeaks			= 0;
	int totalnumreads_chip  = 0;
	int totalnumreads_input = 0;
	int totalpeakheight		= 0;
	int totalpeaksize		= 0;
	int numreadsclonal_chip = 0;
	
	//printf("Number: %lu", sizeof(chroms)/sizeof(*chroms));
	
	//handling lack of arguments
	if (argc < 5) {
		printf("Usage: ChIPseeqer.bin -chipdir DIR [-inputdir DIR] [options]\n\n");
		printf("where:\n");
		printf("-chipdir DIR \t contains the ChIP files (MANDATORY)\n"); 
		printf("-inputdir DIR \t contains the input files\n");
		printf("-outfile STR \t the name of the output file (MANDATORY)\n\n");
		printf("Options:\n");
		printf("-t FLOAT \t specifies the threshold (log p-value ratio). Default is 15 (stringent)\n");
		printf("-fold_t FLOAT \t specifies the fold threshold. Default is 2\n");
		printf("-fraglen INT \t specifies the average size of the DNA fragments whose extremities were sequenced\n");
		printf("-readlen INT \t specifies the average size of the DNA reads\n");
		printf("-format STR \t format of the read files: mit, sam, bowtiesam, bed, or eland. Default is eland\n");
		printf("-minlen STR \t minimum size of targets\n");
		printf("-mindist STR \t minimum distance between targets (merge otherwise)\n");
		printf("-genome STR \t genome version: hg18, hg19. Default is hg18\n");
		printf("-countreads INT \t creates a file that contains the numbers of reads per position\n");
		printf("-uniquereads INT \t 1 merges clonal reads into a single read\n");
		printf("-minpeakheight INT \t minimum peak height (reads count at the peak summit). Default is 0 (no mimimum - the pvalue calculation decides what a peak is).");
		printf("\n");
		exit(0);
	}
	
	
	//initialize variables with the arguments values
	
	if (exist_parameter(argc, argv, "-chipdir"))
		chipdir  = get_parameter(argc, argv, "-chipdir");
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata  = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);
		genome   = 0;
	}
	
	if (exist_parameter(argc, argv, "-inputdir"))
		inputdir  = get_parameter(argc, argv, "-inputdir");
	
	if (exist_parameter(argc, argv, "-t"))
		t  = atoi(get_parameter(argc, argv, "-t"));
	
	if (exist_parameter(argc, argv, "-uniquereads"))
		uniquereads  = atoi(get_parameter(argc, argv, "-uniquereads"));
	
	if (exist_parameter(argc, argv, "-fold_t"))
		fold_t  = atof(get_parameter(argc, argv, "-fold_t"));
	
	if (exist_parameter(argc, argv, "-fold"))
		fold_t  = atof(get_parameter(argc, argv, "-fold"));
	
	if (exist_parameter(argc, argv, "-f"))
		fold_t  = atof(get_parameter(argc, argv, "-f"));
	
	if (exist_parameter(argc, argv, "-ext"))
		ext  = atoi(get_parameter(argc, argv, "-ext"));
	
	if (exist_parameter(argc, argv, "-fraglen"))
		fraglen  = atoi(get_parameter(argc, argv, "-fraglen"));
	
	if (exist_parameter(argc, argv, "-readlen"))
		readlen  = atoi(get_parameter(argc, argv, "-readlen"));
	
	if (exist_parameter(argc, argv, "-format")) {
		formattext = get_parameter(argc, argv, "-format");
		format     = formatToId(formattext);			
	}
	
	if (exist_parameter(argc, argv, "-minlen"))
		minlen  = atoi(get_parameter(argc, argv, "-minlen"));
	
	if (exist_parameter(argc, argv, "-mindist"))
		mindist  = atoi(get_parameter(argc, argv, "-mindist"));
	
	if (exist_parameter(argc, argv, "-genome")) {
		genome_text = get_parameter(argc, argv, "-genome");
		if (strcmp(genome_text, "hg18") == 0)
			genome = 1;
		else if (strcmp(genome_text, "hg19") == 0)
			genome = 2;
	}
	
	if (exist_parameter(argc, argv, "-verbose"))
		verbose = atoi(get_parameter(argc, argv, "-verbose"));
	
	if (exist_parameter(argc, argv, "-countreads"))
		countreads = atoi(get_parameter(argc, argv, "-countreads"));
	
	if (exist_parameter(argc, argv, "-minpeakheight"))
		defmaxheight  = atoi(get_parameter(argc, argv, "-minpeakheight"));
	
	
	//OPEN COUNTREADSFILE
	if (countreads == 1) {
		countReadsFile = fopen ("COUNT_READS.txt","w");
	}
	
	if (exist_parameter(argc, argv, "-outfile")) {
		outfile = get_parameter(argc, argv, "-outfile");
	}
	
	if (outfile == 0) {
		fprintf(stdout, "The -outfile option was not defined. Please try again.\n");
		die("");
	}
	
	if ((verbose == 1) && (chipdir != 0)) {
	    fprintf(stdout, "reading reads from %s\n", chipdir);
	}	
	if ((verbose == 1) && (inputdir != 0)) {
		fprintf(stdout, "reading reads from %s\n", inputdir);
	}        
	
	fprintf(stdout, "Starting peak detection process, with parameters:\n");
	fprintf(stdout, "Reads format: %s\n", formattext);
	fprintf(stdout, "Threshold: %d\n", t);
	fprintf(stdout, "Clonal read collapse: %d\n", uniquereads);
	fprintf(stdout, "Fold threshold: %f\n", fold_t);
	fprintf(stdout, "Fragment length: %d\n", fraglen);
	fprintf(stdout, "Reads length: %d\n", readlen);
	fprintf(stdout, "Minimum size of peak (bp): %d\n", minlen);
	fprintf(stdout, "Minimum distance between peaks (bp): %d\n", mindist);
	if (chrdata == 0) {
	  fprintf(stdout, "Using hg18 chromosome data\n");
	} else {
	  fprintf(stdout, "Using user-provided chromosome data from %s\n", chrdata);
	}	
	
	FILE* outputFile = fopen(outfile, "w");
	
	if (outputFile==NULL)
	{
		fprintf(stdout, "Could not open outfile: %s. Please try again.\n", outfile);
		exit(1);	
	}
	
	 /*FILE* myFile = 0;
	 myFile = fopen(blah);
	 if (myFile == 0) {
	 //something went wrong while opening
	 }
	 
	 FILE * pFile;
	 pFile = fopen ("myfile.txt","w");
	 if (pFile!=NULL)
	 {
	 fputs ("fopen example",pFile);
	 fclose (pFile);
	 }
	 
	*/
	
	//for each chromosome
	for (c=0; c<numchroms; c++) {	
		if (chrdata == 0)
			chrname = (char*)(chroms[c]);
		else // means that chrdata was read from file
			chrname = chrnames[c];
		
		//if(strcmp(chrname, "chr10") != 0) continue;
		
		if (countreads == 0) {
			fprintf(stdout, "Processing chromosome: %s\n", chrname);
		}
		
		if (countReadsFile) {
			fprintf(countReadsFile, "Chromosome: %s\n", chrname);
		}
		
		// chr data read from file
		if (genome == 0) {
			chrlen      = chrlens[c];
			if (chrmap[c] >= 0) {
				mappability = chrmap[c];
			} else {
				mappability = 1.0;
			}
		} else
			if (genome == 1) {
				chrlen = hg18_chrlen(chrname);
				if (chrlen < 0) {
					printf("Problem ... chr %s unrecognized.\n", chrname);
				}
				
				if (exist_parameter(argc, argv, "-mappability"))
					mappability = atof(get_parameter(argc, argv, "-mappability"));
				else 
					mappability = hg18_30mer_mappability(chrname);	  
			}
			else if (genome == 2) {
				chrlen = hg19_chrlen(chrname);
				if (chrlen < 0) {
					printf("Problem ... chr %s unrecognized.\n", chrname);
				}
				
				if (exist_parameter(argc, argv, "-mappability"))
					mappability = atof(get_parameter(argc, argv, "-mappability"));
				else 
					mappability = hg19_30mer_mappability(chrname);
			}
		
		if (chrlen < 0) {
			printf("Problem...chromosome %s unrecognized.\n", chrname);
		}
		
		//alloc the file name pointers
		chipfile	= (char*)calloc(strlen(chipdir) + strlen(chrname) + 9, sizeof(char));
		if (inputdir != 0) {
			inputfile	= (char*)calloc(strlen(inputdir) + strlen(chrname) + 9, sizeof(char));
		}
		
		//get the file names
		sprintf(chipfile, nameformat, chipdir, chrname);
		
		if (inputdir !=  0) {
			sprintf(inputfile, nameformat, inputdir, chrname);
		}
		
		
		if (!file_exists(chipfile)) {
			printf("Error: ChIPfile (%s) not found. Try again.\n", chipfile);
			exit(0);
		}
		
		if (!file_exists(inputfile)) {
			if(inputdir != NULL) {
				printf("Warning: input read file not found (ignored)\n");
			}
			inputfile = 0;
		}
		
		// alloc counter vectors
		chip_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
		if (chip_counts == 0) {
			die("Problem allocating chip_counts.\n");
		}
		input_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
		if (input_counts == 0) {
			die("Problem allocating input_counts.\n");
		}
		
		if (verbose == 1)
			fprintf(stdout, "%s\t%ld\t%f\n", chrname, chrlen, mappability); 
		
		//read the files
		getCountFromReadFile(chipfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &chip_counts, &numreads_chip, &numreadsclonal_chip);
		totalnumreads_chip += numreads_chip;
		
		if (inputfile != 0) {
			getCountFromReadFile(inputfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &input_counts, &numreads_input, &numreadsclonal_chip);
			totalnumreads_input += numreads_input;
		}
		
		
		//computes the number of read counts for each nucleotide
		/*if (countfile != 0) {
		 
		 int* readcounts;
		 readcounts = (int*)calloc(100, sizeof(int));
		 
		 for (i=0; i<chrlen; i++) {
		 if (chip_counts[i] < 100)
		 readcounts[ chip_counts[i] ] ++;
		 }
		 
		 cfp = fopen(countfile, "w");
		 if (!cfp) {
		 die("Could not open countfile.\n");
		 }
		 
		 for (i=0; i<50; i++) {
		 fprintf(cfp, "%ld\t%d\n", i, readcounts[i]);
		 }
		 
		 fclose(cfp);
		 }*/
		
		double lambda_chip;
		double lambda_input;
		
		if (fraglen != 0) {
			lambda_chip  = numreads_chip  * (double)fraglen / ( (double)chrlen * mappability );
			lambda_input = numreads_input * (double)fraglen / ( (double)chrlen * mappability );
		}
		else {
			lambda_chip  = numreads_chip / ( (double)chrlen * mappability );
			lambda_input = numreads_input / ( (double)chrlen * mappability );
		}
		
		double scale_numreads = numreads_chip / (double)numreads_input;
		
		int    cntz  = 0;
		double lpr_t = log10(pow(10,-t));
		
		//"hash" table in order to not have to recalculate too many log10(igam(etc))
		char*   a_cc_i = (char*)calloc(mynmax, sizeof(char));
		double* a_cc_p = (double*)calloc(mynmax, sizeof(double));
		
		char*   a_ci_i = (char*)calloc(mynmax, sizeof(char));
		double* a_ci_p = (double*)calloc(mynmax, sizeof(double));
		
		
		a_int	= (GenInt*)malloc( maxnumint * sizeof(GenInt));
		numint	= 0;
		
		//for each position (nucleotide)
		for (i=0; i<chrlen; i++) {
			
			int cc = chip_counts[i];
			int ci = input_counts[i];
			
			if (cc == 0) {
				cntz++;
			} else {
				
				if (countReadsFile) {
					fprintf(countReadsFile, "%ld\t%d\t%d\n", i, cc, ci);
				}
				double lpc = 0;
				double lpi = 0;
				
				if (cc > 0) {
					if ((cc < mynmax) && (a_cc_i[cc] == 1)) 
						lpc = a_cc_p[cc];
					else {
						lpc = log10(igam((double)cc,lambda_chip));
						if (cc < mynmax) {
							a_cc_p[cc] = lpc;
							a_cc_i[cc] = 1;
						}
					}
				}
				
				if (ci > 0) {
					if ((ci < mynmax) && (a_ci_i[ci] == 1)) 
						lpi = a_ci_p[ci];
					else {
						lpi = log10(igam((double)ci,lambda_input));
						if (ci < mynmax) {
							a_ci_p[ci] = lpi;
							a_ci_i[ci] = 1;
						}
					}
				}
				
				//if (ci > 0)
				//lpi = log10(igam((double)ci,lambda_input));
				
				
				//double lpc = log10(cumpois(cc, lambda_chip));    
				//double lpi = log10(cumpois(ci, lambda_input));    
				
				/*
				 if (lpc < -16.0) {
				 lpc = -16.0;
				 }
				 
				 if (lpi < -16.0) {
				 lpi = -16.0;
				 }
				 */
				
				double lpr = lpc - lpi;
				
				//if the peak does not pass the threshold
				if (lpr < lpr_t) {
					
					//if not already in a peak
					if (t_prev == 0) {
						t_start = i;
						t_end   = i;
						
						t_sumnumreads_input  = ci;
						t_sumnumreads_chip   = cc;
						t_sum_lpr            = lpr;
						t_sum_lpi            = lpi;
						t_sum_lpc            = lpc;
						
					} else {
						t_end ++;    // extend
						t_sumnumreads_input += ci;
						t_sumnumreads_chip  += cc;
						t_sum_lpr           += lpr;
						t_sum_lpi           += lpi;	  
						t_sum_lpc           += lpc;
					}
					
					t_prev = 1;  // set for next pos
					
				} else { //if the peak passes the threshold
					
					if (t_prev == 1) {
						
						len       = t_end - t_start + 1;
						avg_chip  = t_sumnumreads_chip /(double)len;
						avg_input = t_sumnumreads_input/(double)len;
						fold      = -1;
						
						if (avg_input > 0) {
							fold = (avg_chip+lambda_input) / (avg_input*scale_numreads+lambda_input);
						} 
						
						//if the peak passes the fold threshold
						if ((fold < 0) || (fold > fold_t)) {
							
							a_int[numint].c                = strdup(chrname);
							a_int[numint].i                = t_start;
							a_int[numint].j                = t_end;
							a_int[numint].l                = len;
							a_int[numint].avg_chip         = avg_chip;
							a_int[numint].avg_input        = avg_input;
							a_int[numint].avg_input_scaled = avg_input * scale_numreads;
							a_int[numint].fold             = (avg_input>0?fold:0);
							a_int[numint].avglpc           = t_sum_lpc/(double)len;
							a_int[numint].avglpi           = t_sum_lpi/(double)len;
							a_int[numint].avglpr           = t_sum_lpr/(double)len;
							numint++;
						
							if (numint == maxnumint) {
								die("Attention: maxnumint reached ... dying (please modify and recompile)\n");
							}
							
						}
					}
					t_prev = 0;
				}
				
				
				/*if (verbose == 1) {
				 if (cntz > 0) {
				 printf(" -- %d positions with 0 count.\n", cntz);
				 }
				 
				 printf("%d\t%f\t%d\t%f", cc, lpc, ci, lpi);
				 printf("\t%f", lpr);
				 
				 //double newlpi = log10(cumpois(ci, lambda_input));
				 
				 //printf(" lpi=%f = log10(cumpois(ci=%d, lambda_input=%f)) .. %f ", lpi, ci, lambda_input, newlpi);    
				 
				 
				 //if (!isfinite(lpi)) {
				 //  printf(" (inf)");
				 //}
				 //if (lpc < log10(0.0000000001)) {
				 //  printf("\t******");
				 //}
				 printf("\n");
				 }*/
				
				cntz = 0;
			}
			
		}//loop over nucleotides
		
		int    start      = a_int[0].i;
		int    last_end   = a_int[0].j;
		char*  last_chr   = a_int[0].c;
		double bestscore;
		if (inputdir != 0) {
			bestscore  = a_int[0].avg_chip - a_int[0].avg_input_scaled;
		} else {
			bestscore  = a_int[0].avg_chip;
		}
		double bestlpr    = a_int[0].avglpr; 
		double score;
		
		for (k=1; k<numint; k++) {
			
			if (inputdir != 0) {
				score = a_int[k].avg_chip - a_int[k].avg_input_scaled;
			} else {
				score = a_int[k].avg_chip;		    
			}
			
			// if gap between peaks is greater than mindist, output prev interval, create new
			if ((a_int[k].i - last_end) > mindist) {
				
				len = last_end - start + 1;
				
				// choose or not to show it
				if ((minlen < 0) || ((minlen >= 0) && (len >= minlen))) {
					
					// determine max peak height
					posmaxheight = -1;
					maxheight    = 0;
					for (j=start; j<=last_end; j++) {
						if (chip_counts[j] > maxheight) {
							maxheight    = chip_counts[j];
							posmaxheight = j;
						}
					}
					
					int peak_len				= last_end - start;
					int mid						= start + peak_len/2;
					int summit_dist_from_mid	= abs(posmaxheight-mid); 
					
					if (maxheight > defmaxheight) {
						
						fprintf(outputFile, "%s\t%d\t%d\t%5.4f", last_chr, start-ext, last_end+ext, bestlpr);
						
						fprintf(outputFile, "\t%5.4f", bestscore);
						
						fprintf(outputFile, "\t%d\t%d\t%3.1f", posmaxheight, maxheight, 100.0*((double)posmaxheight-start)/(double)len);
						
						fprintf(outputFile, "\t%d\t%d\t%d", peak_len, mid, summit_dist_from_mid);
						
						fprintf(outputFile, "\n");
						
						totalpeakheight += maxheight;
						totalpeaksize	+= peak_len;
						
						numpeaks ++;
					}
					
				}
				
				// in any cases,  new interval
				start      = a_int[k].i;
				bestscore  = score;
				bestlpr    = a_int[k].avglpr;        
				
			} else {
				
				// incorporate into prev interval      
				if (score > bestscore) {
					bestscore = score;
					bestlpr   = a_int[k].avglpr;
				}
			}
			
			// stuff to do in all cases
			last_chr   = a_int[k].c; 
			last_end   = a_int[k].j;
		} // for k=1; k<numint
	
		/* cleanup */
		for (k=0; k<numint; k++)
			free(a_int[k].c);
		
		free(a_int);
		free(a_cc_i);
		free(a_cc_p);
		free(a_ci_i);
		free(a_ci_p);
		free(chip_counts);
		free(input_counts);
		free(chipfile);
		if (inputfile != 0) {
			free(inputfile);
		}
		
		a_int	= 0;
		a_cc_i	= 0;
		a_cc_p	= 0;
		a_ci_i	= 0;
		a_ci_p	= 0;
		chip_counts		= 0;
		input_counts	= 0;
		chipfile		= 0;
		inputfile		= 0;
		fprintf(stdout, "%d peaks detected so far\n", numpeaks);
		
		
	}//loop over chromosomes
	
	//CLOSE COUNTREADSFILE
	if (countReadsFile) {
		fclose(countReadsFile);
	}
	fclose(outputFile);
	
	fprintf(stdout, "%d peaks detected in total\n", numpeaks);	
	
	if(numpeaks != 0) {
		
		float avgpeakheight		= totalpeakheight/numpeaks;
		float avgpeaksize		= totalpeaksize/numpeaks;
		
		fprintf(stdout, "%3.1f Average peak height (reads)\n", avgpeakheight);	
		fprintf(stdout, "%3.1f Average peak size (bp)\n", avgpeaksize);	
		
	}
	fprintf(stdout, "Used %d ChIP reads", totalnumreads_chip);
	if (inputdir != 0) {
		fprintf(stdout, " and %d input reads", totalnumreads_input);
	}
	fprintf(stdout, "\n");

	/* cleanup */
	for (k=0; k<numchroms; k++)
		free(chrnames[k]);
	free(chrnames);
	free(chrlens);
	free(chrmap);
	
	return 0;
}

