#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <search.h>

#include "dataio.h"
#include "statistics.h"
#include "sequences.h"
#include "postscript.h"
#include "hashtable.h"

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

#define  NBPEAKS 40000
#define  NBWINDOWS 1000

//struct defining a Genomic Interval
typedef struct _GenInt {
	int i;
	int j;
	//char* c;
	//char* id;
	float score;
	int del;
	//int strand;
	char* str;
	char* tsname;
	char* genename;
	char* chr;
	int   strand;
} GenInt;

int CmpGenIntStart(const void* _a, const void* _b);

void readFileSumUpColumnPerChrom(char* file, int col, int colchr, int** counts, int** countgenes);
void readRefGenePromotersOrDownstream(char* file, GenInt** a_proms, int* numproms, int up, int dn, int region);

int main(int argc, char** argv) {
	
	// general
	long i;
	long k;
	int c;
	
	unsigned short* chip_counts;	//vector that holds the number of ChIP peaks per position
	
	char*	chipfile		= 0;	//name of the file that contains the ChIP reads
	char*	chipdir         = 0;	//name of the directory that contains the ChIP files
	
	long  chrlen			= 0; 
	int   fraglen			= 0;
	int   readlen			= -1;  //OE 
	int   numreads_chip		= 0;
	int   numreads_input	= 0;
	int   verbose			= 0;
	char* chrname			= NULL;
	int    format			= 0;
	char*  formattext		= NULL;
	int     ws				= 10;
	long    j				= 0.0;
	int     from			= 0;
	int		to				= 0;
	char*   desc			= 0;
	char*   intervals		= 0;
	FILE*   f1				= NULL;
	GenInt* a_int			= NULL;
	int     maxnumint		= 10000000;
	int     numint			= 0;
	int     hasid			= 0;
	int		maxnorm			= 0;
	int     rpkmnorm        = 0;   // RPKM normalzs
	int		numntnorm	    = 0;
	int     totalnumreads   = 0;
	int     numreads		= 0;
	double  normfactor	    = 1.0;
	int     normto			= 10000000;   // 10 M reads
	
	// reading file
	char*	buff			= NULL;
	int		mynmax			= 100000;
	char**	p				= NULL;
	int		m				= 0;
	char*	tmpid			= NULL;
	int		ext				= 0;
	int		idx				= 0;
	
	
	const char*		nameformat	= "%s/reads.%s";
	int				uniquereads	= 1;
	long int		numnt		= 0;
	long int		totalnumnt	= 0;
	
	char*			outfile		= NULL;
	FILE*			fpout		= NULL;
	
	int				numchroms	= 0;		
	char*           chrdata     = NULL;
	char**          chrnames    = NULL;
	int*            chrlens     = NULL;
	float*          chrmap      = NULL;
	
	int             numreadsclonal_chip = 0;	
	
	float**	data		= NULL;
	char*	outepsmap	= NULL;	
	int		nbcond		= 0;
	int		nborfs		= 0;
	int		peak		= 0;
	int		bin			= 0;
	char*	xlabel		= NULL;
	char*	ylabel		= NULL;	
	
	
	int   refgene		= 0;
	int   dn			= 2000;
	int   up			= 2000;
	int   tsregion		= 0;
	
	//handling lack of arguments
	if (argc < 2) {
		die("Usage: ChIPseeqerGetReadDensityProfiles.bin -chipdir DIR -ws INT -intervals FILE\n");
	}
	
	//initialize variables with the arguments values
	if (exist_parameter(argc, argv, "-ws"))
		ws = atoi(get_parameter(argc, argv, "-ws"));
	
	if (exist_parameter(argc, argv, "-w"))
		ws = atoi(get_parameter(argc, argv, "-w"));
	
	if (exist_parameter(argc, argv, "-ext"))
		ext = atoi(get_parameter(argc, argv, "-ext"));
	
	if (exist_parameter(argc, argv, "-uniquereads"))
		uniquereads  = atoi(get_parameter(argc, argv, "-uniquereads"));
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals = get_parameter(argc, argv, "-intervals");
	
	
	// refgene
	if (exist_parameter(argc, argv, "-refgene")) {
		intervals = get_parameter(argc, argv, "-refgene");
		refgene = 1;
	}
	if (exist_parameter(argc, argv, "-dn"))
		dn = atoi(get_parameter(argc, argv, "-dn"));  
	if (exist_parameter(argc, argv, "-up"))
		up = atoi(get_parameter(argc, argv, "-up"));
	
	// TES or TSS
	if (exist_parameter(argc, argv, "-tsregion")) {
		char* tsregiontxt = get_parameter(argc, argv, "-tsregion");
		if (strcmp(tsregiontxt, "TSS") == 0)
			tsregion = 0;
		else if  (strcmp(tsregiontxt, "TES") == 0)
			tsregion = 1;
		else 
			die("-tsregion is either TES or TSS\n");
	}
	// same just in case (region instead of tsregion)
	if (exist_parameter(argc, argv, "-region")) {
		char* tsregiontxt = get_parameter(argc, argv, "-region");
		if (strcmp(tsregiontxt, "TSS") == 0)
			tsregion = 0;
		else if  (strcmp(tsregiontxt, "TES") == 0)
			tsregion = 1;
		else 
			die("-region is either TES or TSS\n");
	}
	
	
	if (exist_parameter(argc, argv, "-chipdir"))
		chipdir  = get_parameter(argc, argv, "-chipdir");
	
	if (exist_parameter(argc, argv, "-chrname"))
		chrname  = get_parameter(argc, argv, "-chrname");
	
	if (exist_parameter(argc, argv, "-fraglen"))
		fraglen  = atoi(get_parameter(argc, argv, "-fraglen"));
	
	if (exist_parameter(argc, argv, "-verbose"))
		verbose = atoi(get_parameter(argc, argv, "-verbose"));
	
	if (exist_parameter(argc, argv, "-desc"))
		desc = get_parameter(argc, argv, "-desc");
	
	if (exist_parameter(argc, argv, "-outfile"))
		outfile = get_parameter(argc, argv, "-outfile");
	
	if (exist_parameter(argc, argv, "-format")) {
		formattext = get_parameter(argc, argv, "-format");
		format     = formatToId(formattext);			
	}
	
	fprintf(stdout, "# format %d\n", format);
	
	if (exist_parameter(argc, argv, "-maxnorm"))
		maxnorm = atoi(get_parameter(argc, argv, "-maxnorm"));
	
	if (exist_parameter(argc, argv, "-rpkmnorm"))
		rpkmnorm = atoi(get_parameter(argc, argv, "-rpkmnorm"));
	
	if (exist_parameter(argc, argv, "-normalize"))
		rpkmnorm = atoi(get_parameter(argc, argv, "-normalize"));
	
	
	if (exist_parameter(argc, argv, "-numntnorm"))
		numntnorm = atoi(get_parameter(argc, argv, "-numntnorm"));
	
	if (exist_parameter(argc, argv, "-from"))
		from = atoi(get_parameter(argc, argv, "-from"));
	else 
		from = 0;
	
	if (exist_parameter(argc, argv, "-to"))
		to = atoi(get_parameter(argc, argv, "-to"));
	else 
		to = chrlen;
	
	if (exist_parameter(argc, argv, "-hasid"))
		hasid = atoi(get_parameter(argc, argv, "-hasid"));
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata  = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);	  
	} else {
		
		get_hg18_chroms(&chrnames, &numchroms, &chrlens);		  
	}
	
	if (exist_parameter(argc, argv, "-outepsmap")) {
		outepsmap =      get_parameter(argc, argv, "-outepsmap");
		fprintf(stdout, "# output to %s\n", outepsmap);
	} else {
		outepsmap = 0;
	}
	
	if (exist_parameter(argc, argv, "-xlabel"))
		xlabel  = get_parameter(argc, argv, "-xlabel");
	
	if (exist_parameter(argc, argv, "-ylabel"))
		ylabel  = get_parameter(argc, argv, "-ylabel");
	
	if(xlabel == NULL) {
		xlabel = "Distance to peak summit";
	}
	
	if(ylabel == NULL) {
		ylabel = "Average read density";
	}
	
	
	if (rpkmnorm == 1) {
		fprintf(stdout, "# RPKM-normalize to %d reads\n", normto);
	} else {
		fprintf(stdout, "# no RPKM-normalization\n");		  
	}
	
	fprintf(stdout, "# ws = %d\n", ws);
	
	// allocate a huge chunk of memory
	data = (float**)malloc(NBPEAKS * sizeof(float*));
	for (i=0; i<NBPEAKS; i++) {
		data[i] = (float*)malloc(NBWINDOWS * sizeof(float));
	}
	
	
	if (refgene == 1) {
		
		readRefGenePromotersOrDownstream(intervals, &a_int,  &numint, up, dn, tsregion);
		fprintf(stdout, "# Loaded refGene data ... found %d regions\n", numint);
		hasid = 1;  // set flag
	} else {
		
		//
		// start reading intervals
		//
		a_int = (GenInt*)malloc(maxnumint * sizeof(GenInt));
		if (a_int == 0) {
			die("Problem allocating a_int.\n");
		}
		
		buff  = (char*)malloc(mynmax * sizeof(char));
		
		// read first set of intervals
		numint	= 0;
		f1		= fopen(intervals, "r");
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
			
			a_int[numint].tsname	= tmpid; 			
			a_int[numint].chr		= strdup(p[idx]);
			a_int[numint].i			= atoi(p[idx+1]);
			a_int[numint].j			= atoi(p[idx+2]);
			if (hasid == 1) {
				a_int[numint].strand  = atoi(p[idx+3]);
			}
			a_int[numint].score = 0; 
			a_int[numint].del	= 0;
			numint++;
			
			
			if (numint == maxnumint) {
				die("Max number of intervals reached.\n");
			}
			free(p);
			//free(tmpid);
		}
		
		//
		// end read intervals
		//
	}
	
	//alloc the file name pointers
	chipfile	= (char*)calloc(1000, sizeof(char));	    
	
	// calc total number of reads for RPKM-style 
	if (rpkmnorm == 1) {
		
		totalnumreads = 0;
		for (c=0; c<numchroms; c++) {
			
			chrname = (char*)(chrnames[c]);
			chrlen  = chrlens[c];
			sprintf(chipfile, "%s/reads.%s", chipdir, chrname); 
			if (!file_exists(chipfile)) {
				fprintf(stdout, "File %s does not exist, skip.\n", chipfile); 
				continue;
			}
			if (uniquereads == 0) 
				numreads = CountReads(chipfile, format);
			else 
				numreads = CountNonClonalReads(chipfile, format, chrlen);
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
		
	} // if normalize by total number of aligned nucleotides
	else if (numntnorm == 1) {
		
		totalnumnt = 0;
		for (c=0; c<numchroms; c++) {
			
			chrname = (char*)(chrnames[c]);
			sprintf(chipfile, "%s/reads.%s", chipdir, chrname); 
			if (!file_exists(chipfile)) {
				fprintf(stdout, "File %s does not exist, skip.\n", chipfile); 
				continue;
			}
			
			numnt = CountAlignedNucleotides(chipfile, format);
			if (verbose == 1) {
				fprintf(stdout, "Found %ld aligned nt in chr %s\n", numnt, chrname);
			}
			totalnumnt += numnt;
			
		}
		if (verbose == 1) {
			fprintf(stdout, "# Found a total of %ld aligned nucleotides\n", totalnumnt);
		}
		
		normfactor =  normto / (double)totalnumnt;
		if (verbose == 1) {
			fprintf(stdout, "# Normalization factor = %4.3f\n", normfactor);
		}
		
	}
	
	
	if (outfile != 0) {
		fpout = fopen(outfile, "w");
		if (!fpout) {
			die("Cannot open outfile\n");
		}
	} else {
		fpout = stdout;
	}
	
	//for each chromosome
	for (c=0; c<numchroms; c++) {	
		
		chrname = (char*)(chrnames[c]);
		chrlen  = chrlens[c];
		
		if (chrlen < 0) {
			printf("Problem ... chromosome %s unrecognized.\n", chrname);
		}
		
		
		//get the file names
		sprintf(chipfile, nameformat, chipdir, chrname);
		if (!file_exists(chipfile)) {
			fprintf(stdout, "File %s does not exist, skip.\n", chipfile); 
			continue;
		}		
		if (verbose == 1) {
			fprintf(stdout, "# Reading reads from %s\n", chipfile);
		}
		
		// alloc counter vectors
		chip_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
		if (chip_counts == 0) {
			die("Problem allocating chip_counts.\n");
		}
		
		//printf("%s\n", chipfile);
		
		getCountFromReadFile(chipfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &chip_counts, &numreads_chip, &numreadsclonal_chip);
		
		
		if (verbose == 1)
			fprintf(stdout, "Number of reads = %d (chip) and %d (input)\n", numreads_chip, numreads_input);
		
		
		float	avg = 0;
		
		peak = 0;
		
		for (i=0; i<numint; i++) {
			
			float max = 0.0000000000001;
			
			if (strcmp(chrname, a_int[i].chr) == 0) {	      
				
				if (verbose == 1) {
					if (refgene == 1) {
						fprintf(stdout, "# ts = %s, strand = %d\n", a_int[i].tsname, a_int[i].strand);
					}
				}
				
				if ((hasid == 1) || (refgene == 1)) {
					fprintf(fpout, "%s", a_int[i].tsname);
				} else {
					fprintf(fpout, "%s-%d-%d", a_int[i].chr, a_int[i].i, a_int[i].j);
				}
				
				
				if(maxnorm == 1) {
					fprintf(stdout, "%s", a_int[i].tsname);
				}
				
				bin = 0;
				
				// go thru interval, ws nt at a time
				if ((hasid == 0) || (a_int[i].strand == 1)) {  // two cases, no gene id so strand does not matter, or 
					
					if (verbose == 1) {
						fprintf(stdout, "# Dealing with this in strand 1\n");
					}
					
					for (j=a_int[i].i-ext; j<a_int[i].j+ext; j+=ws) {
						avg = 0;
						for (k=j; k<min(chrlen,j+ws); k++) {
							//printf("\t + %3.1d", chip_counts[k]);
							avg += chip_counts[k];
						}
						avg /= ws;
						
						// RPKM-style norm
						avg *= normfactor;
						
						fprintf(fpout, "\t%3.1f", avg);
						
						data[peak][bin] = avg;
						
						//printf("1:\t(%d,%d):%3.4f", peak, bin, data[peak][bin]);
						
						if(maxnorm == 1) {
							if(avg > max && max != 0) {
								max = avg;
							}
						}
						bin++;
					}	  
					
					
					
				} else {	// other strand
					
					if (verbose == 1) {
						fprintf(stderr, "# Dealing with this in strand -1\n");
					}
					for (j=a_int[i].j+ext; j>a_int[i].i-ext; j-=ws) {
						avg = 0;
						
						for (k=j-ws; k<min(chrlen,j); k++) {
							//printf("\t - %3.1d", chip_counts[k]);
							avg += chip_counts[k];
						}
						avg /= ws;
						
						// RPKM-style norm
						//if ((normalize == 1) || (numntnorm == 1))
						avg *= normfactor;
						
						fprintf(fpout, "\t%3.1f", avg);
						
						data[peak][bin] = avg;
						
						//printf("1:\t(%d,%d):%3.4f", peak, bin, data[peak][bin]);
						
						if(maxnorm == 1) {
							if(avg > max && max != 0) {
								max = avg;
							}
						}
						bin++;
					}
				}
				
				// 
				if(maxnorm == 1) {
					
					if (a_int[i].strand == 1) {		
						
						for (j=a_int[i].i-ext; j<=a_int[i].j+ext; j+=ws) {
							avg = 0;
							for (k=j; k<min(chrlen,j+ws); k++) 
								avg += chip_counts[k];
							avg /= ws;
							avg /= max;
							fprintf(stdout, "\t%3.4f", avg);
						}	      
						
					} else {
						
						for (j=a_int[i].j+ext; j>=a_int[i].i-ext; j-=ws) {
							avg = 0;
							
							for (k=j-ws; k<min(chrlen,j); k++) 
								avg += chip_counts[k];
							avg /= ws;
							avg /= max;
							fprintf(stdout, "\t%3.4f", avg);
							
						}
					}
				}
				
				fprintf(fpout, "\n");
				if (maxnorm == 1)
					fprintf(stdout, "\n");
				
				peak++;
			} 
			
		} // end for loop over intervals    
		
		free(chip_counts);
		
	} //loop over chromosomes
	
	nborfs = peak;
	nbcond = bin;
	
	//printf("\n%d\t%d\n", peak, bin);
	
	if (outepsmap != 0) {
		printf("Writing EPS map ...");
		colMeans2Dplot(outepsmap, data, nborfs, nbcond, xlabel, ylabel, ws); 
		printf("Done.\n");
	}
	
	/* cleanup */
	if (outfile != 0) {
		fclose(fpout);
		fprintf(stdout, "# created %s\n", outfile);
	}
	fclose(f1);
	for (i=0; i<NBPEAKS; i++)
		free(data[i]);
	free(data);
	for (i=0; i<numint; i++) { 
		if(a_int[i].tsname != NULL)
			free(a_int[i].tsname);
		free(a_int[i].chr);
	}
	free(a_int);	
	for (c=0; c<numchroms; c++)
		free(chrnames[c]);
	free(chrnames);
	free(chrlens);
	free(chipfile);
	free(buff);
	
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



// region = 0 -> promoters, 1 -> downstream
void readRefGenePromotersOrDownstream(char* file, GenInt** a_proms, int* numproms, int up, int dn, int region)
{
	
	char*  buff;
	int    mynmax = 100000;
	char** a;
	int    m;
	FILE*  f;
	//int    cidx = -1;
	//int    i;
	char*  line = 0;
	int    hashret;
	ENTRY  e;
	ENTRY* ep;
	struct my_hsearch_data* hash_genes;
	
	// initialize num rnas per chr
	*numproms = nbLinesInFile(file);
	
	/* build a hash table */
	hash_genes = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
	hashret = my_hcreate_r((*numproms)*2, hash_genes);
	if (hashret == 0) {
		printf("Could not create hash table ...\n");
		exit(0);
	}
	
	// alloc intervals
	*a_proms = (GenInt*)malloc( (*numproms) * sizeof(GenInt));
	if (*a_proms == 0) {
		die("Problem allocating a_rnas.\n");
	}
	
	//
	// read set of intervals
	//
	buff  = (char*)malloc(mynmax * sizeof(char));  
	f = fopen(file, "r");
	if (f == 0) 
		die("Cannot open refGene file.\n");
	*numproms = 0;
	while (fgets(buff, mynmax, f) != 0) {
		chomp(buff);
		
		line = strdup(buff);
		
		split_line_delim(buff, "\t", &a, &m);
		
		char* n      = a[1];
		
		/* query */
		e.key = n;        
		my_hsearch_r(e, FIND, &ep, hash_genes) ;  
		if (ep) {
			/* success, already rgere */
			free(a);
			continue;
		} else {
			
			// add
			/* enter key/value pair into hash */
			e.key   = strdup( n );
			e.data  = (char*)(*numproms);
			hashret = my_hsearch_r(e, ENTER, &ep, hash_genes);
			if (hashret == 0) {
				printf("Could not enter entry into hash table ...\n");
				exit(0);
			}
		}
		
		char* c      = a[2];     //chr
		int strand = 0;
		char fr = a[3][0];
		if (fr == '+') {
			strand    = 1;
		} else if (fr == '-') {
			strand     = -1;
		} else {
			die("Strand unknown\n");
		}
		
		int   tss    = (fr=='+'?atoi(a[4]):atoi(a[5]));
		int   tes    = (fr=='+'?atoi(a[5]):atoi(a[4]));
		char* g      = a[12];
		
		(*a_proms)[ (*numproms) ].str      = line;
		(*a_proms)[ (*numproms) ].tsname    = strdup(n);
		(*a_proms)[ (*numproms) ].chr       = strdup(c);
		(*a_proms)[ (*numproms) ].genename  = strdup(g);
		
		if (region == 0) {
			
			// prom
			if (fr == '+') {
				(*a_proms)[ (*numproms) ].i        = tss-up;
				(*a_proms)[ (*numproms) ].j        = tss+dn;
			} else {
				(*a_proms)[ (*numproms) ].i        = tss-dn;
				(*a_proms)[ (*numproms) ].j        = tss+up;
			}
			
		} else if (region == 1) {
			
			// TES
			if (fr == '+') {
				(*a_proms)[ (*numproms) ].i        = tes-up;
				(*a_proms)[ (*numproms) ].j        = tes+dn;
			} else {
				(*a_proms)[ (*numproms) ].i        = tes-dn;
				(*a_proms)[ (*numproms) ].j        = tes+up;
			}
			
		} else {
			die("Region unrecognized\n");
		}
		
		// frame
		(*a_proms)[ (*numproms)].strand        = strand;
		
		(*numproms) ++;
        
		free(a);
		
	}
	
	
	/* cleanup */
	free(buff);
	fclose(f);
}
