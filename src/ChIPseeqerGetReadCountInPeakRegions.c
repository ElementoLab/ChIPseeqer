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
void readRefGenePromotersOrDownstream(char* file, GenInt** a_proms, int* numproms, int up, int dn, int region);


int main(int argc, char** argv) {
	
	// general
	long i;
	
	// reading file
	char* buff				= NULL;
	int   mynmax				= 100000;
	char** p				= NULL;
	int    m				= 0;
	
	unsigned short*	chip_counts		= NULL;
	unsigned short*	input_counts		= NULL;
	
	char*		chipfile		= 0; 
	char*		inputfile		= 0; 
	long		chrlen			= 0; 
	
	int		fraglen                 = -1;
	int		readlen                 = 0;
	int		numreads_chip           = 0;
	int		numreadsclonal_chip     = 0;
	int		numreads_input          = 0;
	int		numreadsclonal_input    = 0;	
	int		t			= 5;
	double		fold_t			= 2.0;
	
	int		verbose			= 0;
	char*		chrname			= 0;
	double		mappability		= 1.0;
	int		fast			= 1;
	int		format			= 0;
	char*		formattext		= 0;
	
	FILE*   	f1 = NULL;
	GenInt* 	a_int;
	int     	maxnumint		= 10000000;
	int     	numint			= 0;
	int     	hasid			= 0;
	char*   	intervals		= 0;
	int     	ws			= 100;
	int     	j			= 0;
	int	  	k				= 0;
	char* 		chrdata			= 0;
	char** 		chrnames		= 0;
	int* 		chrlens        		= 0;
	int 		numchroms		= 24;
	char*  		chipdir			= 0;
	
	char readfile[1000];
	int  		normalize     		= 1;
	double 		normfactor  		= 0;
	int  		numreads      		= 0;
	int  		totalnumreads 		= 0;
	int  		normto        		= 10000000;
	int  		c					= 0;
	int  		genome        		= 0;
	char* 		inputdir     		= 0;
	int   		ext          		= 0;
	int  		st, en;
	unsigned short int* tmp_rat = NULL;
	char* 		output_text		= 0;
	int 		output			= 0; //avg (max=1)
	float 		out				= 0.0;
	float 		sum				= 0.0;
	char* 		mystr			= NULL;
	int 		idx				= 0.0;
	float 		tmpval			= 0.0;
	int   		motifmatches		= 0;
	int   		randomize		= 0;
	int   		seed			= 1000;
	float* 		chrmap			= 0;
	
	seqI     	si;
	char*    	seq			= NULL;
	int      	size		= 0;
	char*    	name		= NULL;
	char*    	fastadir		= 0;
	int      	checkN			= 0;
	char     	fastafile[1000];
	int      	uniquereads 		= 1;
	int      	refgene			= 0;
	int     	up			= 2000;
	int      	dn			= 2000;
	char*    	outfile			= 0;
	FILE*    	fpout			= stdout;
	Histogram 	h;
	
	char*    	mapdir			= 0;
	unsigned char* 	mapmap;
	char     	mapfile[1000];
	int      	mapt			= -1;
	int      	tsregion   		= 0;
	char*           desc = 0;
	if (argc < 2) {
		die("Usage: ChIPseeqerGetReadCountInPeakRegions [ -refgene FILE -up INT -dn INT -tsregion 0/1 | -intervals FILE ] -chipdir DIR [ -fraglen INT(-1 def) -output max|avg(def) ] -chrdata FILE  -desc STR \n");
	}
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-outfile"))
		outfile = get_parameter(argc, argv, "-outfile");
	
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
	
	if (exist_parameter(argc, argv, "-hasid"))
		hasid = atoi(get_parameter(argc, argv, "-hasid"));
	
	if (exist_parameter(argc, argv, "-ws"))
		ws = atoi(get_parameter(argc, argv, "-ws"));
	
	if (exist_parameter(argc, argv, "-chipfile"))
		chipfile  = get_parameter(argc, argv, "-chipfile");
	
	if (exist_parameter(argc, argv, "-chipdir"))
		chipdir  = get_parameter(argc, argv, "-chipdir");
	
	if (exist_parameter(argc, argv, "-mapdir"))
		mapdir  = get_parameter(argc, argv, "-mapdir");  
	
	if (exist_parameter(argc, argv, "-mapt"))
		mapt  = atoi(get_parameter(argc, argv, "-mapt"));  
	
	if (exist_parameter(argc, argv, "-inputdir"))
		inputdir  = get_parameter(argc, argv, "-inputdir");
	
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
		
	if (exist_parameter(argc, argv, "-desc")) // desc
	  desc = get_parameter(argc, argv, "-desc");
	else {
	  die("You need to specify -desc\n");
	}

	if (exist_parameter(argc, argv, "-format")) {
	  formattext = get_parameter(argc, argv, "-format");
	  format     = formatToId(formattext);				
	}
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata  = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);
		genome   = 0;
	} else {
		die("Please define -chrdata\n");
	}
	
	if (exist_parameter(argc, argv, "-seed")) {
		seed = atoi(get_parameter(argc, argv, "-seed"));  
	}
	default_set_seed(seed);
	
	if (normalize == 1) {
		fprintf(stderr, "# normalize to %d reads\n", normto);
	}
	
	fprintf(stderr, "# tsregion = %d\n", tsregion);
	fprintf(stderr, "# fraglen = %d\n", fraglen);
	fprintf(stderr, "# uniquereads = %d\n", uniquereads);
	fprintf(stderr, "# output = %d\n", output);

	if (refgene == 1) {
		
	    readRefGenePromotersOrDownstream(intervals, &a_int,  &numint, up, dn, tsregion);
		
	} else {
		//
		// start reading intervals
		//
		maxnumint = nbLinesInFile(intervals);
		a_int = (GenInt*)malloc(maxnumint * sizeof(GenInt));
		if (a_int == 0) {
			printf("Problem allocating a_int(%d).\n", maxnumint);
			exit(1);
		}
		
		buff  = (char*)malloc(mynmax * sizeof(char));
		
		// read first set of intervals
		numint = 0;
		f1 = fopen(intervals, "r");
		if (!f1) {
			die("Cannot open f1\n");
		}
		while (fgets(buff, mynmax, f1) != 0) {
			chomp(buff);
			mystr = strdup(buff);
			split_line_delim(buff, "\t", &p, &m);
			
			idx = 0;
			if (hasid == 1) {
				idx = 1;
				a_int[numint].tsname = strdup(p[0]);
			} else {
				a_int[numint].tsname = 0;      
			}
			
			a_int[numint].chr = strdup(p[idx]);
			a_int[numint].i = atoi(p[idx+1]);
			
			if (motifmatches == 1)
				a_int[numint].j =  a_int[numint].i + strlen(p[4]);
			else 
				a_int[numint].j = atoi(p[idx+2]);
			
			a_int[numint].score = 0; 
			a_int[numint].str = mystr;
			free(p);
			numint++;
		}
	}

	if (verbose == 1)
		fprintf(stderr, "# Found %d intervals to evaluate.\n", numint);
	
	//
	// end read intervals
	//
	
	// calc total number of reads for RPKM-style 
	if (normalize == 1) {
		totalnumreads = 0;
		for (c=0; c<numchroms; c++) {
			chrname = chrnames[c];
			chrlen  = chrlens[c];

			sprintf(readfile, "%s/reads.%s", chipdir, chrname); 
			
			if (!file_exists(readfile))
				continue;
			
			if (uniquereads == 0) 
			  numreads = CountReads(readfile, format);
			else 
			  numreads = CountNonClonalReads(readfile, format, chrlen);

			if (verbose == 1) {
			  fprintf(stderr, "# Found %d aligned reads in chr %s (uniquereads = %d)\n", numreads, chrname, uniquereads);
			}
			totalnumreads += numreads;
			
		}
		if (verbose == 1) {
			fprintf(stderr, "# Found a total of %d aligned reads\n", totalnumreads);
		}
		
		normfactor = normto / (double)totalnumreads;
		if (verbose == 1) {
			fprintf(stderr, "# Normalization factor = %4.3f\n", normfactor);
		}
		
	}
	
	if (outfile != 0) {
		fpout = fopen(outfile, "w");
		if (fpout == 0) {
			printf("Cannot open %s for saving\n", outfile);
			exit(0);
		}
		
		Histogram_create(&h, 100, 100);
	}
	
	
	fprintf(fpout, "SEQID\t%s\n", desc);
	
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
		chip_counts = calloc(chrlen, sizeof(unsigned short));
		if (chip_counts == 0) {	    
			printf("Problem in allocating chip_counts (chrlen=%ld).\n", chrlen);
			exit(1);
		}
		
		if (inputdir != 0) {
			input_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
			if (input_counts == 0) {
				die("Problem in allocating input_counts.\n");
			}
		}
		
		sprintf(readfile, "%s/reads.%s", chipdir, chrname); 
		
		if (!file_exists(readfile)) {
			/* cleanup */
			free(chip_counts);
			if (inputdir != 0)
				free(input_counts);
			continue;
		}
		
		getCountFromReadFile(readfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &chip_counts, &numreads_chip, &numreadsclonal_chip);
		
		if (mapdir != 0) {
			sprintf(mapfile, "%s/map.%s", mapdir, chrname);
			mapmap = create_binarized_array(chrlen);
			readMapData(mapfile, &mapmap, mapt);
		}
		
		if (outfile != 0) {
			
			for (i=0; i<chrlen; i++) {
				tmpval = (float)(chip_counts[i]);
				if (normalize == 1)
					tmpval *= normfactor;	
				
				if ((mapdir == 0) || ((mapdir != 0) && (get_entry_in_binarized_array(mapmap, i) > 0))) {
					Histogram_add(&h, tmpval);
				}	
			}
		}
		
		if (mapdir != 0)
			free(mapmap);
		
		if (inputdir != 0) {
			sprintf(readfile, "%s/reads.%s", inputdir, chrname); 
			getCountFromReadFile(readfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &input_counts, &numreads_input, &numreadsclonal_input);
		}
		
		if (verbose == 1)
		  fprintf(stderr, "# Number of reads = %d (chip) and %d (input)\n", numreads_chip, numreads_input);
		
		if ((randomize == 1) && (checkN != 0)) {
			sprintf(fastafile, "%s/%s.fa", fastadir, chrname); 
			
			// working for the human genome (fast)
			seqI_set_max_seqlen(&si, 500000000UL);
			seqI_set_seqlen_inc(&si, 300000000UL);
			seqI_set_fread_chunk(&si, 100000000UL);
			
			if (seqI_open_fast(&si, fastafile) == 0) {
				die("Error opening file\n");
			}
			
			seq = seqI_nextSequence_fast(&si, &name, &size);
			if (verbose == 1)
				fprintf(stderr, "# Loaded %s (size=%d)\n", fastafile, size);
			seqI_close(&si);
		}
		
		
		// now go thru all intervals
		for (i=0; i<numint; i++) {
			
			if ( strcmp(a_int[i].chr, chrname) != 0 )
				continue;
			
			// extend
			st = max(0, a_int[i].i - ext);
			en = min(chrlen-1, a_int[i].j + ext);
			
			// interval length
			int li = en - st + 1;
			
			if (randomize == 1) {
				if (checkN == 0) {
					st = (int)(0.5 + default_rand() * chrlen - li);
					en = st + li;
				} 
				else {
					while (1) {
						st = (int)(0.5 + default_rand() * chrlen - li);
						en = st + li;
						char* ss = substr(seq, st, li);
						char* sc = uc(ss);
						int cntN = 0;
						for (j=0; j<li; j++) {
							if (sc[j] == 'N') 
								cntN ++;	       
						}	
						
						if (verbose == 2)
							fprintf(stderr, "%s\n", sc);
						
						free(ss);
						free(sc);
						if ( (cntN / (float)li) < 0.25) {
							break;
						} 
					} // while (1) 
				} // if checkN
			} // if randomize
			
			// allocate temporary interval
			tmp_rat = (unsigned short int*)calloc(li, sizeof(unsigned short int));
			
			// go thru interval
			for (j=st,k=0; j<=en; j++,k++) {
				
				if (inputfile != 0)
					tmp_rat[k] = max(0, chip_counts[j] - input_counts[j]) ;
				else 
					tmp_rat[k] = chip_counts[j] ;
			}
			
			if (output == 1) { // compute maximum
				out = -10000;
				for (j=0; j<li; j++) {
					tmpval = (float)(tmp_rat[j]);
					if (verbose == 2)
						fprintf(stderr, "%d\n", tmp_rat[j]);
					if (normalize == 1)
						tmpval *= normfactor;	  
					if (tmpval > out) {
						out = tmpval;
					}
				}       
			} else { //compute average
				
				sum = 0;
				for (j=0; j<li; j++) {
					tmpval = (float)(tmp_rat[j]);
					if (verbose == 2)
						fprintf(stderr, "%d\n", tmp_rat[j]);
					if (normalize == 1)
						tmpval *= normfactor;	  
					sum += tmpval;
				}
				sum /= li;
				out  = sum;
			}
			
			free(tmp_rat);
			
			if (refgene == 1) {
				fprintf(fpout, "%s\t", a_int[i].tsname);
			} else {
				fprintf(fpout, "%s\t", a_int[i].str);
			}
			
			fprintf(fpout, "%4.3f\n", out);
		}
		
		/* cleanup */
		if ((randomize == 1) && (checkN != 0))
			free(seq);
		free(chip_counts);
		if (inputdir != 0)
			free(input_counts);
	}  // loop over chr
	
	
	if (outfile != 0) {
		
		char prmfile[1000];
		sprintf(prmfile, "%s.prm", outfile);
		FILE* fprm = fopen(prmfile, "w");
		if (!fprm) { die("Cannot open parameter file\n"); }
		fprintf(fprm, "q0.25\t%f\n", Histogram_quantile(&h, 0.25));
		fprintf(fprm, "q0.50\t%f\n", Histogram_quantile(&h, 0.50));
		fprintf(fprm, "q0.75\t%f\n", Histogram_quantile(&h, 0.75));
		fclose(fprm);
		fclose(fpout);
		
		//Histogram_print(&h);
	}
	
	/* cleanup */
	fclose(f1);	
	for (i=0; i<numint; i++) {
		if (hasid == 1) 
			free(a_int[1].tsname);
		free(a_int[i].chr);
		free(a_int[i].str);
	}
	free(a_int);
	free(buff);
	for (c=0; c<numchroms; c++)
		free(chrnames[c]);
	free(chrnames);
	free(chrlens);
	free(chrmap);
	
	return 0;	
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
		char  fr     = a[3][0];
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
		(*a_proms)[ (*numproms)].fr        = (fr=='+'?1:-1);
		
		(*numproms) ++;
        
		free(a);
		
	}
	free(buff);
	fclose(f);
}
