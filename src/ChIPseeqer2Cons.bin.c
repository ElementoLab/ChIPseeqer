#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>

#include "dataio.h"
#include "statistics.h"
#include "sequences.h"
#include "postscript.h"

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

#define  NBPEAKS 40000
#define  NBWINDOWS 1000

typedef struct _GenInt {
	int i;
	int j;
	char* c;
	char* id;
	float score;
	char* str;
} GenInt;

void readUZcons(char* file, float** cons);
void readGZcons(char* file, float** cons, int chrlen);
void readNucleosomePositions(char* file, float** cons, char* chrname, int format, int verbose);
void readMapScores(char* file, float** cons);

int main(int argc, char** argv) {
	
	// general
	long 	i		= 0;
	long	chrlen		= 0;
	int	numchroms	= 24;			//number of chromosomes
	char*	chroms[]	= {"chrY",		//array that holds the chromosomes names
		"chrX",
		"chr9",
		"chr8",
		"chr7",
		"chr6",
		"chr5",
		"chr4",
		"chr3",
		"chr22",
		"chr21",
		"chr20",
		"chr2",
		"chr19",
		"chr18",
		"chr17",
		"chr16",
		"chr15",
		"chr14",
		"chr13",
		"chr12",
		"chr11",
		"chr10",
		"chr1"};
	
	int	verbose		= 0;
	char*	chrname		= NULL;
	int	format		= 0;
	char* 	formattext	= NULL;
	FILE*	f1		= NULL;
	GenInt*	a_int1		= NULL;
	GenInt*	a_int2		= NULL;
	
	int	maxnumint	= 10000000;
	int	mynmax		= 10000;
	int	hasid		= 0;
	char*	intervals	= 0;
	char*	consdir		= 0;
	char*	category	= 0;
	char*	consfile	= 0;
	char*	method		= 0;
	int	print_rand	= 0;
	int     j		= 0;
	int	c		= 0;
	float*  cons		= NULL;
	char*   buff		= NULL;
	char**	p		= NULL;
	int	m		= 0;
	int     numint		= 0;	
	int	show_profiles	= 0;
	
	int	numpos1		= 0;
	float	avgcons1	= 0;
	float	mincons1	= 1000;
	float	maxcons1	= -1000;
	int	numpos2		= 0;
	float	avgcons2	= 0;
	float	mincons2	= 1000;
	float	maxcons2	= -1000;
	int     ws		= 10;
	int     k		= 0;
	int     has_rand	= 0;
	int     make_rand	= 0;
	int     randist		= 50000;
	char*   outfile		= NULL;
	char*   outrandom	= NULL;
	int     adjcons		= 0;
	FILE*   fpout		= stdout;
	FILE*   fpran		= stdout;
	char*   chr		= NULL;	
	char*   chrdata		= NULL;
	int     addheader	= 0;
	char**  chrnames	= NULL;
	int*    chrlens		= NULL;
	float* chrmap		= NULL;
	int     avgonly		= 0;
	int showalldata		= 0;
	
	float**	data1		= NULL;
	float**	data2		= NULL;
	char*	outepsmap	= NULL;
	int	nbcond		= 0;
	int	nborfs		= 0;
	int	peak1		= 0;
	int	bin1		= 0;
	int	peak2		= 0;
	int	bin2		= 0;
	char*	xlabel		= NULL;
	char*	ylabel		= NULL;
	char*	legend1		= NULL;
	char*	legend2		= NULL;
	char* genome	= NULL;

	if (argc < 2) {
		die("Usage: ChIPseeqer2Cons.bin -peakfile FILE -consdir FILE -species STR -show_profiles INT -has_rand INT -make_rand INT -randist INT  -showalldata INT  \n");
	}
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-peakfile"))
		intervals = get_parameter(argc, argv, "-peakfile");
	
	if (exist_parameter(argc, argv, "-regions"))
		intervals = get_parameter(argc, argv, "-regions");
	
	if (exist_parameter(argc, argv, "-avgonly"))
		avgonly = atoi(get_parameter(argc, argv, "-avgonly"));
	
	if (exist_parameter(argc, argv, "-showalldata")) 		
	        showalldata = atoi(get_parameter(argc, argv, "-showalldata"));

	if (exist_parameter(argc, argv, "-targets"))
		intervals = get_parameter(argc, argv, "-targets");
	
	if (exist_parameter(argc, argv, "-chr"))
		chr = get_parameter(argc, argv, "-chr");
	
	if (exist_parameter(argc, argv, "-ws"))
		ws = atoi(get_parameter(argc, argv, "-ws"));
	
	if (exist_parameter(argc, argv, "-hasid"))
		hasid = atoi(get_parameter(argc, argv, "-hasid"));
	
	if (exist_parameter(argc, argv, "-adjcons"))
		adjcons = atoi(get_parameter(argc, argv, "-adjcons"));
	
	if (exist_parameter(argc, argv, "-show_profiles"))
		show_profiles = atoi(get_parameter(argc, argv, "-show_profiles"));
	
	if (exist_parameter(argc, argv, "-consdir"))
		consdir  = get_parameter(argc, argv, "-consdir");
	
	if (exist_parameter(argc, argv, "-addheader"))
		addheader  = atoi(get_parameter(argc, argv, "-addheader"));
	
	if (exist_parameter(argc, argv, "-outfile"))
		outfile  = get_parameter(argc, argv, "-outfile");
	
	if (exist_parameter(argc, argv, "-outrandom"))
		outrandom  = get_parameter(argc, argv, "-outrandom");
	
	if (exist_parameter(argc, argv, "-format")) {
		formattext = get_parameter(argc, argv, "-format");
		
		if (strcmp(formattext, "scores") == 0)
			format = 0;
		else if (strcmp(formattext, "gzscores") == 0)
			format = 1;
		else if (strcmp(formattext, "nucleosomes") == 0)
			format = 2;
		else if (strcmp(formattext, "mapscores") == 0)
			format = 3;
		else {
			die("Format unrecognized.\n");
		}
	}
	
	if (exist_parameter(argc, argv, "-category"))
		category  = get_parameter(argc, argv, "-category");
	
	if (exist_parameter(argc, argv, "-method"))
		method  = get_parameter(argc, argv, "-method");
	
	if (exist_parameter(argc, argv, "-species"))
		category  = get_parameter(argc, argv, "-species");
	
	if (exist_parameter(argc, argv, "-verbose"))
		verbose = atoi(get_parameter(argc, argv, "-verbose"));
	
	if (exist_parameter(argc, argv, "-print_rand"))
		print_rand = atoi(get_parameter(argc, argv, "-print_rand"));
	
	// does the data have random peaks too ?
	if (exist_parameter(argc, argv, "-has_rand"))
		has_rand = atoi(get_parameter(argc, argv, "-has_rand"));
	
	// if not, maybe we should make random peaks
	if (exist_parameter(argc, argv, "-make_rand")) {
		make_rand = atoi(get_parameter(argc, argv, "-make_rand"));
		has_rand = 1;
	}
	
	// if not, maybe we should make random peaks
	if (exist_parameter(argc, argv, "-randist"))
		randist = atoi(get_parameter(argc, argv, "-randist"));
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata  = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);
		if (verbose == 1) {
			printf("Found %d chr\n", numchroms);
		}
	}
	
	if (exist_parameter(argc, argv, "-outepsmap")) {
		outepsmap =      get_parameter(argc, argv, "-outepsmap");
	} else {
		outepsmap = 0;
	}
	
	if (exist_parameter(argc, argv, "-xlabel"))
		xlabel  = get_parameter(argc, argv, "-xlabel");
	
	if (exist_parameter(argc, argv, "-ylabel"))
		ylabel  = get_parameter(argc, argv, "-ylabel");
	
	if (exist_parameter(argc, argv, "-legend1"))
		legend1  = get_parameter(argc, argv, "-legend1");
	
	if (exist_parameter(argc, argv, "-legend2"))
		legend2  = get_parameter(argc, argv, "-legend2");
	
	if(xlabel == NULL) {
		xlabel = "Distance to peak summit";
	}
	
	if(ylabel == NULL) {
		ylabel = "Average conservation level";
	}
	
	if(legend1 == NULL) {
		legend1 = "TF peaks";
	}
	
	if(legend2 == NULL) {
		legend2 = "Random peaks";
	}
	
	if (exist_parameter(argc, argv, "-genome"))
		genome = get_parameter(argc, argv, "-genome");
	
	if ((outrandom != 0) && (make_rand == 0) && (has_rand == 0)) {
		die("no random data to write to -outrandom\n");
	}
	
	// allocate a huge chunk of memory
	data1 = (float**)malloc(NBPEAKS * sizeof(float*));
	for (i=0; i<NBPEAKS; i++) {
		data1[i] = (float*)malloc(NBWINDOWS * sizeof(float));
	}

	
	// allocate a huge chunk of memory
	data2 = (float**)malloc(NBPEAKS * sizeof(float*));
	for (i=0; i<NBPEAKS; i++) {
		data2[i] = (float*)malloc(NBWINDOWS * sizeof(float));
	}
	
	//
	// start reading intervals
	//
	
	a_int1 = (GenInt*)calloc(maxnumint, sizeof(GenInt));
	if (a_int1 == 0) {
		die("Problem allocating a_int1.\n");
	}
	
	a_int2 = (GenInt*)calloc(maxnumint, sizeof(GenInt));
	if (a_int2 == 0) {
		die("Problem allocating a_int2.\n");
	}
	
	// read set of intervals
	numint = 0;
	f1 = fopen(intervals, "r");
	if (!f1) {
		printf("Cannot open %s\n", intervals);
		exit(1);
	}
	
	buff  = (char*)malloc(mynmax * sizeof(char));
	
	while (!feof(f1)) {
		fgets(buff, mynmax, f1);
		if (feof(f1))
			break; 
		chomp(buff);
		char* line = strdup(buff);
		split_line_delim(buff, "\t", &p, &m);
		
		int idx = 0;
		if (hasid == 1) {
			idx = 1;
			a_int1[numint].id = strdup(p[0]);
		} else {
			a_int1[numint].id = NULL;      
		}
		a_int1[numint].c = strdup(p[idx]);
		a_int1[numint].i = atoi(p[idx+1]);
		a_int1[numint].j = atoi(p[idx+2]);
		a_int1[numint].score = 0; 
		a_int1[numint].str = line;

		if ((has_rand == 1) && (make_rand == 0)) { // if has random interval and not asked to make new ones, read them
			a_int2[numint].c = strdup(p[idx+3]);
			a_int2[numint].i = atoi(p[idx+4]);
			a_int2[numint].j = atoi(p[idx+5]);
			a_int2[numint].score = 0; 
		}
		
		free(p);
		free(line);
		numint++;
		
		if (numint == maxnumint) {
			die("Max number of intervals reached.\n");
		}
		
	}
	fclose(f1);
	//
	// end read intervals
	//
	
	if (outrandom != 0) {
		fpran = fopen(outrandom, "w");
		if (fpran == 0) {
			printf("cannot open %s\n", outrandom);
			exit(1);
		}
	}
	
	if (outfile != 0) {
		fpout = fopen(outfile, "w");
		if (fpout == 0) {
			printf("cannot open %s\n", outfile);
			exit(1);
		}
	}
	
	
	if ((addheader == 1) && (show_profiles == 0)) {
		fprintf(fpout, "chr\tstart\tend\tavg");
		if (avgonly == 0) {
			fprintf(fpout, "\tmin\tmax");
		}
		fprintf(fpout, "\n");		  
	}
	
	
	for (c=0; c<numchroms; c++) {	
		
		if (chrdata == 0) { // hg18 defaults
			chrname = (char*)(chroms[c]);
			chrlen = hg18_chrlen(chrname);
		} else {
			chrname = chrnames[c];
			chrlen  = chrlens[c];
		}
		
		if (verbose == 1) 
			printf("Chr %d = %s\n", c, chrname);
		
		
	    if ((chr != 0) && (strcmp(chr, chrname) != 0))
			continue;
		
		printf("Processing chromosome: %s\n", chrname);
		
	    consfile = (char*)calloc(100000, sizeof(char)); // strlen(consdir) + strlen(chrname) + strlen(category) + 32, sizeof(char));
		
		//
		// init cons vector
		//
		cons = (float*)calloc(chrlen, sizeof(float));
		if (cons == 0) {
			die("Cannot allocate memory for conservation vector.\n");
		}
		
		// important ... conservation scores with no data are set to -1
		for (i=0; i<chrlen; i++) 
			cons[i] = -1;
		
		if (format == 0) {
			if(strcmp(genome, "hg18") == 0)
				sprintf(consfile, "%s/%s.%s44way.%s.wigFix", consdir, chrname, method, category);
			else if (strcmp(genome, "mm9") == 0)
				sprintf(consfile, "%s/%s.%s30way.wigFix", consdir, chrname, method);
			readUZcons(consfile, &cons);
		}
		else if (format == 1) {
			if(strcmp(genome, "hg18") == 0)
				sprintf(consfile, "%s/%s.%s44way.%s.wigFix.gz", consdir, chrname, method, category);
			else if (strcmp(genome, "mm9") == 0){
				sprintf(consfile, "%s/%s.%s30way.wigFix.gz", consdir, chrname, method);
			}
			readGZcons(consfile, &cons, chrlen);
		}
		else if (format == 2) {
			sprintf(consfile, "%s/RestingNucleosomes-%s-scores.vstep", consdir, chrname);
			readNucleosomePositions(consfile, &cons, chrname, 1, verbose);
		}
		else if (format == 3) {
			sprintf(consfile, "%s/map.%s", consdir, chrname);
			readMapScores(consfile, &cons);
		}
		
		// now go through all intervals
		// printf("chr-start-end\tavgcons\tmincons\tmaxcons\n");
		
		peak1 = 0;
		peak2 = 0;
		
		for (i=0; i<numint; i++) {
			
			// if not right chr, skip
			if ( strcmp(a_int1[i].c, chrname) != 0 ) 
				continue;
			
			if (make_rand == 1) {
			    // get random segments on the spot
			    a_int2[i].c = a_int1[i].c;
			    a_int2[i].i = rand_chr_interval_startpoint_around_interval(a_int1[i].c, a_int1[i].i, a_int1[i].j, randist);
			    a_int2[i].j = a_int2[i].i + ( a_int1[i].j - a_int1[i].i + 1); // just add len			    
			}
			
			//
			// MAIN OUTPUT: interval, avg cns, min, max, etc (that is, not a profile)
			//			
			if (show_profiles == 0) {
				
				numpos1  = 0;
				avgcons1 = 0;
				mincons1 = 1000;
				maxcons1 = -1000;			  
				numpos2  = 0;
				avgcons2 = 0;
				mincons2 = 1000;
				maxcons2 = -1000;
				
				// compute avg conservation for peak intervals
				for (j=a_int1[i].i; j<=a_int1[i].j; j++) {
					if (cons[j] >= 0) {
						avgcons1 += cons[j];
						numpos1++;
						if (cons[j] > maxcons1)
							maxcons1 = cons[j];
						
						if (cons[j] < mincons1)
							mincons1 = cons[j];
					}
				}
				
				if (has_rand == 1) {
					// compute avg conservation for random peak intervals			
					for (j=a_int2[i].i; j<=a_int2[i].j; j++) {
						if (cons[j] >= 0) {
							avgcons2 += cons[j];
							numpos2++;
							if (cons[j] > maxcons2)
								maxcons2 = cons[j];
							
							if (cons[j] < mincons2)
								mincons2 = cons[j];
						}
					}			    
				}
				
				// print results
				if (showalldata == 1) // show entire line		
				  fprintf(fpout, "%s", a_int1[i].str); 		
				else 
				  fprintf(fpout, "%s\t%d\t%d", a_int1[i].c, a_int1[i].i, a_int1[i].j);			  

				if (numpos1 > 0) { 
					avgcons1 /= numpos1;
					fprintf(fpout, "\t%4.3f", avgcons1);
				}
				else {
					fprintf(fpout, "\tNA");
				}			 
				if (avgonly == 0) {
					fprintf(fpout, "\t%4.3f\t%4.3f", mincons1, maxcons1);
				}
				
				
				if ((make_rand == 1) && (fpran != 0)) {
					fprintf(fpran, "%s\t%d\t%d", a_int2[i].c, a_int2[i].i, a_int2[i].j);			
					if (numpos2 > 0) { 
						avgcons2 /= numpos2;
						fprintf(fpran, "\t%4.3f", avgcons2);
					}
					else {
						fprintf(fpran, "\tNA");
					}
					
					if (avgonly == 0) {
						fprintf(fpran, "\t%4.3f\t%4.3f", mincons2, maxcons2);
					}
					
					fprintf(fpran, "\n");
					
				} else if (adjcons == 1) {
					// NEW ADDITION: SHOW CONS iN ADJACENT REGIONS
					
					// get left segment
					int ls = a_int1[i].j - a_int1[i].i + 1;
					free(a_int2[i].c);
					a_int2[i].c = a_int1[i].c;
					a_int2[i].i = a_int1[i].i - ls;
					a_int2[i].j = a_int1[i].j - ls;
					avgcons2 = 0;
					numpos2  = 0;
					mincons2 = 1000;
					maxcons2 = -1000;
					for (j=a_int2[i].i; j<=a_int2[i].j; j++) {
						if (cons[j] >= 0) {
							avgcons2 += cons[j];
							numpos2++;
							if (cons[j] > maxcons2)
								maxcons2 = cons[j];
							
							if (cons[j] < mincons2)
								mincons2 = cons[j];
						}
					}	
					if (numpos2 > 0) { 
						avgcons2 /= numpos2;
						fprintf(fpout, "\t%4.3f", avgcons2);
					} else {
						fprintf(fpout, "\tNA");
					}
					
					// get right segment
					a_int2[i].c = a_int1[i].c;
					a_int2[i].i = a_int1[i].i + ls;
					a_int2[i].j = a_int1[i].j + ls;
					float avgcons3 =     0;
					int numpos3    =     0;
					float mincons3 =  1000;
					float maxcons3 = -1000;
					for (j=a_int2[i].i; j<=a_int2[i].j; j++) {
						if (cons[j] >= 0) {
							avgcons3 += cons[j];
							numpos3++;
							if (cons[j] > maxcons3)
								maxcons3 = cons[j];				
							if (cons[j] < mincons3)
								mincons3 = cons[j];
						}
					}	
					if (numpos3 > 0) { 
						avgcons3 /= numpos3;
						fprintf(fpout, "\t%4.3f", avgcons3);
					} else {
						fprintf(fpout, "\tNA");
					}
					
					if ((numpos1 > 0) && (numpos2 > 0) && (numpos3 > 0)) {
						fprintf(fpout, "\t%4.3f", 2.0*avgcons1/(avgcons2+avgcons3+0.002));
					} else {
						fprintf(fpout, "\tNA");
					}
					
				} // if adcons == 1
				
				// end of line in any case
				fprintf(fpout, "\n");				
				
			} else {
				
				//
				// ALTERNATIVE OUTPUT: PROFILES
				//
				
				// smoothe then show profile
				fprintf(fpout, "%s-%d-%d", a_int1[i].c, a_int1[i].i, a_int1[i].j);
				
				bin1 = 0;
				
				for (j=a_int1[i].i; j<=a_int1[i].j; j+=ws) {
					avgcons1 = 0;
					numpos1  = 0;
					for (k=j; k<min(chrlen,j+ws); k++) {
						if (cons[k] >= 0) {
							avgcons1 += cons[k];
							numpos1++;
						}
					}
					if (numpos1 == 0) 
						avgcons1 = 0.0;
					else 
						avgcons1 /= numpos1;
					
					fprintf(fpout, "\t%3.1f", avgcons1);
					
					data1[peak1][bin1] = avgcons1;
					
					bin1++;
				}
				fprintf(fpout, "\n");
				
				if ((make_rand == 1) && (fpran != 0)) {
					// same for random peaks
					fprintf(fpran, "%s-%d-%d", a_int2[i].c, a_int2[i].i, a_int2[i].j);
					
					bin2 = 0;

					for (j=a_int2[i].i; j<=a_int2[i].j; j+=ws) {
						avgcons2 = 0;
						numpos2  = 0;
						for (k=j; k<min(chrlen,j+ws); k++) {
							if (cons[k] >= 0) {
								avgcons2 += cons[k];
								numpos2++;
							}
						}
						if (numpos2 == 0)
							avgcons2 = 0;
						else 
							avgcons2 /= numpos2;			    
						
						fprintf(fpran, "\t%3.1f", avgcons2);
						
						data2[peak2][bin2] = avgcons2;

						bin2++;
					}
					fprintf(fpran, "\n");
				} // if need to putput ran
				
			}
		
			peak1++;
			peak2++;
			
		} // end for loop over intervals
		
		free(consfile);
		free(cons);
	} //loop over chromosomes
	
	nborfs = peak1;
	nbcond = bin1;
	
	//printf("\n%d\t%d\n", peak, bin);
	
	if (outepsmap != 0) {
		printf("Writing EPS map ...");
		colMeans2Ddoubleplot(outepsmap, data1, data2, nborfs, nbcond, xlabel, ylabel, ws, legend1, legend2); 
		printf("Done.\n");
	}	
	
	
	/* cleanup */
	for (i=0; i<NBPEAKS; i++)
		free(data1[i]);
	for (i=0; i<NBPEAKS; i++)
		free(data2[i]);
	free(data1);
	free(data2);
	for (i=0; i<numint; i++) {
		if (hasid == 1)
			free(a_int1[i].id);
		if (a_int1[i].c == a_int2[i].c)
			a_int2[i].c = NULL;
		free(a_int1[i].c);
	}
	for (i=0; i<numint; i++) {
		if (a_int2[i].c != NULL)
			free(a_int2[i].c);
	}
	free(a_int1);
	free(a_int2);
	free(buff);
	if (outrandom != 0)
		fclose(fpran);
	if (outfile != 0) 
		fclose(fpout);
	
	return 0;
}

//
// start reading gzipped conservation file
//
void readGZcons(char* file, float** cons, int chrlen)
{
	
	gzFile zf;
	char* buff	= NULL;
	int   len	= 1000;
	int   j		= 0;
	int   i		= -1;
	
	buff = (char*)malloc(1024 * sizeof(char));
	
	//char * gzgets  (gzFile file, char * buf, int len);
	zf = gzopen (file , "rb" );
	while (!gzeof(zf)) {
		
		gzgets(zf, buff, len);
		
		if (gzeof(zf))
			break;
		
		// process this type of line: fixedStep chrom=chr22 start=14430001 step=1
		
		if (buff[0] == 'f') {
			
			int lb = strlen(buff);
			int ce = 0;
			j = 0;
			// advance to start=
			while ((j < lb) && (ce < 2)) {
				//printf("%c\n", buff[j]);
				if (buff[j] == '=')
					ce ++;
				j++;
			}
			
			char pos_start[100];
			int  j_d = 0;
			while (buff[j] != ' ') {
				//printf("%c\n", buff[j]);
				pos_start[j_d] = buff[j];
				j_d ++;
				j   ++;
			}
			pos_start[j_d] = '\0';
			
			//fflush(stdout);      
			//printf("Got pos_start = %s\n", pos_start);
			
			i = atoi(pos_start);
			
			continue;
		}
		
		// other than that, just store conservation score
		
		if(i>chrlen) { 		
			break; 		
		} 		
		else {
			(*cons)[ i ] = atof(buff);
		
			i++;
		}
	}
	
	/* cleanup */
	gzclose(zf);
	free(buff);
}  


//
// start reading gzipped conservation file
//
void readUZcons(char* file, float** cons)
{
	
	FILE* fp;
	char* buff = 0;
	int   len  = 1000;
	int   j;
	int   i = -1;
	
	buff = (char*)malloc(1024 * sizeof(char));
	
	
	//char * gzgets  (gzFile file, char * buf, int len);
	fp = fopen (file , "r" );
	while (!feof(fp)) {
		
		fgets(buff, len, fp);
		
		if (feof(fp))
			break;
		
		// process this type of line: fixedStep chrom=chr22 start=14430001 step=1
		
		if (buff[0] == 'f') {
			
			int lb = strlen(buff);
			int ce = 0;
			j = 0;
			// advance to start=
			while ((j < lb) && (ce < 2)) {
				//printf("%c\n", buff[j]);
				if (buff[j] == '=')
					ce ++;
				j++;
			}
			
			char pos_start[100];
			int  j_d = 0;
			while (buff[j] != ' ') {
				//printf("%c\n", buff[j]);
				pos_start[j_d] = buff[j];
				j_d ++;
				j   ++;
			}
			pos_start[j_d] = '\0';
			
			//fflush(stdout);      
			//printf("Got pos_start = %s\n", pos_start);
			
			i = atoi(pos_start);
			
			continue;
		}
		
		// other than that, just store conservation score
		(*cons)[ i ] = atof(buff);
		
		i++;
	}
	
	/* cleanup */
	fclose(fp);
	free(buff);
}  

void readNucleosomePositions(char* file, float** cons, char* chrname, int format, int verbose)
{
	FILE* f2;
	char* buff = 0;
	
	int   mynmax = 100000;
	char** p;
	int    m;
	
	//
	// start reading conservation file
	//
	
	float avg = -1;
	int   npo = -1;
	int   st  = -1;
	int   en  = -1;
	int   i;
	
	f2 = fopen(file, "r");
	if (!f2) {
		printf("Cannot open file %s\n", file);
		exit(1);
	}
	
	buff  = (char*)malloc(mynmax * sizeof(char));
	
	if (format == 1) {
		// nucl, remove top lines
		fgets(buff, mynmax, f2);
		fgets(buff, mynmax, f2);    
	}
    
	
	while (!feof(f2)) {
		
		// read cons fragment
		fgets(buff, mynmax, f2);
		if (feof(f2))
			break; 
		chomp(buff);
		split_line_delim(buff, "\t", &p, &m);
		
		// conservation scores
		if (format == 0) {
			
			if (strcmp(p[1], chrname) != 0)
				die("Chrname does not match.\n");
			
			avg = atof(p[12]);
			npo = atoi(p[11]);
			
			avg /= npo;
			
			if (verbose == 1) 
				printf("%f\t%d\n", avg, npo);
			
			st  = atoi(p[2]);
			en  = atoi(p[3]);
			
			for (i=st; i<=en; i++) {
				(*cons)[i] = avg;
			}
			
		} else if (format == 1) {
			// nucleosome scores
			
			st  = atoi(p[0]);
			avg = atof(p[1]);
			for (i=st; i<=st+9; i++) {
				(*cons)[i] = avg;
			}
		}
		
		free(p);
	}
	
	
	/* cleanup */
	//flose(f2);
	//free(buff);
	
	if (verbose == 1)
		printf("read all nucleosome positions.\n");
}

void readMapScores(char* file, float** cons)
{
	
	FILE* f2;
	char* buff = 0;
	
	int   mynmax = 100000;
	
	int   i;
	
	f2 = fopen(file, "r");
	if (!f2) {
		printf("Cannot open file %s\n", file);
		exit(1);
	}
	
	buff  = (char*)malloc(mynmax * sizeof(char));
	
	i = 0;
	while (fgets(buff, mynmax, f2) != 0) {
		chomp(buff);
		(*cons)[i] = atof(buff);	  
		i++;
	}
	
	/* cleanup */
	//flose(f2);
	//free(buff);
}
