//
// This program returns average normalized read density data for a given set of genes
// prom, first exon, first intron, second exon, second intron, exons, introns, downtream
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
#include <math.h>


#define MAX 0
#define RPKM 1

//#include "protos.h"
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
	int** exons;
	int   numexons;
	
} GenInt;

void readFileSumUpColumnPerChrom(char* file, int col, int colchr, int** counts, int** countgenes);
void readRefGenePromoters(char* file, GenInt** a_proms, int* numproms, int up, int dn);

void QuantizeReadCountsOverInterval(int i, int j, int q, unsigned short* counts, double** qavgcounts);
void readGenesAndExonDataFromRefGene(char* file, GenInt** a_genes, int* numgenes);
void readGeneListIntoHash(char* file, struct my_hsearch_data** hash_genes);

int verbose = 0;

int main(int argc, char** argv) {
	
	// general
	long i;
	
	unsigned short*	chip_counts;
	//unsigned short*	input_counts;
	
	long			chrlen			= 0; 	
	int				fraglen         = 170;
	int				readlen         = 0;
	int				numreads_chip   = 0;
	int				numreads_input  = 0;
	int				t				= 5;
	char*			chrname			= 0;
	double			mappability		= 1.0;
	int				format			= 0;
	char*			formattext		= 0;
	
	char*			intervals		= 0;
	//char*			targets			= 0;

	int				j;
	char*           chrdata         = 0;
	char**          chrnames        = 0;
	int*			chrlens         = 0;
	float*			chrmap			= 0;
	int				numchroms		= 24;
	char*			chipdir			= 0;
	
	char			readfile[1000];
	int			normalize		= MAX;

	double			normfactor		= 0;
	int				numreads		= 0;
	int				totalnumreads	= 0;
	int				normto			= 10000000;
	int				c;
	int				genome			= 0;
	char*			inputdir		= 0;
	int				ext				= 0;
	char*			output_text		= 0;
	int				output			= 0;
	
	int				uniquereads		= 1;
	int				q				= 5; 
	double*			qavgcounts		= 0;
	GenInt*			a_genes;
	int				numgenes;
	int				nq;
	int				promlen			= 3000;
	int				efnq;
	char*			onlychr			        = 0;
	char*			onlyts			        = 0;
	int				posexon			= 0;
	char*			genelist		        = 0;
	char*			labelsALL[8]		        = {"Promoter", "Exon 1", "Intron 1", "Exon 2", "Intron 2", "Other Exons", "Other Introns", "Downstream"};
	char*                   labelsBO[3]                     = {"Upstream", "Gene Body", "Downtream"};
	struct			my_hsearch_data* hash_genes     = 0;
	ENTRY			e;
	ENTRY*			ep                              = 0;
	int			le=0, li=0;
	int			numbins                         = 10;  // will be used for geneparts bodies (q ignored)
	int			numreadsclonal_chip          	= 0;
	//int			numreadsclonal_input	        = 0;
	char*                   geneparts                       = "all";
	int                     numgp                           = 0;
	char*                   outdata                         = 0;
	FILE*                   fpod                            = 0;
	char*                   lab                             = 0;

	if (argc < 2) {
		die("Usage: ChIPseeqerMakeGeneEnrichmentData -intervals FILE -chipdir DIR [ -geneparts [all|body] -numbins INT ]\n");
	}
	
	if (exist_parameter(argc, argv, "-intervals"))
		intervals = get_parameter(argc, argv, "-intervals");
	
	if (exist_parameter(argc, argv, "-genelist"))
		genelist = get_parameter(argc, argv, "-genelist");
	
	if (exist_parameter(argc, argv, "-onlychr"))
		onlychr = get_parameter(argc, argv, "-onlychr");
	
	if (exist_parameter(argc, argv, "-onlyts"))
		onlyts = get_parameter(argc, argv, "-onlyts");
	
	if (exist_parameter(argc, argv, "-q"))
		q = atoi(get_parameter(argc, argv, "-q"));
	
	if (exist_parameter(argc, argv, "-chipdir"))
		chipdir  = get_parameter(argc, argv, "-chipdir");
	
	if (exist_parameter(argc, argv, "-inputdir"))
		inputdir  = get_parameter(argc, argv, "-inputdir");
	
	if (exist_parameter(argc, argv, "-promlen"))
		promlen  = atoi(get_parameter(argc, argv, "-promlen"));
	
	if (exist_parameter(argc, argv, "-chrname"))
		chrname  = get_parameter(argc, argv, "-chrname");
	
	// type of plot
	if (exist_parameter(argc, argv, "-geneparts"))
		geneparts  = get_parameter(argc, argv, "-geneparts");
	if (exist_parameter(argc, argv, "-numbins"))
	  numbins  = atoi(get_parameter(argc, argv, "-numbins"));

	

	if (exist_parameter(argc, argv, "-uniquereads"))
		uniquereads  = atoi(get_parameter(argc, argv, "-uniquereads"));
	
	if (exist_parameter(argc, argv, "-t"))
		t  = atoi(get_parameter(argc, argv, "-t"));
	
	if (exist_parameter(argc, argv, "-normalize")) {
	  char* normalizetxt  = get_parameter(argc, argv, "-normalize");
	  if (strcmp(normalizetxt, "rpkm") == 0)
	    normalize = RPKM;
	  else if (strcmp(normalizetxt, "max") == 0)
	    normalize = MAX;
	  else 
	    die("Don't know this norm scheme\n");
	}
	
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
	
	if (exist_parameter(argc, argv, "-format")) {
		formattext = get_parameter(argc, argv, "-format");
		if (strcmp(formattext, "mit") == 0)
			format = 1;
		else if (strcmp(formattext, "bed") == 0)
			format = 2;
		else if (strcmp(formattext, "sam") == 0)
			format = 3;
		else if (strcmp(formattext, "bowtiesam") == 0)
			format = 4;
	}
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata  = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);
		genome   = 0;
	}

	// where to store the profiles (1 per gene)
	if (exist_parameter(argc, argv, "-outdata")) {
	  outdata  = get_parameter(argc, argv, "-outdata");
	  fpod     = fopen(outdata, "w");
	  if (!fpod) {
	    printf("cannot open %s\n", outdata);
	    exit(1);
	  }
	}	
	if (verbose == 1) {
		printf("About to start reading genes \n");
	}
	
	// read gene data from annotation (does it really work with non-refgene?) 
	readGenesAndExonDataFromRefGene(intervals, &a_genes, &numgenes);
	
	if (verbose == 1) {
		printf("About to start reading genes into hash \n");
	}
	
	//  read genelist
	if (genelist != 0) {
		readGeneListIntoHash(genelist, &hash_genes);
	}
	
	// let's estimate the number of parts
	if (strcmp(geneparts, "all") == 0)
	  numgp = 8;
	else if (strcmp(geneparts, "body") == 0)
	  numgp = 3;
	else {
	  die("Unknown plot type\n");
	}

	if (normalize == RPKM)
	  fprintf(stderr, "# norm=RPKM\n");
	else 
	  fprintf(stderr, "# norm=MAX\n");
	

	// estimate average sizes of all elements
	int* sizes = (int*)calloc(numgp, sizeof(int));
	int* nums  = (int*)calloc(numgp, sizeof(int));
	int* bins  = (int*)calloc(numgp, sizeof(int));
	
	for (i=0; i<numgenes; i++) {
		
	  // may skip if genelist specified and not in list
		if (genelist != 0) {
			
			/* query */
			e.key = a_genes[i].tsname;        
			my_hsearch_r(e, FIND, &ep, hash_genes) ;  
			if (!ep)
				continue;
		}
		
		int numexons   = a_genes[i].numexons;
		int st         = a_genes[i].fr;
		int** a_exons  = a_genes[i].exons;
		int tst        = a_genes[i].i; // transcript (not stranded)
		int ten        = a_genes[i].j;

		sizes[0] += promlen; // prom (always 0)
		nums [0] ++;
		
		sizes[numgp-1] += promlen; // downstream (numgp-1)
		nums [numgp-1] ++;
		
		// type 1, all gene structures considered
		if (strcmp(geneparts, "all") == 0) {
		  
		  // here only we need to loop over exons
		  for (j=0; j<numexons; j++) {
		    
		    // get pos exon
		    posexon = j;
		    if (st == -1)  // reverse if reverse strand of course
		      posexon = numexons - j - 1;
		    
		    // start intron if num exon > 0
		    if (j > 0) {
		      li =  a_exons[j][0] - a_exons[j-1][1] + 1;	  
		      if (posexon == 1) {
			sizes[2] += li;
			nums [2] ++ ;
		      } else if (posexon == 2) {
			sizes[4] += li;
			nums [4] ++;
		      } else {
			sizes[6] += li;
			nums [6] ++;
		      }	  
		    }
		    
		    // check whether exon long enough to matter
		    le = a_exons[j][1] - a_exons[j][0] + 1;	
		    if (posexon == 0) {
		      sizes[1] += le;
		      nums [1] ++;
		    } else if (posexon == 1) {
		      sizes[3] += le;
		      nums [3] ++;
		    } else {
		      sizes[5] += le;
		      nums [5] ++;
		    }
		  } // loop over exons
		}  // if all gene structures
		else if (strcmp(geneparts, "body") == 0) {
		  // only gene bodies
		  le = ten - tst + 1;
		  sizes[1] += le;
		  nums[1] ++;		    
		  
		} // end gene body
	
	} // loop over genes    
	
	// average out and plot
	for (i=0; i<numgp; i++) {
	  fprintf(stderr, "AVG=%f\n", sizes[i]/(double)(nums[i]));
	  sizes[i] = (int)( 0.5 + sizes[i] / (double)(nums[i]) );
	  
	  if (strcmp(geneparts, "all") == 0) { 
	    bins[i]  = (int)( sizes[i] / (double)q );    // average size of region divided size of quantized region (q) ==> number of bins for that region
	  } else 
	    bins[i] = numbins;			\
	}
	
	for (i=0; i<numgp; i++) {
		fprintf(stderr, "BINS: %d\n", bins[i]);
	}
	
	if (verbose == 1)
		fprintf(stderr, "Found %d genes to evaluate.\n", numgenes);
	
	//
	// end read intervals
	//
	
	// calc total number of reads for RPKM-style (not used right now)
	if (normalize == RPKM) {
	  totalnumreads = 0;
	  for (c=0; c<numchroms; c++) {
	    chrname = chrnames[c];
	    sprintf(readfile, "%s/reads.%s", chipdir, chrname); 
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
	
	double** intdata    = (double**)calloc(numgp, sizeof(double*));
	double** tmpintdata = (double**)calloc(numgp, sizeof(double*));
	int**    numdata    = (int**)   calloc(numgp, sizeof(int*));
	int**    tmpnumdata = (int**)   calloc(numgp, sizeof(int*));
	for (i=0; i<numgp; i++) {
	  // this will receive local (gene) and global (average) data coutnts
	  intdata[i]    = (double*)calloc(bins[i], sizeof(double));
	  tmpintdata[i] = (double*)calloc(bins[i], sizeof(double));
	  numdata[i]    = (int*)   calloc(bins[i], sizeof(int   ));
	  tmpnumdata[i] = (int*)   calloc(bins[i], sizeof(int   ));
	}
	
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
		
		
		if ((onlychr != 0) && (strcmp(onlychr, chrname) != 0))
			continue;
		
		if (verbose == 1)
			printf("loading data for %s\n", chrname);
		
		// alloc counte vectors
		chip_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
		if (chip_counts == 0) {	    
			die("Problem allocating chip_counts.\n");
		}

		/*
		if (inputdir != 0) {
			input_counts = (unsigned short*)calloc(chrlen, sizeof(unsigned short));
			if (input_counts == 0) {
				die("Problem allocating chip_counts.\n");
			}
		}
		*/
		
		// build a chr-wide read count profile
		sprintf(readfile, "%s/reads.%s", chipdir, chrname); 
		getCountFromReadFile(readfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &chip_counts, &numreads_chip, &numreadsclonal_chip);
		
		/*
		if (inputdir != 0) {
			sprintf(readfile, "%s/reads.%s", inputdir, chrname); 
			getCountFromReadFile(readfile, format, chrname, chrlen, readlen, fraglen, uniquereads, 1, &input_counts, &numreads_input, &numreadsclonal_input);
			
		}
		*/
		
		if (verbose == 1)
			fprintf(stderr, "Number of reads = %d (chip) and %d (input)\n", numreads_chip, numreads_input);
		
		
		// now go thru all genes
		for (i=0; i<numgenes; i++) {
			
			if ( strcmp(a_genes[i].chr, chrname) != 0 )
				continue;
			
			if ((onlyts != 0) && (strcmp(onlyts, a_genes[i].tsname) != 0))
				continue;
			
			if (genelist != 0) {
				
				/* query */
				e.key = a_genes[i].tsname;        
				my_hsearch_r(e, FIND, &ep, hash_genes) ;  
				if (!ep)
					continue;
			}
			
			if (verbose == 1)
				printf("%s\n", a_genes[i].tsname);
			
			// reset tmp data
			for (j=0; j<numgp; j++) {
				for (nq=0; nq<bins[j]; nq++) {
					tmpintdata[j][nq] = 0.0;
					tmpnumdata[j][nq] = 0;
				}
			}
			
			// transcript
			int t_st       = a_genes[i].i;
			int t_en       = a_genes[i].j;
			int numexons   = a_genes[i].numexons;
			int** a_exons  = a_genes[i].exons;
			int st         = a_genes[i].fr;
			
			// define where promoter starts and end
			int p_st = -1;
			int p_en = -1;
			if (st == 1) {
				p_st = t_st - promlen;
				p_en = t_st;
			} else {
				p_st = t_en;
				p_en = t_en + promlen;
			}	
						
			//
			// quantize promoter (same for all geneparts, bins[0] ha been set to numbins )
			//
			QuantizeReadCountsOverInterval(p_st, p_en, bins[0], chip_counts, &qavgcounts);
			for (nq=0; nq<bins[0]; nq++) {
				// reverse if strand asks for it
				efnq = nq;
				if (st == -1) 
					efnq = bins[0] - nq - 1;
				tmpintdata[0][efnq] = qavgcounts[nq];
				tmpnumdata[0][efnq] ++;	
			}	
			free(qavgcounts);
	
			// here depends on gp
			// start with "all"
			if (strcmp(geneparts, "all") == 0) {
			  
			  //
			  // quantize exons
			  //
			  if (verbose == 1)
			    printf("numexons = %d\n", numexons);
			  
			  for (j=0; j<numexons; j++) {
			    
			    // get pos exon
			    posexon = j;
			    if (st == -1)
			      posexon = numexons - j - 1;
			    
			    // start intron if num exon > 0
			    if (j > 0) {
			      
			      if (verbose == 1)
				printf("Intron %d\n", j-1);
			      
			      // determine numbins for Quantize function
			      if (posexon == 1)
				numbins = bins[2];
			      else if (posexon == 2)
				numbins = bins[4];
			      else 
				numbins = bins[6];
			      
			      // check whether intron long enough to matter
			      le =  a_exons[j][0] - a_exons[j-1][1] + 1;
			      if (le / numbins > 0) {
				
				QuantizeReadCountsOverInterval(a_exons[j-1][1], a_exons[j][0], numbins, chip_counts, &qavgcounts);
				
				for (nq=0; nq<numbins; nq++) {	    
				  // reverse if strand asks for it
				  efnq = nq;
				  if (st == -1) 
				    efnq = numbins - nq - 1;
							
				  // if intron 1, add to avg intron 1
				  if (posexon == 1) {
				    tmpintdata[2][efnq] = qavgcounts[nq];
				    tmpnumdata[2][efnq] ++;
				    // sam if intron 2
				  } else if (posexon == 2) {
				    tmpintdata[4][efnq] = qavgcounts[nq];
				    tmpnumdata[4][efnq] ++;
				    // all other introns
				  } else {
				    tmpintdata[6][efnq] += qavgcounts[nq];
				    tmpnumdata[6][efnq] ++;
				  }
				  
				} // end for (nq
				// end intron
				
				free(qavgcounts);
			      } // if intron long enough
			      
			    } // if (j>0)
			    
			    if (verbose == 1)
			      printf("Exon %d: %d-%d \n", j, a_exons[j][0], a_exons[j][1]);
			    
			    if (posexon == 0)
			      numbins = bins[1];
			    else if (posexon == 1)
			      numbins = bins[3];
			    else 
			      numbins = bins[5];
			    
			    // check whether exon long enough to matter
			    le = a_exons[j][1] - a_exons[j][0] + 1;
			    if (le / numbins == 0) {
			      continue;
			    }
			    
			    QuantizeReadCountsOverInterval(a_exons[j][0], a_exons[j][1], numbins, chip_counts, &qavgcounts);
			    
			    for (nq=0; nq<numbins; nq++) {
			      
			      // reverse if strand asks for it
			      efnq = nq;
			      if (st == -1) 
				efnq = numbins - nq - 1;
			      
			      // if exon 1, add to avg exon 1
			      if (posexon == 0) {
				tmpintdata[1][efnq] = qavgcounts[nq];
				tmpnumdata[1][efnq] ++;
				// sam if exon 2
			      } else if (posexon == 1) {
				tmpintdata[3][efnq] = qavgcounts[nq];
				tmpnumdata[3][efnq] ++;
				// all other exons
			      } else {
				tmpintdata[5][efnq] += qavgcounts[nq];
				tmpnumdata[5][efnq] ++;
			      }
					
			    } // end for (nq
			    // end exon
			    
			    free(qavgcounts);
			    
			  } // master for loop over exons (j=idx)
			} // end case 1, genepaers == "all"
			else {
			  
			  // quantize gene bodies (1=gb idx)
			  numbins = bins[1];
			  QuantizeReadCountsOverInterval(t_st, t_en, numbins, chip_counts, &qavgcounts);
			  for (nq=0; nq<numbins; nq++) {
			    // reverse if strand asks for it
			    efnq = nq;
			    if (st == -1) 
			      efnq = numbins - nq - 1;
			    tmpintdata[1][efnq] = qavgcounts[nq];
			    tmpnumdata[1][efnq] ++;	
			  }	
			  free(qavgcounts);  

			} // else gp == "body"
	  			  
			if (verbose == 1)
			  printf("Downstream\n");
			
			// find downstream boundaries
			int d_st = -1;
			int d_en = -1;
			if (a_genes[i].fr == 1) {
			  d_st = t_en;
			  d_en = t_en + promlen;
			} else {
			  d_st = t_st - promlen;
			  d_en = t_st;
			}
			// quantize downstream
			numbins = bins[numgp-1];
			QuantizeReadCountsOverInterval(d_st, d_en, numbins, chip_counts, &qavgcounts);
			for (nq=0; nq<numbins; nq++) {
			  // reverse if strand asks for it
			  efnq = nq;
			  if (st == -1) 
			    efnq = numbins - nq - 1;
			  tmpintdata[numgp-1][efnq] = qavgcounts[nq];
			  tmpnumdata[numgp-1][efnq] ++;	
			}	
			free(qavgcounts);
			
			
			if (strcmp(geneparts, "all") == 0) {
			  // average all introns, and all exons
			  for (nq=0; nq<bins[5]; nq++) {
			    if (tmpnumdata[5][nq] > 0) 
			      tmpintdata[5][nq] /= tmpnumdata[5][nq];
			  }
			  for (nq=0; nq<bins[6]; nq++) {
			    if (tmpnumdata[6][nq] > 0) 
			      tmpintdata[6][nq] /= tmpnumdata[6][nq];
			  }
			}
			
			
			// find max density (read count) across the entire gene
			double maxdensity = 0;
			if (verbose == 1)
				printf("NUM");
			// loop over gene parts
		
			for (j=0; j<numgp; j++) {
			  // loop over bins
			  for (nq=0; nq<bins[j]; nq++) {
			    if (verbose == 1)
			      printf("\t%d", tmpnumdata[j][nq]);
			    if (tmpintdata[j][nq] > maxdensity) {
			      maxdensity = tmpintdata[j][nq];
			    }	  
			  }
			  if (verbose == 1)
			    printf("\n\n");	
			}
			
			if (verbose == 1)
				printf("\nmaxd=%f\n", maxdensity);						
			
			// normalize and add
			for (j=0; j<numgp; j++) {
			  for (nq=0; nq<bins[j]; nq++) {
			    if (tmpnumdata[j][nq] > 0) {
			      
			      if (normalize == MAX) { // local max based normalization 
				if ( fabs(tmpintdata[j][nq]) < 1e-3)	
				  tmpintdata[j][nq] = 0.0;
				else 
				  tmpintdata[j][nq] /= maxdensity;			      			      
				if (tmpintdata[j][nq] > 5) {
				  die("wht ?\n");
				}
			      } else if (normalize == RPKM) { // RPKM normalization
				tmpintdata[j][nq] *= normfactor;
			      }
			      
			      if (verbose == 1)
				printf("\t%3.2f", tmpintdata[j][nq]);
			      intdata[j][nq] += tmpintdata[j][nq];
			      numdata[j][nq] ++;
			    } else {
			      if (verbose == 1)
				printf("\t%3.2f", -1.0);
			    }
			  }
			  if (verbose == 1)
			    printf("\n\n");
			  
			} // loop over gp
		
			// if outdata defined, write gene profile
			if (outdata != 0) {
			  for (j=0; j<numgp; j++) {
			    // select which label to show
			    if (strcmp(geneparts, "all") == 0)
			      lab = labelsALL[j];
			    else 
			      lab = labelsBO[j];
			    fprintf(fpod, "%s\t%s", a_genes[i].tsname, lab);

			    for (nq=0; nq<bins[j]; nq++) {
			      fprintf(fpod, "\t%3.2f",  tmpintdata[j][nq]);
			    }
			    fprintf(fpod, "\n");
			  } // loop over gp
			} // if outdata
			      if (verbose == 1)
			  printf("\n");
			
		} // loop over genes
		
		free(chip_counts);
		
		// loop over chr
	}
	
	// normalize and 
	
	if (verbose == 1)
		printf("FINAL\n");
       
	for (j=0; j<numgp; j++) {

	  // select which label to show
	  if (strcmp(geneparts, "all") == 0)
	    lab = labelsALL[j];
	  else 
	    lab = labelsBO[j];

	  // if internal, show average
	  if ((j > 0) && (j < numgp-1)) {	   
	    printf("%s (avg=%d bp)", lab, sizes[j]);
	  } else { // otherwise show fixed size
	    printf("%s (%d bp)", lab, sizes[j]);
	  }
	  // output actual profile
	  for (nq=0; nq<bins[j]; nq++) {
	    printf("\t%3.2f(%d)", intdata[j][nq] / numdata[j][nq], numdata[j][nq]);
	  }
	  printf("\n");
	}
	
	if (outdata != 0)
	  fclose(fpod);

	return 0;
	
}

void readGeneListIntoHash(char* file, struct my_hsearch_data** hash_genes)
{
	
	char*  buff;
	int    mynmax = 100000;
	
	char** a;
	int    m;
	FILE*  f;
	//int    cidx = -1;
	//int    i;
	//char*  line = 0;
	int    hashret;
	ENTRY  e;
	ENTRY* ep;
	
	//char** a_e_st;
	//char** a_e_en;
	//char** a_e_fr;
	//int    m_e_st;
	//int    m_e_en;
	//int    m_e_fr;
	//int**  a_myexons;
	//int    i, ii;
	
	// initialize num rnas per chr
	int numgenes = nbLinesInFile(file);
	
	/* build a hash table */
	*hash_genes = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
	hashret = my_hcreate_r( numgenes *2, *hash_genes);
	if (hashret == 0) {
		printf("Could not create hash table ...\n");
		exit(0);
	}
	
	// read set of intervals
	buff  = (char*)malloc(mynmax * sizeof(char));  
	f = fopen(file, "r");
	if (f == 0) 
		die("Cannot open refGene file.\n");
	
	while (fgets(buff, mynmax, f) != 0) {
		chomp(buff);
		
		
		split_line_delim(buff, "\t", &a, &m);
		
		char* n      = a[0];
		
		/* query */
		e.key = n;        
		my_hsearch_r(e, FIND, &ep, *hash_genes) ;  
		if (ep) {
			/* success, already rgere */
			free(a);
			continue;
		} else {
			
			// add
			/* enter key/value pair into hash */
			e.key   = strdup( n );
			e.data  = (char*)(1);
			hashret = my_hsearch_r(e, ENTER, &ep, *hash_genes);
			if (hashret == 0) {
				printf("Could not enter entry into hash table ...\n");
				exit(0);
			}
		}
		free(a);
	}
}


//
// read exon data and promoters
//
void readGenesAndExonDataFromRefGene(char* file, GenInt** a_genes, int* numgenes)
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
	char** a_e_st;
	char** a_e_en;
	char** a_e_fr;
	int    m_e_st;
	int    m_e_en;
	int    m_e_fr;
	int**  a_myexons;
	int    i, ii;
	//int*   sizes;
	
	// initialize num rnas per chr
	*numgenes = nbLinesInFile(file);
	
	/* build a hash table */
	hash_genes = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
	hashret = my_hcreate_r((*numgenes)*2, hash_genes);
	if (hashret == 0) {
		printf("Could not create hash table ...\n");
		exit(0);
	}
	
	// alloc intervals
	*a_genes = (GenInt*)malloc( (*numgenes) * sizeof(GenInt));
	if (*a_genes == 0) {
		die("Problem allocating a_rnas.\n");
	}
	
	// read set of intervals
	buff  = (char*)malloc(mynmax * sizeof(char));  
	f = fopen(file, "r");
	if (f == 0) 
		die("Cannot open refGene file.\n");
	*numgenes = 0;
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
			free(line);
			free(a);
			continue;
		} else {
			
		        // add
			/* enter key/value pair into hash */
			e.key   = strdup( n );
			e.data  = (char*)(*numgenes);
			hashret = my_hsearch_r(e, ENTER, &ep, hash_genes);
			if (hashret == 0) {
				printf("Could not enter entry into hash table ...\n");
				exit(0);
			}
			
		}
		
		char* c      = a[2];     //chr
		char  fr     = a[3][0];
		//int   tss    = (fr=='+'?atoi(a[4]):atoi(a[5]));
		int   t_st   = atoi(a[4]);
		int   t_en   = atoi(a[5]);
		char* g      = a[12];
		char* e_st   = a[9];
		char* e_en   = a[10];
		char* e_fr   = a[15];
		//int   cds_st = atoi(a[6]);
		//int   cds_en = atoi(a[7]);
		
		(*a_genes)[ (*numgenes) ].str       = line;
		(*a_genes)[ (*numgenes) ].tsname    = strdup(n);
		(*a_genes)[ (*numgenes) ].chr       = strdup(c);
		(*a_genes)[ (*numgenes) ].genename  = strdup(g);
		(*a_genes)[ (*numgenes) ].i         = t_st;
		(*a_genes)[ (*numgenes) ].j         = t_en;
		
		// frame
		(*a_genes)[ (*numgenes)].fr        = (fr=='+'?1:-1);
		
		
		//
		// EXONS
		//
		
		// chomp
		int  l_e_st = strlen(e_st);
		if (e_st[l_e_st-1] == ',')
			e_st[l_e_st-1] = '\0';    
		int  l_e_en = strlen(e_en);
		if (e_en[l_e_en-1] == ',')
			e_en[l_e_en-1] = '\0';    
		int  l_e_fr = strlen(e_fr);
		if (e_fr[l_e_fr-1] == ',')
			e_fr[l_e_fr-1] = '\0';
		
		// split
		split_line_delim(e_st, ",", &a_e_st, &m_e_st);
		split_line_delim(e_en, ",", &a_e_en, &m_e_en);
		split_line_delim(e_fr, ",", &a_e_fr, &m_e_fr);
		
		// exon store
		a_myexons = (int**)malloc(m_e_st * sizeof(int*));
		for (i=0; i<m_e_st; i++) {
			
			// take into account strand
			//if (fr == '+')
			ii = i;
			//else 
			//ii = m_e_st - i - 1;
			
			a_myexons[i]    = (int*)calloc(4, sizeof(int));
			a_myexons[i][0] = atoi(a_e_st[ii]);
			a_myexons[i][1] = atoi(a_e_en[ii]);
			a_myexons[i][2] = a_myexons[i][0];
			a_myexons[i][3] = a_myexons[i][1];
			
		}
		
		(*a_genes)[ (*numgenes) ].exons     = a_myexons;
		(*a_genes)[ (*numgenes) ].numexons  = m_e_st;
		//
		// END EXONS
		//
		free(a_e_st);
		free(a_e_en);
		free(a_e_fr);
		(*numgenes) ++;
		
		//if (*numgenes == 10)
		//  break;
        
		free(a);
		
	}
	
	free(buff);
	fclose(f);
	
}

//
// this function calculates 
//
void QuantizeReadCountsOverInterval(int i, int j, int q, unsigned short* counts, double** qavgcounts)
{
	
	int le = j - i + 1;
	int lq = le / q;  // floored on purpose
	int iq = 0;
	int jq = 0;
	int nq = 0;
	
	double avgcount = 0;
	int k = 0;
	
	if (j <= i) {
	  die("j <= i ! shouldn't happen\n");
	}

	if (lq == 0) {
		die("lq = 0, why ?\n");
	}
	
	*qavgcounts = (double*)calloc(q, sizeof(double));  
	
	for (nq=0; nq<q; nq++) {
		iq = i + nq * lq;
		jq = iq + lq;
		
		if (j < jq)
			jq = j;
		
		// get average read count in that interval
		avgcount = 0;
		for (k=iq; k<=jq; k++) {
			avgcount += counts[k]; 
		}
		avgcount /= lq;
		
		//if (verbose == 1)
		//  printf("Avg over [%d;%d] = %f\n", iq, jq, avgcount);
		
		(*qavgcounts)[nq] = avgcount;
	}
	
}

