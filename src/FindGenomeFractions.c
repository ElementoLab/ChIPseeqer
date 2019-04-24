/*
 * This program:
 * (A)	(Default) estimates the fraction of the genome in: (1) promoters,
 *		(2) downstream extremities (around TES), (3) exons, (4) introns,
 *		(5) distal and (6) intergenic.
 * (B)	If show_distal=1, prints the distal regions of the genome.
 *		If show_intergenic=1, prints the intergenic regions of the genome.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <search.h>
#include "dataio.h"
#include "statistics.h"
#include "sequences.h"
#include "hashtable.h"

#define FREE(p) do { free(p); p = NULL; } while(0)

typedef struct _GenInt {
	int   i;
	int   j;
	int   ci;
	int   cj;
	char* c;
	char* chr;
	
	char* id;
	char* dbsnp;
	int   numexons;
	int*  exon_starts;
	int*  exon_ends;
	int   strand;
	char* tsname;
	
	char* genename;
	int   fr;
	int   len;
	int   cdsst; 
	int   cdsen;
	int   ef; // exon frame
	int   num; // exon number
	char  utrtype;
	int   mut[4];
	int   N;
	int   k;
	float score;
	char  nc;
	float entropy;
	int   tes;
	int   tss;
	
} GenInt;


void readGenesFromAnnotation(char* file, GenInt** a_genes, int* numgenes, HASH* hash_genes, int is_refgene);

int main(int argc, char** argv) {
	
	char*  chrdata		= 0;
	int    numchroms	= 0;
	int*   chrlens		= 0;
	char** chrnames		= 0;
	float* chrmap		= 0;
	char*  chrname		= 0;
	int    c			= 0;
	int    chrlen		= 0;
	int    up        = 2000;
	int    dn        = 2000;
	int    maxdist   = 50000;
	
	unsigned char*  a_exons		= 0;
	unsigned char*  a_introns	= 0;
	unsigned char*  a_prom		= 0;
	unsigned char*  a_down		= 0;
	unsigned char*  a_extgenes  = 0;
	
	int	i	= 0;
	int	j	= 0;
	int k	= 0;
	char*  annotation	= 0;
	GenInt* a_genes;
	int     numgenes	= 0;
	HASH hc;
	
	int extgene_start	= 0;
	int extgene_end		= 0;
	
	double totalnt			= 0.0;
	long int exonnt			= 0;
	long int intronnt		= 0;
	long int downnt			= 0;
	long int promnt			= 0;
	long int distalnt		= 0;
	long int intergenicnt	= 0;
	
	int show_distal			= 0;
	int show_intergenic		= 0;
	
	int is_refgene			= 1;
	
	if (exist_parameter(argc, argv, "-chrdata")) {
		chrdata = get_parameter(argc, argv, "-chrdata");
		readChrData(chrdata, &chrnames, &chrlens, &chrmap, &numchroms);
	} else {
		die("Please provide -chrdata\n");
	}
	
	if (exist_parameter(argc, argv, "-annotation"))
		annotation = get_parameter(argc, argv, "-annotation");
	
	if (exist_parameter(argc, argv, (char*)"-show_distal"))
		show_distal = atoi(get_parameter(argc, argv, (char*)"-show_distal"));
	
	if (exist_parameter(argc, argv, (char*)"-show_intergenic"))
		show_intergenic = atoi(get_parameter(argc, argv, (char*)"-show_intergenic"));
	
	if (exist_parameter(argc, argv, (char*)"-is_refgene"))						/* It should be set to 0 if annotation file is not RefSeq. Otherwise, returns segmentation fault */
		is_refgene = atoi(get_parameter(argc, argv, (char*)"-is_refgene"));
	
	readGenesFromAnnotation(annotation, &a_genes, &numgenes, &hc, is_refgene);
	
	for (c=0; c<numchroms; c++) {
		//	for (c=4; c<7; c++) {
		
		chrlen		= chrlens[c];
		chrname		= chrnames[c];
		
		totalnt		+= (double)chrlen;
		
		a_exons		= create_binarized_array(chrlen);
		a_introns	= create_binarized_array(chrlen);
		a_prom		= create_binarized_array(chrlen);
		a_down		= create_binarized_array(chrlen);
		a_extgenes  = create_binarized_array(chrlen);
		
		for (i=0; i<numgenes; i++) {
			
			if (strcmp(chrname, a_genes[i].chr) != 0)
				continue;
			
			int numexons = a_genes[i].numexons;
			
			// paint exons      
			for (j=0; j<numexons; j++) {	
				for (k=a_genes[i].exon_starts[j]; k<=a_genes[i].exon_ends[j]; k++) {
					set_entry_in_binarized_array(a_exons, k);
				}	
			}
			
			// paint introns      
			for (j=1; j<numexons; j++) {	
				for (k=a_genes[i].exon_ends[j-1]; k<=a_genes[i].exon_starts[j]; k++) {
					set_entry_in_binarized_array(a_introns, k);
				}
			}
			
			/* Make sure that the promoter and dw extremity does not fall out of the gene (otherwise error in sac.cer) */
			int prom_up		= -1;
			int prom_down	= -1;
			
			if (a_genes[i].tss+dn > chrlen) {
				prom_down = chrlen;
			}
			else {
				prom_down = a_genes[i].tss+dn;
			}

			if (a_genes[i].tss-up < 0) {
				prom_up = 0; 
			}
			else {
				prom_up = a_genes[i].tss-up;
			}

			// paint prom
			for (k=prom_up; k<prom_down; k++) {
				set_entry_in_binarized_array(a_prom, k);
			}
			
			int dw_up		= -1;
			int dw_down		= -1;
			
			if (a_genes[i].tes+dn > chrlen) {
				dw_down = chrlen;
			}
			else {
				dw_down = a_genes[i].tes+dn;
			}
			
			if (a_genes[i].tes-up < 0) {
				dw_up = 0; 
			}
			else {
				dw_up = a_genes[i].tes-up;
			}			
			
			// paint down
			for (k=dw_up; k<a_genes[i].tes+dn; k++) {
				set_entry_in_binarized_array(a_down, k);
			}
			
			extgene_start	= a_genes[i].tss-maxdist;
			extgene_end		= a_genes[i].tes+maxdist;
			
			// paint intergenic
			if (extgene_end > chrlen) {
				extgene_end = chrlen;
			}
			if (extgene_start < 0) {
				extgene_start = 0; 
			}
			
			/*for (k=extgene_start; k<extgene_end; k++) {
			 set_entry_in_binarized_array(a_extgenes, k);
			 }*/
			
			for (k=extgene_start; k<a_genes[i].tss; k++) {
				set_entry_in_binarized_array(a_extgenes, k);
			}
			
			for (k=a_genes[i].tss; k<extgene_end; k++) {
				set_entry_in_binarized_array(a_extgenes, k);
			}
			
		} // for 
		
		int dist_start			= 0;	/* start position of distal region */
		int dist_end			= 0;	/* end position of distal region */
		int dist_has_previous	= 0;	/* set to 1 if we are already in distal region */
		int dist_last_pos		= -1;	/* last position of distal region in the array */
		int print_dist_end		= -1;	/* set to 1 if we have printed the end position */
		
		int inter_start			= 0;
		int inter_end			= 0;
		int inter_has_previous	= 0;
		int print_inter_end		= -1;	/* set to 1 if we have printed the end position */
		
		/* Store the last position of distal region in the array */
		int l;
		for (l=chrlen-1; l>=0; l--) {
			if(get_entry_in_binarized_array(a_extgenes, l) == 1) {
				dist_last_pos = l;
				break;
			}
		}
		
		for (i=0; i<chrlen; i++) {
			if (get_entry_in_binarized_array(a_prom, i) == 1) {
				promnt ++;
			} else if  (get_entry_in_binarized_array(a_down, i) == 1) {
				downnt ++;
			} else if (get_entry_in_binarized_array(a_exons, i) == 1) {
				exonnt ++;
			} else if (get_entry_in_binarized_array(a_introns, i) == 1) {
				intronnt ++;
			} else if (get_entry_in_binarized_array(a_extgenes, i) == 1) {
				distalnt ++;
				
				if(show_distal == 1) {					/* The following code is when we want to extract distal regions */
					if(dist_has_previous == 0) {		/* That is the start, goes here only the first time */
						dist_start			= i;		/* set start = current position */
						dist_end			= i;		/* set end = current position */
						dist_has_previous	= 1;		/* set dist_has_previous = TRUE */
						print_dist_end		= 0;		/* set print_dist_end = FALSE, we haven't printed end position yet */
						
						fprintf(stdout, "%s\t%d\t", chrname, dist_start);	/* print chromosome and start position */
					}
					else {								/* If dist_has_previous = TRUE */
						if (i-dist_end == 1) {			/* If we are in the exact next position, we update */
							dist_end		= i;		/* set end = current position */
							
							if(dist_end - dist_last_pos == 0) {		/* If we are at the last distal position in the chromosome, we need to print end position */
								fprintf(stdout, "%d\n", dist_end);	/* that's it, we print end position */
								print_dist_end = 1;					/* set print_end = TRUE, we printed end position */
							}
							
						}
						else {									/* If we are far away from the previous position, we print and start a new region */
							fprintf(stdout, "%d\n", dist_end);	/* that's it, we print end position */
							
							dist_start			= i;			/* set start = current position */
							dist_end			= i;			/* set end = current position */
							dist_has_previous	= 1;			/* set dist_has_previous = TRUE */
							
							fprintf(stdout, "%s\t%d\t", chrname, dist_start);	/* print chromosome and start position */
							print_dist_end		= 0;							/* set print_end = FALSE, we haven't printed end position yet */
						}
					}
				}
				
			} else {
				intergenicnt ++;
				
				if(show_intergenic == 1) {				/* The following code is when we want to extract distal regions */
					if(inter_has_previous == 0) {		/* That is the start, goes here only the first time */
						inter_start			= i;		/* set start = current position */
						inter_end			= i;		/* set end = current position */
						inter_has_previous	= 1;		/* set inter_has_previous = TRUE */
						print_inter_end		= 0;		/* set print_inter_end = FALSE, we haven't printed end position yet */
						
						fprintf(stdout, "%s\t%d\t", chrname, inter_start);	/* print chromosome and start position */
					}
					else {								/* If inter_has_previous = TRUE */
						if (i-inter_end == 1) {			/* If we are in the exact next position, we update */
							inter_end		= i;		/* set end = current position */
							
							print_inter_end	= 0;		/* set print_inter_end = FALSE, we haven't printed end position yet */
						}
						else {									/* If we are far away from the previous position, we print and start a new region */
							fprintf(stdout, "%d\n", inter_end);	/* that's it, we print end position */
							
							inter_start		= i;				/* set start = current position */
							inter_end		= i;				/* set end = current position */
							inter_has_previous	= 1;			/* set inter_has_previous = TRUE */
							
							fprintf(stdout, "%s\t%d\t", chrname, inter_start);	/* print chromosome and start position */
							print_inter_end		= 0;		/* set print_inter_end = FALSE, we haven't printed end position yet */
						}
					}
				}
			}
		}
		
		/* Before we go to the next chromosome, we check if the last end was printed. If not we print. */
		/* This prevents losing the end position of distal, when we have other regions -e.g., promoters- after the last distal. */
		if(show_distal == 1) {					
			if(print_dist_end == 0) {
				fprintf(stdout, "%d\n", dist_end);
			}
		}
		
		if(show_intergenic == 1) {					
			if(print_inter_end == 0) {
				fprintf(stdout, "%d\n", inter_end);
			}
		}
		
		FREE(a_exons);
		FREE(a_introns);
		FREE(a_prom);
		FREE(a_down);
		FREE(a_extgenes);
		
		fprintf(stderr, "# Done processing %s\n", chrname);
		
		//printf("Tot\t%e\n",  totalnt);
		
		/*
		 printf("Promoters\t%d\n",  100.0*(double)promnt/(double)totalnt);
		 printf("Downstream\t%d\n", 100.0*(double)downnt/(double)totalnt);
		 printf("Exons\t%d\n",      100.0*(double)exonnt/(double)totalnt);
		 printf("Introns\t%d\n",   100.0*(double)intronnt/(double)totalnt);
		 printf("Distal\t%d\n",   100.0*(double)distalnt/(double)totalnt);
		 */
		/*
		 printf("Promoters\t%3.2f\n",  100.0*(double)promnt/(double)totalnt);
		 printf("Downstream\t%3.2f\n", 100.0*(double)downnt/(double)totalnt);
		 printf("Exons\t%3.2f\n",      100.0*(double)exonnt/(double)totalnt);
		 printf("Introns\t%3.2f\n",   100.0*(double)intronnt/(double)totalnt);
		 printf("Distal\t%3.2f\n",   100.0*(double)distalnt/(double)totalnt);
		 */
		//break;
	}
	
	/* print the fractions of the genome */
	if(show_intergenic == 0 && show_distal == 0) {
		printf("Promoters\t%3.2f\n", 100.0*promnt/(double)totalnt);
		printf("Downstream\t%3.2f\n", 100.0*downnt/(double)totalnt);
		printf("Exons\t%3.2f\n", 100.0*exonnt/(double)totalnt);
		printf("Introns\t%3.2f\n", 100.0*intronnt/(double)totalnt);
		printf("Distal\t%3.2f\n", 100.0*distalnt/(double)totalnt);
		printf("Intergenic\t%3.2f\n", 100.0*intergenicnt/(double)totalnt);
	}
	
	/* cleanup */
	free(a_genes);
	
	
	return 0;
}


void readGenesFromAnnotation(char* file, GenInt** a_genes, int* numgenes, HASH* hc, int is_refgene)
{
	char*  buff;
	int    mynmax = 100000;  
	char** a;
	int    m;
	FILE*  f;
	
	char** a_e_st;
	char** a_e_en;
	char** a_e_fr;
	int    m_e_st;
	int    m_e_en;
	int    m_e_fr;
	int    i, j, l;
	int*   a_exon_starts;
	int*   a_exon_ends; 
	
	int    relposinexon;
	int    numexons       = 0;
	int    poscdsstingene = 0;
	int    poscdseningene = 0;
	int    posingene = 0;
	int    gidx = 0;
	
	if (is_refgene != 0 && is_refgene != 1) {
		die("Set the -is_refgene variable to 1, if the annotation file is RefGene. Otherwise set it to 0.\n");
	}
	
	// estimate line numbers
	*numgenes = nbLinesInFile(file);
	
	*a_genes = (GenInt*)calloc(*numgenes, sizeof(GenInt));
	if (*a_genes == 0) {
		die("Cannot alloc data for a_genes\n");
	}
	
	HASH_init(hc, (*numgenes) + 10000);
	
	// read genes
	*numgenes = 0;
	f = fopen(file, "r"); 
	if (f == 0) 
		die("Cannot open annotation file.\n");
	buff = (char*)calloc(mynmax, sizeof(char));
	while (!feof(f)) {
		fgets(buff, mynmax, f);
		if (feof(f))
			break; 
		chomp(buff);
		
		split_line_delim(buff, "\t", &a, &m);
		
		char* n      = a[1];  // refseq id
		char* c      = a[2];    
		char  fr     = a[3][0];
		char* e_st   = a[9];
		char* e_en   = a[10];
		int   cds_st = atoi(a[6]);
		int   tss    = (fr=='+'?atoi(a[4]):atoi(a[5]));
		int   tes    = (fr=='+'?atoi(a[5]):atoi(a[4]));
		
		(*a_genes)[ *numgenes ].tss = tss;
		(*a_genes)[ *numgenes ].tes = tes;
		
		int   cds_en = atoi(a[7])- 1;		// debug: fix 1-based coordinate
		
		// lookup
		if (HASH_find(hc, n, &gidx)) {
			free(a);
			continue;
		} 
		
		// else
		HASH_enter(hc, n, *numgenes);
		
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
		
		// split
		split_line_delim(e_st, ",", &a_e_st, &m_e_st);
		split_line_delim(e_en, ",", &a_e_en, &m_e_en);
		
		a_exon_starts   = (int*)malloc(m_e_st * sizeof(int));
		a_exon_ends     = (int*)malloc(m_e_st * sizeof(int));
		
		(*a_genes)[ *numgenes ].len = 0;
		
		// calculate ts length
		for (i=0; i<m_e_st; i++) {
			a_exon_starts[i] = atoi(a_e_st[i]);
			a_exon_ends[i]   = atoi(a_e_en[i])-1;		// fix 1-based coordinate
			int l = a_exon_ends[i] - a_exon_starts[i] + 1;
			(*a_genes)[ *numgenes ].len += l;
		}
		
		(*a_genes)[ *numgenes ].tsname      = strdup(n);
		(*a_genes)[ *numgenes ].chr         = strdup(c);
		(*a_genes)[ *numgenes ].fr          = (fr=='+'?1:-1);
		(*a_genes)[ *numgenes ].exon_starts = a_exon_starts;
		(*a_genes)[ *numgenes ].exon_ends   = a_exon_ends;
		(*a_genes)[ *numgenes ].numexons    = m_e_st;
		

		/* The following happens only when te annotation file is RefGene. Otherwise we get a segmentatio fault (e.g., no a[15] in Ensemblee). */
		if(is_refgene == 1) {
			char* e_fr   = a[15];
			char* g      = a[12];
			char* unk    = a[13];  // unk for ncRNA
			
			int  l_e_fr = strlen(e_fr);
			
			if (e_fr[l_e_fr-1] == ',')
				e_fr[l_e_fr-1] = '\0';
			
			split_line_delim(e_fr, ",", &a_e_fr, &m_e_fr);
			
			(*a_genes)[ *numgenes ].genename    = strdup(g);
			
			if (strcmp(unk, "unk") == 0) { // ncRNA
				(*a_genes)[ *numgenes ].nc = 1;
			} else {
				(*a_genes)[ *numgenes ].nc = 0;
			}
		}
		
		// calculate pos cds in 
		numexons       = m_e_st;
		poscdsstingene = 0;
		poscdseningene = 0;
		
		if (fr == '+') {
			
			posingene = 0;      
			// traverse exons
			for (j=0; j<numexons; j++) {
				
				// if cds_st in exon
				if ((a_exon_starts[j] <= cds_st) && (cds_st <= a_exon_ends[j])) {
					
					// relative positon wrt begining of exon
					// relposinexon = cds_st - a_exon_starts[j];
					relposinexon = cds_st - a_exon_starts[j] + 1;		// length (
					
					// pos cds wrt begining of gene
					poscdsstingene = posingene + relposinexon;
					
				}
				
				// if cds_en in exon
				if ((a_exon_starts[j] <= cds_en) && (cds_en <= a_exon_ends[j])) {
					
					// relative positon wrt begining of exon
					// relposinexon = cds_en - a_exon_starts[j];
					relposinexon = cds_en - a_exon_starts[j] + 1;  	// length (
					
					// pos cds wrt begining of gene
					poscdseningene = posingene + relposinexon;
					
				}
				
				l = a_exon_ends[j] - a_exon_starts[j] + 1;
				posingene += l;
				
			}
			poscdsstingene--;	// debug
			poscdseningene--;	// debug
		} // if (fr == +	
		else {
			// fr == -
			
			posingene = 0;      
			// traverse exons
			for (j=numexons-1; j>=0; j--) {
				//printf("exon %d\n", j);
				// if cds_st in exon
				if ((a_exon_starts[j] <= cds_st) && (cds_st <= a_exon_ends[j])) {
					
					// relative positon wrt begining of exon
					relposinexon   = a_exon_ends[j] - cds_st +1 ;  //leading length (
					
					// pos cds en wrt begining of gene (cds_st is actually cds_en of course)
					poscdseningene = posingene + relposinexon;
					
				}
				
				// if cds_en in exon
				if ((a_exon_starts[j] <= cds_en) && (cds_en <= a_exon_ends[j])) {
					//printf("found cds_en\n");
					// relative positon wrt begining of exon
					relposinexon = a_exon_ends[j] - cds_en + 1;  // trailing length
					
					// pos cds st wrt begining of gene
					poscdsstingene = posingene + relposinexon;
					
				}
				
				l = a_exon_ends[j] - a_exon_starts[j] + 1;
				posingene += l;
				//printf("posingene becomes %d\n", posingene);
			}
			
			//poscdsstingene++;
			poscdsstingene--;	//debug: length -> position 
			poscdseningene--;
		} // else fr == -
		
		(*a_genes)[ *numgenes ].cdsst = poscdsstingene;
		(*a_genes)[ *numgenes ].cdsen = poscdseningene;
		
		
		(*numgenes) ++;
		
		free(a_e_st);
		free(a_e_en);
		
		if(is_refgene == 1) {
			free(a_e_fr);
		}
		
		free(a);
		
	}
	
	/* cleanup */
	HASH_destroy(hc);
	free(a_exon_starts);
	free(a_exon_ends);
	free(buff);
	fclose(f);
}



