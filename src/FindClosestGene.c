#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "dataio.h"
#include "statistics.h"
#include "sequences.h"

typedef struct _GenInt {
	int i;			//start position
	int j;			//end position
	char* c;		//chromosome
	char* id;		//id
	float score;	//score
	int strand;		//strand
} GenInt;


int main(int argc, char** argv) {
	
	// general
	long i	= 0;
	long j	= 0;
	long k	= 0;
	
	int w	= 0;

	// files to compare
	FILE* f1;
	FILE* f2;
	
	// reading file
	char* buff			= NULL;
	int   mynmax		= 100000;
	char** p			= NULL;
	int    m			= 0;
	char* intervals1	= 0;
	char* intervals2	= 0;	
	int   verbose		= 0;	
	int   maxnumint		= 1000000000;
	
	// store intervals
	GenInt* a_int1;
	GenInt* a_int2;
	int     numint1 = 0;
	int     numint2 = 0;	
	int     ext1	= 0;
	int     ext2	= 0;
	int     hasid1	= 0;
	int     hasid2	= 0;
	int     rand2	= 0;
	
	int		num_close_genes = 1; //can be either 1 or 2
	
	int     showscores			= 0;
	int     showprofile			= 0;
	int     show_ov_int			= 0;
	int		showmin				= 0;
	int		show_gene_only		= 0;
	int		show_upstream_gene	= 0;
	int     showunion			= 0;
	char*   featname			= 0;
	
	int*    a_ov_int = 0;
	
	int     showpos  = 0;
	int     iswig    = 0;
	int     minx;
	int     maxx;
	int     hasheader1	= 0;
	int     hasheader2	= 0;
	
	int		distance_min		= 0;	//mimimum distance between intervals
	int		distance_max		= 0;	//maximum distance between intervals
	long	distance_estim		= 0;	//holds the computed distance
	long	orig_distance_estim = 0;
	
	//variables for minimum 1
	int		min1_flag = 0;
	int		min1	  = 0;
	char*	min1_id	  = 0;
	char*	min1_c	  = 0;
	int		min1_i	  = 0;
	int		min1_j	  = 0;
	int		min1_str  = 0;
	int		min1_dir  = 0;
	
	//variables for minimum 2	
	int		min2_flag = 0;
	int		min2	  = 0;
	char*	min2_id	  = 0;
	char*	min2_c	  = 0;
	int		min2_i	  = 0;
	int		min2_j	  = 0;
	int		min2_str  = 0;
	int		min2_dir  = 0;
	
	//variables for upstream gene	
	int		ups_flag = 0;
	int		ups		 = 0;
	char*	ups_id	 = 0;
	char*	ups_c	 = 0;
	int		ups_i	 = 0;
	int		ups_j	 = 0;
	int		ups_str  = 0;
	//int		ups_dir  = 0;
	
	int		dir		  = 0;	//directionality
	char*	peak_pos  = 0;
	
	int showstrand = 0;
	int  addheader = 0;

	if (argc < 2) {
		die("Usage: FindClosestGene -intervals1 FILE -ext1 INT -hasid1 INT -intervals2 FILE -ext2 INT -hasid2 INT -showprofile INT -distance_min INT -distance_max INT -num_close_genes INT \n");
	}
	
	if (exist_parameter(argc, argv, "-intervals1"))
		intervals1 = get_parameter(argc, argv, "-intervals1");

	if (exist_parameter(argc, argv, "-addheader"))
		addheader = atoi(get_parameter(argc, argv, "-addheader"));
	
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
	
	if (exist_parameter(argc, argv, "-distance_max"))
		distance_max = atoi(get_parameter(argc, argv, "-distance_max"));
	
	if (exist_parameter(argc, argv, "-num_close_genes"))
		num_close_genes = atoi(get_parameter(argc, argv, "-num_close_genes"));
	
	if (exist_parameter(argc, argv, "-showmin"))
		showmin = atoi(get_parameter(argc, argv, "-showmin"));
	
	if (exist_parameter(argc, argv, "-show_gene_only"))
		show_gene_only = atoi(get_parameter(argc, argv, "-show_gene_only"));
	
	if (exist_parameter(argc, argv, "-show_upstream_gene"))
		show_upstream_gene = atoi(get_parameter(argc, argv, "-show_upstream_gene"));
	
	default_set_seed(time(NULL));
	
	// alloc counter vectors
	int estnumint1 = nbLinesInFile(intervals1) + 1000;
	
	a_int1 = (GenInt*)malloc(estnumint1 * sizeof(GenInt));
	if (a_int1 == 0) {
		die("Problem allocating a_int1.\n");
	}
	int estnumint2 = nbLinesInFile(intervals2) + 1000;
	a_int2 = (GenInt*)malloc(estnumint2 * sizeof(GenInt));
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
		//char* line = strdup(buff);
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
		if (showstrand == 1) {
			a_int1[numint1].strand = atoi(p[idx+3]);
		}
		if (showscores == 1) {
			a_int1[numint1].score = atoi(p[idx+4]);
		}
		free(p);
		numint1++;
		
		if (numint1 == maxnumint) {
			die("Max number of intervals reached.\n");
		}
		
	}
	
	//printf("\n");
	
	
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
		a_int2[numint2].strand = atoi(p[idx+3]);
		
		//printf("%s\t%s\t%d\t%d\n", a_int2[numint2].id, a_int2[numint2].c, a_int2[numint2].i, a_int2[numint2].j);
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
	
	a_ov_int = (int*)calloc(100000, sizeof(int));

	if (addheader == 1) {
	  printf("chr\tstart\tend\tnumgenes\tgenename\n");
	}
	
	int numinter = 0;
	// for every interval of the first file
	for (i=0; i<numint1; i++) {
		if (hasid1 == 1) {
			printf("%s\t", a_int1[i].id);
		}
		if (showprofile == 0) {
			printf("%s\t%d\t%d\t", a_int1[i].c, a_int1[i].i, a_int1[i].j);
		}
		if (showstrand == 1) {
			printf("%d\t", a_int1[i].strand);
		}
		
		if (showscores == 1) {
			printf("%f\t", a_int1[i].score);
		}
		
		int cnt		= 0;
		min1_flag	= 0;
		min2_flag	= 0;
		peak_pos	= 0;
		
		//initialize min1 and min2
		min1 = 0;
		for (k=0; k<numint2; k++) {
			if(a_int2[k].i-ext2 > min1) {
				min1 = a_int2[k].i-ext2;
			}
			if(a_int2[k].j+ext2 > min1) {
				min1 = a_int2[k].j+ext2;
			}
		}
		min2 = min1;
		
		// for every interval of the second file
		for (j=0; j<numint2; j++) {
					
			minx = a_int1[i].i;
			maxx = a_int1[i].j;
			
			//count the distance
			orig_distance_estim = sequencesDistance(a_int1[i].i-ext1, a_int1[i].j+ext1, a_int2[j].i-ext2, a_int2[j].j+ext2);
			distance_estim		= abs(orig_distance_estim);
			
			if(j==0 && distance_estim!=0) {
				min1 = distance_estim+2;
			}
			
			if (orig_distance_estim > 0) {
				dir = 1;
			}
			else {
				dir = -1;
			}
			
			if ((strcmp(a_int1[i].c, a_int2[j].c) == 0) && (distance_estim > distance_min)) {
				
				//if maximum distance is set
				if (distance_max != 0) {

					//compare estimated distance with max distance
					if (distance_estim < distance_max) {
						//compare estimated distance with current min1 distance
						if (distance_estim < min1) {
							
							/*if(show_upstream_gene == 1) {
							 if ((a_int2[j].strand == 1) && (a_int1[i].i+ext1 > a_int2[j].i-ext2)) {
							 continue;
							 }
							 if ((a_int2[j].strand == -1) && (a_int1[i].j+ext1 < a_int2[j].j+ext2)) {
							 continue;
							 }
							 }*/
							
							//if min1 exists for that interval, set min2=min1
							if (min1_flag) {
								if (dir == min1_dir) {
								}
								else {
									min2_flag	= 1;
									min2		= min1;
									min2_c		= min1_c;
									min2_i		= min1_i;
									min2_j		= min1_j;
									min2_str	= min1_str;
									min2_dir	= min1_dir;
									if (hasid2 == 1) {
										min2_id = min1_id;
									}
								}
							}
							//update min1
							min1_flag = 1;
							min1	  = distance_estim;
							min1_c	  = a_int2[j].c;
							min1_i	  = a_int2[j].i;
							min1_j	  = a_int2[j].j;
							min1_str  = a_int2[j].strand;
							min1_dir  = dir;
							if (hasid2 == 1) {
								min1_id = a_int2[j].id;
							}
						}
						//compare estimated distance with current min2 distance
						else if (distance_estim < min2) {
							//update min2
							if (dir == min1_dir) {
							}
							else {
								min2_flag = 1;
								min2	  = distance_estim;
								min2_c	  = a_int2[j].c;
								min2_i	  = a_int2[j].i;
								min2_j	  = a_int2[j].j;
								min2_str  = a_int2[j].strand;
								min2_dir  = dir;
								if (hasid2 == 1) {
									min2_id = a_int2[j].id;
								}									
							}								
						}
					}
				}
				//if no need to compare with maximum distance
				else {
					//compare estimated distance with current min1 distance
					if (distance_estim < min1) {
						
						/*if(show_upstream_gene == 1) {
						 if ((a_int2[j].strand == 1) && (a_int1[i].i+ext1 > a_int2[j].i-ext2)) {
						 continue;
						 }
						 if ((a_int2[j].strand == -1) && (a_int1[i].j+ext1 < a_int2[j].j+ext2)) {
						 continue;
						 }
						 }*/
						
						//if min1 exists for that interval, set min2=min1
						if (min1_flag) {
							if (dir == min1_dir) {
							}
							else {
								min2_flag	= 1;
								min2		= min1;
								min2_c		= min1_c;
								min2_i		= min1_i;
								min2_j		= min1_j;
								min2_str	= min1_str;
								min2_dir	= min1_dir;
								if (hasid2 == 1) {
									min2_id = min1_id;
								}
							}
						}
						//update min1
						min1_flag = 1;
						min1	  = distance_estim;
						min1_c	  = a_int2[j].c;
						min1_i	  = a_int2[j].i;
						min1_j	  = a_int2[j].j;
						min1_str  = a_int2[j].strand;
						min1_dir  = dir;
						if (hasid2 == 1) {
							min1_id = a_int2[j].id;
						}
					}
					//compare estimated distance with current min2 distance
					else if (distance_estim < min2) {
						//update min2
						if (dir == min1_dir) {
						}
						else {
							min2_flag = 1;
							min2	  = distance_estim;
							min2_c	  = a_int2[j].c;
							min2_i	  = a_int2[j].i;
							min2_j	  = a_int2[j].j;
							min2_str  = a_int2[j].strand;
							min2_dir  = dir;
							if (hasid2 == 1) {
								min2_id = a_int2[j].id;
							}
						}
					}
				}
				
				if(show_upstream_gene == 1) {
					if (((min1_str == 1) && (a_int1[i].j+ext1 < min1_i-ext2)) || ((min1_str == -1) && (a_int1[i].i-ext1 > min1_j+ext2))) {
						ups_flag = 1;
						ups		= min1;
						ups_c	= min1_c;
						ups_i	= min1_i;
						ups_j	= min1_j;
						ups_str = min1_str;
						if (hasid2 == 1) {
							ups_id = min1_id;
						}
					}
					else if (((min2_str == 1) && (a_int1[i].j+ext1 < min2_i-ext2)) || ((min2_str == -1) && (a_int1[i].i-ext1 > min2_j+ext2))) {
						ups_flag = 1;
						ups		= min2;
						ups_c	= min2_c;
						ups_i	= min2_i;
						ups_j	= min2_j;
						ups_str = min2_str;
						if (hasid2 == 1) {
							ups_id = min2_id;
						}
					}
					else {
						ups_flag = 0;
					}
				}
				
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
		
		//printf("%d", mycnt);
		
		/*if (show_ov_int == 1) {
		 for (j=0; j<cnt; j++) {
		 printf("\t%s-%d-%d", 
		 a_int2[ a_ov_int[j] ].c,
		 a_int2[ a_ov_int[j] ].i,
		 a_int2[ a_ov_int[j] ].j);
		 if (hasid2 == 1) {
		 printf(":D-%s", a_int2[ a_ov_int[j] ].id);
		 }
		 }
		 }*/
		
		if(show_upstream_gene == 1) {
			if (ups_flag) {
				printf("%d", 1);
				if(show_gene_only == 1) {
					printf("\t%s", ups_id);
				}
				else {					
					printf("\t%s-%d-%d:D-%s", ups_c, ups_i, ups_j, ups_id);
					if(showmin) {
						printf("-%d", ups);
					}
				}
			}
			else {
				printf("0\t-");
			}
		}
		else if (num_close_genes == 1) {
			if (min1_flag) {
				printf("%d", 1);
				
				if(show_gene_only == 1) {
					printf("\t%s", min1_id);
				}
				else {
					printf("\t%s-%d-%d:D-%s", min1_c, min1_i, min1_j, min1_id);
					if(showmin) {
						printf("-%d", min1);
					}
				}
				
				if (showunion == 1) {
					
					if (min1_i < minx) {
						minx = min1_i;
					}
					
					if (min1_j > maxx)
						maxx = min1_j;
					
					printf("\t%s-%d-%d", min1_c, minx, maxx);
				}
				
			}
			else {
				printf("0");
			}
		}
		else if (num_close_genes == 2) {
			if (min1_flag && min2_flag) {
				
				printf("%d", 2);
				
				if(show_gene_only == 1) {
					printf("\t%s\t%s", min1_id, min2_id);
				}
				else {
					printf("\t%s-%d-%d:D-%s", min1_c, min1_i, min1_j, min1_id);
					if(showmin) {
						printf("-%d", min1);
					}
					
					printf("\t%s-%d-%d:D-%s", min2_c, min2_i, min2_j, min2_id);
					if(showmin) {
						printf("-%d", min2);
					}
				}
			}
			else {
				printf("0");
				
			}
		}
		
		printf("\n");
		if (cnt > 0) 
			numinter ++;
		
	}
	
	/* cleanup */
	
	for(w=0; w<numint1; w++) {
		free(a_int1[w].c);
		free(a_int1[w].id);
	}

	for(w=0; w<numint2; w++) {
		free(a_int2[w].c);
		free(a_int2[w].id);
	}
	
	fclose(f1);
	fclose(f2);
	free(a_int1);
	free(a_int2);
	free(buff);
	free(a_ov_int);
	
	return 0;
}
