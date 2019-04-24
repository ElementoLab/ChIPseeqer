/*
 * This program implements the Self Organizing Maps
 * (Kohonen) algorithm.
 */

#define _GNU_SOURCE
#include <search.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <sys/types.h>
#include <unistd.h>
#include <values.h>

#define  NBGENES 40000
#define  NBCONDS 1000

char* get_parameter(int argc, char** argv, char* param);
int exist_parameter(int argc, char** argv, char* param);

float alpha (int iter, int length, float alp0) ;
float neigh (float d, float s);
float radius (int iter, int length, float rad0, float radf);
float distance(float* ind1, float* ind2, int n);
float pearson(float* data1, float* data2, int m); 
float distOnMap (int x1, int y1, int x2, int y2); 
void write2Dmap(char* psname, float** data, int nborfs, int nbcond, int xdim, int ydim, int* affectx, int* affecty, int shade); 
void write2Dmap_general(char* psname, float** data, int nborfs, int nbcond, int xdim, int ydim, int* affectx, int* affecty, int shade); 
float estimate_tightness(float** data, int nborfs, int nbcond, int xdim, int ydim, int* affectx, int* affecty, char** names, int distance_method);

int verbose = 0;

int main (int argc, char** argv) {
	
	char*	expfile	= NULL;
	char*	buff	= NULL;
	FILE*	fp1		= NULL;
	int		nbcond	= 0;
	int		nborfs	= 0;
	float** data	= NULL;
	char**	names	= NULL;
	int		xdim	= 0;
	int		ydim	= 0;
	int		iter	= 0;
	int		init_iter = 1;
	
	int		i		= 0;
	int		j		= 0;
	int		r		= 0;
	int		k		= 0;
	char*	s		= NULL;
	int cnt			= 0;
	int exp_map		= 1;
	int show_shade	= 1;
	
	int* used		= NULL;
	char* outfile	= NULL;
	float*** neuron	= NULL;	// tri-dimensional array
	
	// training
	int		t		= 0;
	int		xwin	= 0;
	int		ywin	= 0;
	float	mindist	= 0.0;
	float	dist	= 0.0;
	
	float	rad0	= 2.5;
	float	radf	= 0.0; 
	float	alpha0	= 0.1;
	float	alp		= 0.0;
	float	dwi		= 0.0;
	float	hwi		= 0.0;
	float	rad		= 0.0;
	int		l		= 0;
	FILE*	fp		= NULL;
	int		tlog	= 0;
	int		tnor	= 0;
	int		tpos	= 0;
	int		distances = 0;
	int*	affectx = NULL;
	int*	affecty = NULL;
	int*	min_affectx = NULL;
	int*	min_affecty = NULL;
	float	sum_x	= 0.0;
	float	sum_x2	= 0.0;
	float	avg		= 0.0;
	float	dev		= 0.0;
	char* outepsmap = NULL;
	
	int distance_method		= 1;	/* 1 for euclidean, 0 for pearson */
	
	float* quality_array	= NULL;
	int   mynmax			= 100000;
	float min_qual			= MAXFLOAT;
	float tmp				= MAXFLOAT;
	
	if (argc < 4) {
		printf("Usage : kohonen -exp FILE -outclusters FILE [ -xdim INT -ydim INT -iter INT  -outepsmap FILE -tlog INT -tnor INT ]\n");
		printf("Explanations:\n");
		printf("All parameters between [ ] are optional.\n");
		printf("Mandatory parameters:\n");
		printf("-exp indicates the tab-delimited gene x conditions table (25000 genes and 1000 conditions maximum by default).\n");
		printf("see data/gasch_stress_microarray_data.txt for an example (gene names and condition names are required)\n");
		printf("-outclusters points the clustering partition output file (to be created).\n");
		printf("\n");
		printf("Optional parameters:\n");
		printf("-xdim and -ydim indicate the dimension of the map (number of nodes is -xdim times -ydim). Default is 10 x 10.\n");
		printf("-iter is the number of iterations. Default is 20.\n"); 
		printf("-outepsmap is the Postscript format map file. Note created by default.\n");
		printf("-tlog specifies whether to log transform the values. Choices are 0 (no) and 1 (yes). Default is 0.\n"); 
		printf("-tnor specifies whether to variance-normalize the rows. Choices are 0 (no) and 1 (yes). Default is 0.\n"); 
		printf("-expmap specifies whether to create a 2D map for expression data (will also print min/max/avg plots). Choices are 0 (no) and 1 (yes). Default is 1.\n"); 
		printf("-shade shows coloring shades for the nodes, according to the size of the nodes. Choices are 0 (no) and 1 (yes). Default is 1.\n"); 
		
		exit(0);
	}
	
	if (exist_parameter(argc, argv, "-exp")) {
		expfile =      get_parameter(argc, argv, "-exp");
	} else {
		printf("Please specify -exp\n");
		exit(0);
	}
	
	if (exist_parameter(argc, argv, "-outclusters")) {
		outfile =      get_parameter(argc, argv, "-outclusters");
	} else {
		printf("Please specify -outclusters\n");
		exit(0);
	}
	
	if (exist_parameter(argc, argv, "-outepsmap")) {
		outepsmap =      get_parameter(argc, argv, "-outepsmap");
	} else {
		outepsmap = 0;
	}
	
	if (exist_parameter(argc, argv, "-xdim")) {
		xdim =      atoi(get_parameter(argc, argv, "-xdim"));
	} else {
		xdim = 10;
	}
	
	if (exist_parameter(argc, argv, "-ydim")) {
		ydim =      atoi(get_parameter(argc, argv, "-ydim"));
	} else {
		ydim = 10;
	}
	
	if (exist_parameter(argc, argv, "-iter")) {
		iter =      atoi(get_parameter(argc, argv, "-iter"));
	} else {
		iter = 20;
	}
	
	if (strlen(get_parameter(argc, argv, "-verbose")) > 0) {
		verbose    = atoi(get_parameter(argc, argv, "-verbose"));
	} else {
		verbose    = 0;
	}
	
	tlog = 0;
	if (strcmp(get_parameter(argc, argv, "-tlog"), "1") == 0)
		tlog = 1;
	
	tnor = 0;
	if (strcmp(get_parameter(argc, argv, "-tnor"), "1") == 0)
		tnor = 1;
	
	tpos = 0;
	if (strcmp(get_parameter(argc, argv, "-tpos"), "1") == 0)
		tpos = 1;
	
	if (strcmp(get_parameter(argc, argv, "-distances"), "1") == 0)
		distances = 1;
	
	if (exist_parameter(argc, argv, "-expmap")) {
		exp_map = atoi(get_parameter(argc, argv, "-expmap"));
	}
	
	if (exist_parameter(argc, argv, "-shade")) {
		show_shade = atoi(get_parameter(argc, argv, "-shade"));
	}
	
	if (exist_parameter(argc, argv, "-init_iter")) {
		init_iter =      atoi(get_parameter(argc, argv, "-init_iter"));
	} else {
		init_iter = 1;
	}
	
	if (exist_parameter(argc, argv, "-distance_method")) {
		distance_method = atoi(get_parameter(argc, argv, "-distance_method"));
	}
	
	// input data
	printf("Making a %d x %d map.\n", xdim, ydim);
	printf("%d iterations\n", iter);
	
	fp1 = fopen(expfile, "r");
	if (!fp1) {
		printf("could not open expression data %s\n", expfile);
	}
	
	// array to store the values for each gene
	data = (float**)malloc(NBGENES * sizeof(float*));
	for (i=0; i<NBGENES; i++) {
		data[i] = (float*)malloc(NBCONDS * sizeof(float));
	}
	
	// array to store the names of the genes
	names = (char**)malloc(NBGENES * sizeof(char*));
	
	// create a line buffer
	buff = (char*)malloc(mynmax * sizeof(char));
	
	quality_array = (float*)malloc(init_iter * sizeof(float));
	
	fgets(buff, mynmax, fp1);
	
	i = 0;
	nborfs = 0;
	
	//read the input file
	while (!feof(fp1)) {
		fgets(buff, mynmax, fp1);
		
		if (feof(fp1)) {
			break;
		}
		
		if (verbose)
			printf("Read line %s .. \n", buff);
		
		j = strlen(buff);
		while ((j >= 0) && (buff[j] != '\n')) j--;
		
		buff[j] = '\0';
		
		s = strtok(buff, "\t");
		
		// store the gene name
		names[i] = strdup(s);
		
		j = 0;
		nbcond = 0;
		while ((s = strtok(0, "\t"))) {
			
			// store the data values for each gene
			data[i][j] = atof(s);
			
			if (tlog && (data[i][j] > 0.0))
				data[i][j] = log(data[i][j]);
			
			j++;
			
			nbcond ++;
		}
		nborfs++;
		i++;
		
		if (i == NBGENES) {
			printf("NBGENES reached .. dying.\n");
			exit(0);
		}
		
	}
	
	fclose(fp1);
	
	// reduce - center ?
	if (tnor) {
		
		for (k=0; k<nborfs; k++) {
			
			sum_x  = 0.0;
			sum_x2 = 0.0;
			
			// add all values
			for (l=0; l<nbcond; l++) {
				sum_x  += data[k][l];
				sum_x2 += data[k][l] * data[k][l];
			}
			
			// divide by the number of conditions
			avg  = sum_x / nbcond;
			dev  = ( sum_x2 - ( sum_x * sum_x / nbcond ) ) / ( nbcond - 1 );
			
			// reduce - center
			for (l=0; l<nbcond; l++) {
				
				data[k][l] = data[k][l] - avg;
				
				if (fabs(data[k][l]) > 0.000001) {
					data[k][l] = data[k][l] / sqrt(dev);
				}
			}			
		}
	}
	
	if (verbose == 1) {
		printf("nbcond = %d\n", nbcond);
	}
	
	/* how many times we want the algorithm to run, to average the tightness of the clusters */
	int q = 0;
	
	for (q=0; q<init_iter; q++) {
		
		neuron = (float***)malloc(xdim * sizeof(float**));
		for (i=0; i<xdim; i++) {
			neuron[i] = (float**)malloc(ydim * sizeof(float*));
			
			for (j=0; j<ydim; j++) {
				neuron[i][j] = (float*)malloc(nbcond * sizeof(float));			
			}
		}
		
		/* init the pseudo-random generator */
		srandom((unsigned)(getpid() + time(NULL)));
		
		// allocate an array for the genes that have been used
		used = (int*)malloc(NBGENES * sizeof(int));
		for (i=0; i<nborfs; i++) {
			used[i] = 0;
		}  
		
		for (i=0; i<xdim; i++) {
			
			for (j=0; j<ydim; j++) {
				
				r = 1+ (int) ((double)nborfs * random() / (RAND_MAX+1.0));
				
				while (used[r] == 1) {
					r = 1+ (int) ((double)nborfs * random() / (RAND_MAX+1.0));				
				}
				
				used[r]     = 1;
				
				for (k=0; k<nbcond; k++) {
					neuron[i][j][k] = data[r][k];
				}
			}
		} 
		
		// train Kohonen
		for (t=0; t<iter; t++) {
			
			printf("Iteration %d ..\n", t);
			
			cnt = 0;
			
			alp = alpha (t, iter, alpha0);
			rad = radius(t, iter, rad0, radf);
			
			//printf("alp=%5.4f, rad=%5.4f\n", alp, rad);
			
			for (k=0; k<nborfs; k++) {
				
				// printf("orf %s ..\n", names[k]);
				mindist = 100000.0;
				
				// determine the distance for each gene 
				for (i=0; i<xdim; i++) {
					for (j=0; j<ydim; j++) {
						
						dist = distance(data[k], neuron[i][j], nbcond);
						cnt++;
						
						if (dist < mindist) {
							xwin    = i;
							ywin    = j;
							mindist = dist;
						}
					}
				}
				
				for (i=0; i<xdim; i++) {
					for (j=0; j<ydim; j++) {
						
						dwi = distOnMap(i, j, xwin, ywin);
						hwi = neigh(dwi, rad);
						
						for (l=0; l<nbcond; l++) {
							cnt++;						
							neuron[i][j][l] = neuron[i][j][l] + alp * hwi * (data[k][l] - neuron[i][j][l]);
						}
					}
				}
			}
		}
		
		// allocate a huge chunk of memory
		affectx = (int*)malloc(nborfs * sizeof(int));
		affecty = (int*)malloc(nborfs * sizeof(int));
		
		// compute the cluster affected to each orf
		for (k=0; k<nborfs; k++) {
			mindist = 100000.0;
			// determine the distance for each gene		
			xwin = -1;
			ywin = -1;
			for (i=0; i<xdim; i++) {
				for (j=0; j<ydim; j++) {
					dist = distance(data[k], neuron[i][j], nbcond);
					cnt++;
					if (dist < mindist) {
						xwin    = i;
						ywin    = j;
						mindist = dist;
					}
				}
			}
			
			// store the coordinates of the winner cluster     
			affectx[k] = xwin;
			affecty[k] = ywin;
			
			if (verbose == 1) {
				printf("orf %s affected to (%d,%d) ..\n", names[k], xwin, ywin);
			}
		}
		/* call the function that calculates and stores the average st. dev for the genes of each cluster */
		
		tmp = estimate_tightness(data, nborfs, nbcond, xdim, ydim, affectx, affecty, names, distance_method);
				
		printf("%f\n", tmp);
		
		if(tmp < min_qual) {
			min_qual = tmp;
			
			if(min_affectx != NULL)
				free(min_affectx);
			
			if(min_affecty != NULL)
				free(min_affecty);
			
			min_affectx = affectx;
			min_affecty	= affecty;
		}
		else {				/* cleanup */
			for (i=0; i<xdim; i++) {		
				for (j=0; j<ydim; j++) {
					free(neuron[i][j]);
				}
				free(neuron[i]);
			}
			free(neuron);
			free(used);
			free(affectx);
			free(affecty);
		}
	} //end of while
	
	printf("Quality (tightness of clusters): %f\n", min_qual);
	
    /* export the clusters to outfile */	
	cnt = 0; 
	fp = fopen(outfile, "w");
	fprintf(fp, "GENE\tCLUSTER\n");
	
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			for (k=0; k<nborfs; k++) {
				if ( (min_affectx[k] == i) && (min_affecty[k] == j) ) {
					fprintf(fp, "%s\t%d\n", names[k], cnt);
				}
			}
			cnt ++;
		}
	}
	
	if (outepsmap != 0) {
		printf("Writing EPS map ...");
		if(exp_map == 1) {
			write2Dmap(outepsmap, data, nborfs, nbcond, xdim, ydim, min_affectx, min_affecty, show_shade); 
		}
		else {
			write2Dmap_general(outepsmap, data, nborfs, nbcond, xdim, ydim, min_affectx, min_affecty, show_shade); 
		}
		printf("Done.\n");
	}
	
	/* cleanup */
	fclose(fp);
	
	for (i=0; i<NBGENES; i++) {
		free(data[i]);
	}
	
	for (i=0; i<nborfs; i++) {
		free(names[i]);
	}
	
	free(data);
	//free(neuron);
	free(names);
	free(buff);
	//free(used);
	free(min_affectx);
	free(min_affecty);
	free(quality_array);
	
	return 0;
}

char* get_parameter(int argc, char** argv, char* param)
{
	int i = 0;
	while ((i < argc) && (strcmp(param, argv[i]))) 
		i++;
	if (i<argc)
		return (argv[i+1]);
	else
		return "";
}

int exist_parameter(int argc, char** argv, char* param) {
	if (strlen(get_parameter(argc, argv, param)) > 0)
		return 1;
	else
		return 0;
}

// remove the last \n, wherever it is ..
void chomp(char* s) {
	
	int j;
	
	j = strlen(s);
	while (j >= 0) {
		if (s[j] == '\n') {
			s[j] = '\0';
			return;
		}
		j--;
	}
}

// function alpha
float alpha (int iter, int length, float alp0) 
{
    return (alp0 * (length - iter) / length);
}

// function de rayon "Gaussian"
float neigh (float d, float s) 
{	
    return (exp(- (d * d) / (2.0 * s * s)));	
}

// function de diminution du rayon
float radius (int iter, int length, float rad0, float radf) 
{    
    return (radf + (rad0 - radf) * (length - iter) / length);
}

// Euclidean distance between two neurons
float distance(float* ind1, float* ind2, int n) 
{	
	int i;
	
	float sum = 0.0;
	
	for (i=0; i<n; i++) {
		sum += pow((ind2[i] - ind1[i]),2);
	}
	
	return sqrt(sum);	
}

// Pearson correlation distance
float pearson(float* data1, float* data2, int m) 
{
	
	float sum1   = 0.0;
	float sum2   = 0.0;
    
	float sum1_2 = 0.0;
	float sum2_2 = 0.0;
	
	float avg1   = 0.0;
	float avg2   = 0.0;
    
	float a, b, c;
	
	int i;
	
	int m_actual = 0;
	
	for (i=0; i<m; i++) {
		
		if (!isnan(data1[i]) && !isnan(data2[i])) {
			sum1   += data1[i];
			sum2   += data2[i];
			
			sum1_2 += data1[i] * data1[i];
			sum2_2 += data2[i] * data2[i];
			m_actual ++;
		}
	}
    
	avg1 = sum1 / m_actual;
	avg2 = sum2 / m_actual;
	
	// calc the Pearson correlation
    
	a = 0.0;
	b = 0.0;
	c = 0.0;
	
	for (i=0; i<m; i++) {
		
		if (!isnan(data1[i]) && !isnan(data2[i])) {      
			a += (data1[i] - avg1) * (data2[i] - avg2);
			b += (data1[i] - avg1) * (data1[i] - avg1);
			c += (data2[i] - avg2) * (data2[i] - avg2);
		}
	}
	
	if (fabs(a) < DBL_EPSILON)
		return 0.0;
	else    
		return a / sqrt ( b * c );
	
}


// calculate distance between two neurons
float distOnMap (int x1, int y1, int x2, int y2) 
{
    return sqrt ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

//make general 2D map
void write2Dmap_general(char* psname, float** data, int nborfs, int nbcond, int xdim, int ydim, int* affectx, int* affecty, int shade) 
{
	int i			= 0;
	int j			= 0;
	int k			= 0;
	int x			= 0;
	int	y			= 0;
	int l			= 0;
	int maxsize		= 0;
	
	int** count		= NULL;
	float*** cavg	= NULL;
	float*** cmin	= NULL;
	float*** cmax	= NULL;
	
	float maxvalue	= 0.0;
	float minvalue	= 0.0;
	
	float width		= 800 / (float)xdim;
	float height	= 553 / (float)ydim;
	float scalefactor = 0.0;
	
	FILE* fp		= NULL;
	
	count = (int**)malloc(xdim * sizeof(int*));
	for (i=0; i<xdim; i++) {
		count[i] = (int*)malloc(ydim * sizeof(int));
		for (j=0; j<ydim; j++) {
			count[i][j] = 0;
		}
	}
	
	//get maximum number of data vector per class to display the discs
	for (k=0; k<nborfs; k++) {     
		count[affectx[k]][affecty[k]] ++;
	}
	
	maxsize = 0;
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			
			if (count[i][j] > maxsize)
				maxsize = count[i][j];
		}
	}
	
	// compute the cavg profiles
	cavg = (float***)malloc(xdim * sizeof(float**));
	cmin = (float***)malloc(xdim * sizeof(float**));
	cmax = (float***)malloc(xdim * sizeof(float**));
    
	for (i=0; i<xdim; i++) {
		cavg[i] = (float**)malloc(ydim * sizeof(float*));
		cmin[i] = (float**)malloc(ydim * sizeof(float*));
		cmax[i] = (float**)malloc(ydim * sizeof(float*));
		
		for (j=0; j<ydim; j++) {
			cavg[i][j] = (float*)malloc(nbcond * sizeof(float));
			cmin[i][j] = (float*)malloc(nbcond * sizeof(float));
			cmax[i][j] = (float*)malloc(nbcond * sizeof(float));
			
			for (l=0; l<nbcond; l++) {
				cavg[i][j][l] = 0.0; 
				cmin[i][j][l] = 100000.0; 
				cmax[i][j][l] = -100000.0; 
				
			}
			
		}
	}
	
	maxvalue = FLT_MIN;
	minvalue = FLT_MAX;
	
	// add all the members of the same class
	for (k=0; k<nborfs; k++) {
		x  =  affectx[k];
		y  =  affecty[k];
		
		for (l=0; l<nbcond; l++) {
			
			// avg
			cavg[x][y][l] += data[k][l];
			
			// min for this condition
			if (cmin[x][y][l] > data[k][l]) {
				cmin[x][y][l] = data[k][l];
			}
			
			// max for this condition
			if (cmax[x][y][l] < data[k][l]) {
				cmax[x][y][l] = data[k][l];
			}
			
			// overall min
			if (minvalue > data[k][l]) {
				minvalue = data[k][l];
			}
			
			// overall max
			if (maxvalue < data[k][l]) {
				maxvalue = data[k][l];
			}
		}
	}
	
	// divide by the number of members
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			for (l=0; l<nbcond; l++) {
				cavg[i][j][l] /= count[i][j];
			}
		}
	}
	
	//printf("max value = %3.2f, min value = %3.2f\n", maxvalue, minvalue);
	
	fp = fopen(psname, "w");
	if (!fp) {
		printf("pb opening %s ..\n", psname);
		exit(0);
	}
	
	fprintf(fp, "%%!PS-Adobe-2.0 EPSF-2.0\n");
	fprintf(fp, "%%%%Title: MAP\n%%%%Creator: myself\n");
	fprintf(fp, "%%%%BoundingBox: 20 20 573 850\n");
	fprintf(fp, "%%%%Orientation: Landscape\n");
	fprintf(fp, "%%%%Pages: 1\n%%%%EndComments\n");
	fprintf(fp, "20 20 translate\n");
	
	if (minvalue > 0.0)
		minvalue = 0.0; //precaution
	
	if (-minvalue > maxvalue)
		maxvalue = -minvalue;
	
	// a few adjustements
	if (width/6.0 > height/4.0)
		scalefactor = height/4.0;
	else
		scalefactor = width/6.0;
	
	for (i=0; i<=xdim; i++) {
		fprintf(fp, "newpath\n");
		fprintf(fp, "0 %d moveto\n",   i * 800 / xdim);
		fprintf(fp, "553 %d lineto\n", i * 800 / xdim);;
		fprintf(fp, "stroke\n");
	}
	
	for (i=0; i<=ydim; i++) {
		fprintf(fp, "newpath\n");
		fprintf(fp, "%d 0 moveto\n", i * 553 / ydim);
		fprintf(fp, "%d 800 lineto\n", i * 553 / ydim);
		fprintf(fp, "stroke\n");
	}
	
	int cluster_num = 1;
	
	for (j=0; j<ydim; j++) {
		for (i=0; i<xdim; i++) {
			
			float colshade;
			if(shade == 1) {
				colshade = (float)count[i][j]/(float)maxsize;
			}
			else {
				colshade = 0.5;
			}
			
			fprintf(fp, "gsave\n");
			fprintf(fp, "0 %3.2f 0 setrgbcolor\n", colshade);
			fprintf(fp, "%3.2f %3.2f translate\n", height / 2.0 + j * 553 / (float)ydim, width / 2.0 + i * 800 / (float)xdim);
			fprintf(fp, "newpath\n");
			fprintf(fp, "0 0 %3.2f 0 360 arc\n", count[i][j] * scalefactor / maxsize);
			fprintf(fp, "fill\n");
			
			fprintf(fp, "0 0 0.5 setrgbcolor\n");
			fprintf(fp, "/Times-Roman findfont\n");
			fprintf(fp, "%3.2f scalefont\n", 60/(float)xdim);
			fprintf(fp, "setfont\n");
			fprintf(fp, "newpath\n");
			
			fprintf(fp, "-%3.2f -%3.2f moveto\n", height/2.4, width/2.2);
			
			fprintf(fp, "(Cluster %d: %d) show\n", cluster_num, count[i][j]);
			
			fprintf(fp, "grestore\n");
			cluster_num++;
		}
	}
	
	fprintf(fp,  "showpage\n");
	
	fclose(fp);
	
	/* cleanup */
	
	for (i=0; i<xdim; i++) {
		free(count[i]);
	}
	
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			free(cavg[i][j]);
			free(cmin[i][j]);
			free(cmax[i][j]);
		}
		free(cavg[i]);
		free(cmin[i]);
		free(cmax[i]);
	}
	
	free(count);
	free(cavg);
	free(cmin);
	free(cmax);
}

//make 2D map for expression data
void write2Dmap(char* psname, float** data, int nborfs, int nbcond, int xdim, int ydim, int* affectx, int* affecty, int shade) 
{
	int i			= 0;
	int j			= 0;
	int k			= 0;
	int x			= 0;
	int	y			= 0;
	int l			= 0;
	int maxsize		= 0;
	
	int** count		= NULL;
	float*** cavg	= NULL;
	float*** cmin	= NULL;
	float*** cmax	= NULL;
	
	float maxvalue	= 0.0;
	float minvalue	= 0.0;
	
	float width		= 800 / (float)xdim;
	float height	= 553 / (float)ydim;
	float scalefactor = 0.0;
	
	FILE* fp		= NULL;
	
	count = (int**)malloc(xdim * sizeof(int*));
	for (i=0; i<xdim; i++) {
		count[i] = (int*)malloc(ydim * sizeof(int));
		for (j=0; j<ydim; j++) {
			count[i][j] = 0;
		}
	}
	
	//get maximum number of data vector per class to display the discs
	for (k=0; k<nborfs; k++) {     
		count[affectx[k]][affecty[k]] ++;
	}
	
	maxsize = 0;
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			
			if (count[i][j] > maxsize)
				maxsize = count[i][j];
		}
	}
	
	// compute the cavg profiles
	cavg = (float***)malloc(xdim * sizeof(float**));
	cmin = (float***)malloc(xdim * sizeof(float**));
	cmax = (float***)malloc(xdim * sizeof(float**));
    
	for (i=0; i<xdim; i++) {
		cavg[i] = (float**)malloc(ydim * sizeof(float*));
		cmin[i] = (float**)malloc(ydim * sizeof(float*));
		cmax[i] = (float**)malloc(ydim * sizeof(float*));
		
		for (j=0; j<ydim; j++) {
			cavg[i][j] = (float*)malloc(nbcond * sizeof(float));
			cmin[i][j] = (float*)malloc(nbcond * sizeof(float));
			cmax[i][j] = (float*)malloc(nbcond * sizeof(float));
			
			for (l=0; l<nbcond; l++) {
				cavg[i][j][l] = 0.0; 
				cmin[i][j][l] = 100000.0; 
				cmax[i][j][l] = -100000.0; 
				
			}
			
		}
	}
	
	maxvalue = FLT_MIN;
	minvalue = FLT_MAX;
	
	// add all the members of the same class
	for (k=0; k<nborfs; k++) {
		x  =  affectx[k];
		y  =  affecty[k];
		
		for (l=0; l<nbcond; l++) {
			
			// avg
			cavg[x][y][l] += data[k][l];
			
			// min for this condition
			if (cmin[x][y][l] > data[k][l]) {
				cmin[x][y][l] = data[k][l];
			}
			
			// max for this condition
			if (cmax[x][y][l] < data[k][l]) {
				cmax[x][y][l] = data[k][l];
			}
			
			// overall min
			if (minvalue > data[k][l]) {
				minvalue = data[k][l];
			}
			
			// overall max
			if (maxvalue < data[k][l]) {
				maxvalue = data[k][l];
			}
		}
	}
	
	// divide by the number of members
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			for (l=0; l<nbcond; l++) {
				cavg[i][j][l] /= count[i][j];
			}
		}
	}
	
	//printf("max value = %3.2f, min value = %3.2f\n", maxvalue, minvalue);
	
	fp = fopen(psname, "w");
	if (!fp) {
		printf("pb opening %s ..\n", psname);
		exit(0);
	}
	
	fprintf(fp, "%%!PS-Adobe-2.0 EPSF-2.0\n");
	fprintf(fp, "%%%%Title: MAP\n%%%%Creator: myself\n");
	fprintf(fp, "%%%%BoundingBox: 20 20 573 850\n");
	fprintf(fp, "%%%%Orientation: Landscape\n");
	fprintf(fp, "%%%%Pages: 1\n%%%%EndComments\n");
	fprintf(fp, "20 20 translate\n");
	
	if (minvalue > 0.0)
		minvalue = 0.0; //precaution
	
	if (-minvalue > maxvalue)
		maxvalue = -minvalue;
	
	// a few adjustements
	if (width/6.0 > height/4.0)
		scalefactor = height/4.0;
	else
		scalefactor = width/6.0;
	
	for (i=0; i<=xdim; i++) {
		fprintf(fp, "newpath\n");
		fprintf(fp, "0 %d moveto\n",   i * 800 / xdim);
		fprintf(fp, "553 %d lineto\n", i * 800 / xdim);
		fprintf(fp, "stroke\n");
	}
	
	for (i=0; i<=ydim; i++) {
		fprintf(fp, "newpath\n");
		fprintf(fp, "%d 0 moveto\n", i * 553 / ydim);
		fprintf(fp, "%d 800 lineto\n", i * 553 / ydim);
		fprintf(fp, "stroke\n");
	}
	
	int cluster_num = 1;
	
	for (j=0; j<ydim; j++) {
		for (i=0; i<xdim; i++) {
			
			float colshade;
			if(shade == 1) {
				colshade = (float)count[i][j]/(float)maxsize;
			}
			else {
				colshade = 0.5;
			}
			
			//float colshade = 0.5;
			//float colshade = (float)count[i][j]/(float)maxsize;
			
			fprintf(fp, "gsave\n");
			fprintf(fp, "0 %3.2f 0 setrgbcolor\n", colshade);
			fprintf(fp, "%3.2f %3.2f translate\n", height / 2.0 + j * 553 / (float)ydim, width / 2.0 + i * 800 / (float)xdim);
			fprintf(fp, "newpath\n");
			fprintf(fp, "0 0 %3.2f 0 360 arc\n", count[i][j] * scalefactor / maxsize);
			fprintf(fp, "fill\n");
			
			fprintf(fp, "0 0 0.5 setrgbcolor\n");
			fprintf(fp, "/Times-Roman findfont\n");
			fprintf(fp, "%3.2f scalefont\n", 60/(float)xdim);
			fprintf(fp, "setfont\n");
			fprintf(fp, "newpath\n");
			
			fprintf(fp, "-%3.2f -%3.2f moveto\n", height/2.4, width/2.2);
			
			fprintf(fp, "(Cluster %d: %d) show\n", cluster_num, count[i][j]);
			
			fprintf(fp, "grestore\n");
			cluster_num++;
			
			if (count[i][j] > 0) {
				
				// AVG
				fprintf(fp, "gsave\n");
				fprintf(fp, "%3.2f setlinewidth\n", (ydim * (maxvalue - minvalue)) / 553);
				fprintf(fp, "%3.2f %3.2f translate\n",  height/2.0 + j * 553 / (float)ydim,  i * 800 / (float)xdim);
				fprintf(fp, "%3.2f %3.2f scale\n", (553 / (ydim * (maxvalue - minvalue))) / 2.0, 800 / ((float)xdim * (nbcond-1)));
				//fprintf(fp, "0 0 1 setrgbcolor\n");
				fprintf(fp,  "newpath\n");
				
				fprintf(fp,  "%3.2f 0 moveto\n", maxvalue - cavg[i][j][0]);
				
				for (k=1; k<nbcond; k++)
					fprintf(fp, "%3.2f %d lineto\n", maxvalue - cavg[i][j][k], k);
				fprintf(fp,  "stroke\n");
				fprintf(fp,  "grestore\n");
				
				
				// MAX
				fprintf(fp, "gsave\n");
				fprintf(fp, "1 0 0 setrgbcolor\n");
				
				fprintf(fp, "%3.2f setlinewidth\n", (ydim * (maxvalue - minvalue)) / 553);
				fprintf(fp, "%3.2f %3.2f translate\n",  height/2.0 + j * 553 / (float)ydim,  i * 800 / (float)xdim);
				fprintf(fp, "%3.2f %3.2f scale\n", (553 / (ydim * (maxvalue - minvalue))) / 2.0, 800 / ((float)xdim * (nbcond-1)));
				fprintf(fp,  "newpath\n");
				
				fprintf(fp,  "%3.2f 0 moveto\n", maxvalue - cmax[i][j][0]);
				
				for (k=1; k<nbcond; k++)
					fprintf(fp, "%3.2f %d lineto\n", maxvalue - cmax[i][j][k], k);
				fprintf(fp,  "stroke\n");
				fprintf(fp,  "grestore\n");
				
				
				// MIN
				fprintf(fp, "gsave\n");
				fprintf(fp, "0 1 0 setrgbcolor\n");
				
				fprintf(fp, "%3.2f setlinewidth\n", (ydim * (maxvalue - minvalue)) / 553);
				fprintf(fp, "%3.2f %3.2f translate\n",  height/2.0 + j * 553 / (float)ydim,  i * 800 / (float)xdim);
				fprintf(fp, "%3.2f %3.2f scale\n", (553 / (ydim * (maxvalue - minvalue))) / 2.0, 800 / ((float)xdim * (nbcond-1)));
				fprintf(fp,  "newpath\n");
				
				fprintf(fp,  "%3.2f 0 moveto\n", maxvalue - cmin[i][j][0]);
				
				for (k=1; k<nbcond; k++)
					fprintf(fp, "%3.2f %d lineto\n", maxvalue - cmin[i][j][k], k);
				fprintf(fp,  "stroke\n");
				fprintf(fp,  "grestore\n");
			}
		}
	}
	
	fprintf(fp,  "showpage\n");
	
	fclose(fp);
	
	/* cleanup */
	
	for (i=0; i<xdim; i++) {
		free(count[i]);
	}
	
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			free(cavg[i][j]);
			free(cmin[i][j]);
			free(cmax[i][j]);
		}
		free(cavg[i]);
		free(cmin[i]);
		free(cmax[i]);
	}
	
	free(count);
	free(cavg);
	free(cmin);
	free(cmax);
}

float estimate_tightness(float** data, int nborfs, int nbcond, int xdim, int ydim, int* affectx, int* affecty, char** names, int distance_method) {
	
	int i			= 0;
	int j			= 0;
	int k			= 0;
	int l			= 0;
	
	int** count				= NULL;
	float*** cavg			= NULL;
	float sum_sq_dist		= 0.0;
	
	float root_sq_dist		= 0.0;
	float avg_sq_dist		= 0.0;
	float total_quality		= 0.0;
	
	/* get the number of genes per cluster */
	count = (int**)malloc(xdim * sizeof(int*));
	for (i=0; i<xdim; i++) {
		count[i] = (int*)malloc(ydim * sizeof(int));
		for (j=0; j<ydim; j++) {
			count[i][j] = 0;
		}
	}	
	
	for (k=0; k<nborfs; k++) {     
		count[affectx[k]][affecty[k]] ++;
		
	}
	
	/* compute the cavg profiles */
	cavg = (float***)calloc(xdim, sizeof(float**));
	
	for (i=0; i<xdim; i++) {
		cavg[i] = (float**)calloc(ydim, sizeof(float*));
		for (j=0; j<ydim; j++) {
			cavg[i][j] = (float*)calloc(nbcond, sizeof(float));			
		}
	}
	
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			
			for (k=0; k<nborfs; k++) {
				if ( (affectx[k] == i) && (affecty[k] == j) ) {
					
					for (l=0; l<nbcond; l++) {			
						cavg[i][j][l] += data[k][l];
					}
				}
			}
		}
	}	
	
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			for (l=0; l<nbcond; l++) {
				cavg[i][j][l] /= count[i][j];
			}
		}
	}
	
	/* estimate euclidean distance of each gene from average profile */
	/* (average square distance for the whole cluster) */
	
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			
			avg_sq_dist = 0.0;
			//printf("Cluster(%d %d)\n", i, j);
			
			for (k=0; k<nborfs; k++) {
				if ( (affectx[k] == i) && (affecty[k] == j) ) {
					//printf("%s\t", names[k]);
					
					sum_sq_dist		= 0.0;
					root_sq_dist	= 0.0;
					
					/* choose the distance_method */
					if(distance_method)
						root_sq_dist = distance(cavg[i][j], data[k], nbcond);
					else
						root_sq_dist = (1.0 - pearson(cavg[i][j], data[k], nbcond));

					/* printf("E_dist:%f\t", distance(cavg[i][j], data[k], nbcond));
					printf("Pearson_dist:%f\t", (1.0 - pearson(cavg[i][j], data[k], nbcond)));					
					printf("root\t%f\n", root_sq_dist); */
					
					avg_sq_dist += root_sq_dist;
					
					//printf("\n");
				}
			}
			
			if(count[i][j] != 0) {
				avg_sq_dist /= count[i][j];
				total_quality += avg_sq_dist;
			}
			
			//printf("AVG(%d, %d)\t%5.4f\n", i, j, avg_sq_dist);
			//printf("COUNT\t%d\n", count[i][j]);
		}
	}
	
	/* cleanup */	
	for (i=0; i<xdim; i++) {
		free(count[i]);
	}
	
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			free(cavg[i][j]);
		}
		free(cavg[i]);
	}
	
	free(cavg);	
	free(count);
	
	return (total_quality/(xdim*ydim));
}
