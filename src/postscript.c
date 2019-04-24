#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <sys/stat.h>
#include <float.h>

#include "postscript.h"

void colMeans2Dplot(char* psname, float** data, int nbpeaks, int nbbins, char* xlabel, char* ylabel, int binsize) 
{
	int i, k, l;
	float* cavg;
	float* cmin;
	float* cmax;
	
	float maxvalue;
	float minvalue;

	float avgmaxvalue;
	float avgminvalue;
	
	float width  = 500;
	float height = 400;
	float scalefactor;
	
	FILE* fp;
	
	if(xlabel == NULL) {
		xlabel = "x axis";
	}
	
	if(ylabel == NULL) {
		ylabel = "y axis";
	}

	//printf("NB PEAKS: %d\n", nbpeaks);
	//printf("NB BINS: %d\n", nbbins);
	//printf("WS: %d\n", binsize);
	
	// compute the cavg profiles
	cavg = (float*)malloc(nbbins * sizeof(float));
	cmin = (float*)malloc(nbbins * sizeof(float));
	cmax = (float*)malloc(nbbins * sizeof(float));
    
	for (i=0; i<nbbins; i++) {
		cavg[i] = 0.0;
		cmin[i] = 100000.0;
		cmax[i] = -100000.0;		
	}
	
	maxvalue = FLT_MIN;
	minvalue = FLT_MAX;

	avgmaxvalue = FLT_MIN;
	avgminvalue = FLT_MAX;
	
	// add all the members of the same class
	for (k=0; k<nbbins; k++) {
		
		for (l=0; l<nbpeaks; l++) {
			
			// avg
			cavg[k] += data[l][k];
			
			// overall min
			if (minvalue > data[l][k]) {
				minvalue = data[l][k];
			}
			
			// overall max
			if (maxvalue < data[l][k]) {
				maxvalue = data[l][k];
			}
		}
		
		cavg[k] /= nbpeaks;
		
		// avg min
		if (avgminvalue > cavg[k]) {
			avgminvalue = cavg[k];
		}
		
		// avg max
		if (avgmaxvalue < cavg[k]) {
			avgmaxvalue = cavg[k];
		}
		
		/*printf("k=%d : %3.4f\t", k, cavg[k]);
		 
		 cavg[k] /= nbpeaks;
		 
		 printf("%3.4f\t %d\n", cavg[k], nbpeaks);*/
	}

	fp = fopen(psname, "w");
	if (!fp) {
		printf("pb opening %s ..\n", psname);
		exit(0);
	}
	
	fprintf(fp, "%%!PS-Adobe-2.0 EPSF-2.0\n");
	fprintf(fp, "%%%%Title: MAP\n%%%%Creator: myself\n");
	fprintf(fp, "%%%%BoundingBox: 0 0 %f %f\n", width, height);
	fprintf(fp, "%%%%Orientation: Landscape\n");
	fprintf(fp, "%%%%Pages: 1\n%%%%EndComments\n");
	
	if (minvalue > 0.0)
		minvalue = 0.0; //precaution
	
	if (-minvalue > maxvalue)
		maxvalue = -minvalue;
	
	// a few adjustements
	if (width/6.0 > height/4.0)
		scalefactor = height/4.0;
	else
		scalefactor = width/6.0;
		
	//make axes
	fprintf(fp, "1.5 setlinewidth\n");
	
	fprintf(fp, "newpath\n");
	
	fprintf(fp, "45 45 moveto\n");
	fprintf(fp, "45 380 lineto\n");
	
	fprintf(fp, "45 45 moveto\n");
	fprintf(fp, "480 45 lineto\n");
	
	fprintf(fp, "0 0 0 setrgbcolor\n");	
	fprintf(fp, "stroke\n");
	
	
	//write axes labels
	fprintf(fp, "0 0 0 setrgbcolor\n");	
	
	fprintf(fp, "/Times-Roman findfont\n");	
	fprintf(fp, "15 scalefont\n");	
	fprintf(fp, "setfont \n");
	
	//x axis
	fprintf(fp, "newpath\n");
	fprintf(fp, "190 10 moveto\n");
	fprintf(fp, "(%s) show\n", xlabel);	
	
	//y axis
	fprintf(fp, "gsave\n");
	fprintf(fp, "newpath\n");
	fprintf(fp, "15 180 moveto\n");
	fprintf(fp, "90 rotate\n");
	fprintf(fp, "(%s) show\n", ylabel);
	fprintf(fp,  "grestore\n");
	
	//y axis min and max
	fprintf(fp, "/Times-Roman findfont\n");	
	fprintf(fp, "10 scalefont\n");	
	fprintf(fp, "setfont \n");
	
	fprintf(fp, "newpath\n");
	fprintf(fp, "20 40 moveto\n");
	fprintf(fp, "(0.0) show\n");	
		
	//print the scaled ymin/ymax labels on the axis
	fprintf(fp, "newpath\n");
	fprintf(fp, "20 %3.2f moveto\n", (height/avgmaxvalue)*1.85*avgminvalue);
	fprintf(fp, "(%3.2f) show\n", avgminvalue);
	
	fprintf(fp, "newpath\n");
	fprintf(fp, "20 %3.2f moveto\n", (height/avgmaxvalue)*0.90*avgmaxvalue);
	fprintf(fp, "(%3.2f) show\n", avgmaxvalue);

	//x axis min and max
	fprintf(fp, "newpath\n");
	fprintf(fp, "45 30 moveto\n");
	fprintf(fp, "(-%d) show\n", nbbins*binsize);

	fprintf(fp, "newpath\n");
	fprintf(fp, "470 30 moveto\n");
	fprintf(fp, "(%d) show\n", nbbins*binsize);
	
	fprintf(fp, "newpath\n");
	fprintf(fp, "260 30 moveto\n");
	fprintf(fp, "(0) show\n");
	
	//make plot
	fprintf(fp, "gsave\n");
	fprintf(fp, "0 0 1 setrgbcolor\n");
	fprintf(fp, "0.0000001 setlinewidth\n");
	
	fprintf(fp,  "newpath\n");
	fprintf(fp, "%3.2f %3.2f translate\n", width/10, height/13);
	fprintf(fp, "%3.2f %3.2f scale\n", (binsize*nbbins)/1000 + 0.12, (height/avgmaxvalue)*0.85);
	
	fprintf(fp, "0 %3.2f moveto\n", cavg[0]);
	
	for (k=1; k<nbbins; k++) {
		
		fprintf(fp, "%d %3.2f lineto\n", k, cavg[k]);
		
	}
	
	fprintf(fp,  "stroke\n");
	fprintf(fp,  "grestore\n");
	fprintf(fp,  "showpage\n");
	
	fclose(fp);
}


void colMeans2Ddoubleplot(char* psname, float** data1, float** data2, int nbpeaks, int nbbins, char* xlabel, char* ylabel, int binsize, char* legend1, char* legend2) 
{
	int i, k, l;
	float* cavg1;
	float* cavg2;
	float* cmin;
	float* cmax;
	
	float maxvalue;
	float minvalue;
	
	float avgmaxvalue1;
	float avgminvalue1;
	
	float avgmaxvalue2;
	float avgminvalue2;
	
	float width  = 500;
	float height = 400;
	float scalefactor;
	
	FILE* fp;
	
	if(xlabel == NULL) {
		xlabel = "x axis";
	}
	
	if(ylabel == NULL) {
		ylabel = "y axis";
	}
	
	/*printf("NB PEAKS: %d\n", nbpeaks);
	printf("NB BINS: %d\n", nbbins);
	printf("WS: %d\n", binsize);*/
	
	// compute the cavg profiles
	cavg1 = (float*)malloc(nbbins * sizeof(float));
	cavg2 = (float*)malloc(nbbins * sizeof(float));
	cmin = (float*)malloc(nbbins * sizeof(float));
	cmax = (float*)malloc(nbbins * sizeof(float));
    
	for (i=0; i<nbbins; i++) {
		cavg1[i] = 0.0;
		cavg2[i] = 0.0;
		cmin[i] = 100000.0;
		cmax[i] = -100000.0;		
	}
	
	maxvalue = FLT_MIN;
	minvalue = FLT_MAX;
	
	avgmaxvalue1 = FLT_MIN;
	avgminvalue1 = FLT_MAX;

	avgmaxvalue2 = FLT_MIN;
	avgminvalue2 = FLT_MAX;
	
	// add all the members of the same class
	for (k=0; k<nbbins; k++) {
		
		for (l=0; l<nbpeaks; l++) {
			
			// avg
			cavg1[k] += data1[l][k];
			cavg2[k] += data2[l][k];

			
			// overall min
			if (minvalue > data1[l][k]) {
				minvalue = data1[l][k];
			}
			
			// overall max
			if (maxvalue < data1[l][k]) {
				maxvalue = data1[l][k];
			}
		}
		
		cavg1[k] /= nbpeaks;
		cavg2[k] /= nbpeaks;

		
		// avg min
		if (avgminvalue1 > cavg1[k]) {
			avgminvalue1 = cavg1[k];
		}
		
		// avg max
		if (avgmaxvalue1 < cavg1[k]) {
			avgmaxvalue1 = cavg1[k];
		}
		
		if (avgminvalue2 > cavg2[k]) {
			avgminvalue2 = cavg2[k];
		}
		
		// avg max
		if (avgmaxvalue2 < cavg2[k]) {
			avgmaxvalue2 = cavg2[k];
		}
		
		/*printf("k=%d : %3.4f\t", k, cavg[k]);
		 
		 cavg[k] /= nbpeaks;
		 
		 printf("%3.4f\t %d\n", cavg[k], nbpeaks);*/
	}
	
	fp = fopen(psname, "w");
	if (!fp) {
		printf("pb opening %s ..\n", psname);
		exit(0);
	}
	
	fprintf(fp, "%%!PS-Adobe-2.0 EPSF-2.0\n");
	fprintf(fp, "%%%%Title: MAP\n%%%%Creator: myself\n");
	fprintf(fp, "%%%%BoundingBox: 0 0 %f %f\n", width, height);
	fprintf(fp, "%%%%Orientation: Landscape\n");
	fprintf(fp, "%%%%Pages: 1\n%%%%EndComments\n");
	
	if (minvalue > 0.0)
		minvalue = 0.0; //precaution
	
	if (-minvalue > maxvalue)
		maxvalue = -minvalue;
	
	// a few adjustements
	if (width/6.0 > height/4.0)
		scalefactor = height/4.0;
	else
		scalefactor = width/6.0;
		
	//make axes
	fprintf(fp, "1.5 setlinewidth\n");
	
	fprintf(fp, "newpath\n");
	
	fprintf(fp, "45 45 moveto\n");
	fprintf(fp, "45 380 lineto\n");
	
	fprintf(fp, "45 45 moveto\n");
	fprintf(fp, "480 45 lineto\n");
	
	fprintf(fp, "0 0 0 setrgbcolor\n");	
	fprintf(fp, "stroke\n");
	
	
	//write axes labels
	fprintf(fp, "0 0 0 setrgbcolor\n");	
	
	fprintf(fp, "/Times-Roman findfont\n");	
	fprintf(fp, "15 scalefont\n");	
	fprintf(fp, "setfont \n");
	
	//x axis
	fprintf(fp, "newpath\n");
	fprintf(fp, "190 10 moveto\n");
	fprintf(fp, "(%s) show\n", xlabel);	
	
	//y axis
	fprintf(fp, "gsave\n");
	fprintf(fp, "newpath\n");
	fprintf(fp, "15 180 moveto\n");
	fprintf(fp, "90 rotate\n");
	fprintf(fp, "(%s) show\n", ylabel);
	fprintf(fp,  "grestore\n");
	
	//y axis min and max
	fprintf(fp, "/Times-Roman findfont\n");	
	fprintf(fp, "10 scalefont\n");	
	fprintf(fp, "setfont \n");
	
	fprintf(fp, "newpath\n");
	fprintf(fp, "20 40 moveto\n");
	fprintf(fp, "(0.0) show\n");	
	
	//print the scaled ymin/ymax labels on the axis
	fprintf(fp, "newpath\n");
	fprintf(fp, "20 %3.2f moveto\n", (height/avgmaxvalue1)*0.90*avgminvalue1);
	fprintf(fp, "(%3.2f) show\n", avgminvalue1);
	
	fprintf(fp, "newpath\n");
	fprintf(fp, "20 %3.2f moveto\n", (height/avgmaxvalue1)*0.90*avgmaxvalue1);
	fprintf(fp, "(%3.2f) show\n", avgmaxvalue1);
	
	//x axis min and max
	fprintf(fp, "newpath\n");
	fprintf(fp, "45 30 moveto\n");
	fprintf(fp, "(-%d) show\n", nbbins*binsize);
	
	fprintf(fp, "newpath\n");
	fprintf(fp, "470 30 moveto\n");
	fprintf(fp, "(%d) show\n", nbbins*binsize);
	
	fprintf(fp, "newpath\n");
	fprintf(fp, "260 30 moveto\n");
	fprintf(fp, "(0) show\n");
	
	//make first plot
	fprintf(fp, "gsave\n");
	fprintf(fp, "1 0 0 setrgbcolor\n");
	fprintf(fp, "0.0000001 setlinewidth\n");
	
	fprintf(fp,  "newpath\n");
	fprintf(fp, "%3.2f %3.2f translate\n", width/10, height/20);
	fprintf(fp, "%3.2f %3.2f scale\n", (binsize*nbbins)/1000 + 0.12, (height/avgmaxvalue1)*0.85);
	
	fprintf(fp, "0 %3.4f moveto\n", cavg1[0]);
	
	for (k=1; k<nbbins; k++) {
		
		fprintf(fp, "%d %3.4f lineto\n", k, cavg1[k]);
		
	}
	
	fprintf(fp,  "stroke\n");
	fprintf(fp,  "grestore\n");
	
	
	//make second plot
	fprintf(fp, "gsave\n");
	fprintf(fp, "0 1 0 setrgbcolor\n");
	fprintf(fp, "0.0000001 setlinewidth\n");
	
	fprintf(fp,  "newpath\n");
	fprintf(fp, "%3.2f %3.2f translate\n", width/10, height/20);
	fprintf(fp, "%3.2f %3.2f scale\n", (binsize*nbbins)/1000 + 0.12, (height/avgmaxvalue2)*0.85);
	
	fprintf(fp, "0 %3.4f moveto\n", cavg2[0]);
	
	for (k=1; k<nbbins; k++) {
		
		fprintf(fp, "%d %3.4f lineto\n", k, cavg2[k]);
		
	}
	
	fprintf(fp,  "stroke\n");
	fprintf(fp,  "grestore\n");
	
	//write legends	
	
	//make frame for legends
	//fprintf(fp, "48 50 moveto\n");
	//fprintf(fp, "0 -50 rlineto\n");
	//fprintf(fp, "50 0 rlineto\n");
	//fprintf(fp, "0 50 rlineto\n");
	//fprintf(fp, "-50 0 rlineto\n");
	//fprintf(fp, "closepath\n");
	
	fprintf(fp, "/Times-Roman findfont\n");	
	fprintf(fp, "15 scalefont\n");	
	fprintf(fp, "setfont \n");
	
	//legend1
	fprintf(fp, "1 0 0 setrgbcolor\n");	
	fprintf(fp, "newpath\n");
	fprintf(fp, "49 380 moveto\n");
	fprintf(fp, "(%s) show\n", legend1);	
	
	//legend2
	fprintf(fp, "0 1 0 setrgbcolor\n");	
	fprintf(fp, "newpath\n");
	fprintf(fp, "49 360 moveto\n");
	fprintf(fp, "(%s) show\n", legend2);	
	
	//print page
	fprintf(fp,  "showpage\n");
	
	/* cleanup */
	fclose(fp);
	free(cavg1);
	free(cavg2);
	free(cmin);
	free(cmax);
}