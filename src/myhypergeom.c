/*
 * This program computes the hypergeometric distribution probability 
 * It takes as input:
 *	N	(e.g., all genes)
 *	s1	(e.g., genes with a concrete pathway)
 *	s2	(e.g., genes with peaks)
 *	i	(e.g., genes with peaks and pathway)
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

char* get_parameter(int argc, char** argv, char* param);

double	cumhyper(int i, int s1, int s2, int N); 
double	factln(int n);
double	gammln(double xx);
void	nrerror(char str[]);
double	hypergeom(int i, int s1, int s2, int N);

int		verbose = 0;

int main(int argc, char** argv) {
	
	int s1		= atoi(get_parameter(argc, argv, "-s1"));
	int s2		= atoi(get_parameter(argc, argv, "-s2"));
	int i		= atoi(get_parameter(argc, argv, "-i"));
	int N		= atoi(get_parameter(argc, argv, "-N"));
	
	double p	= cumhyper(i, s1, s2, N);
	
	//printf("i=%d, s1=%d, s2=%d, N=%d, p=%5.4e, log=%5.4f\n", i, s1, s2, N, p, log(p));
	printf("p=%5.4e, log=%5.4f\n", p, log(p));
	
	return 0;
}

int fillOccurrenceTable(int* table, char* scoreFile, double t) 
{
	FILE* fp;
	int	n		= 0; 
	float s		= 0.0;
	int d		= 0;
	
	fp = fopen(scoreFile, "r");
	
	while (1) {
		fscanf(fp, "%d\t%f\t%d\n", &n, &s, &d);
		
		if (s > t) {
			table[n]++;
		} else {
			break; // no need to go further
		}		
		if (feof(fp)) {
			break;
		}
	}
	
	fclose(fp);
	return 0;
}

double cumhyper(int i, int s1, int s2, int N) 
{
	
	int min		= (s1<s2?s1:s2);
	double prod = 0.0;
    int a		= 0;
    double tmp	= 0.0;
	
	for (a=i; a<=min; a++) {
		tmp = (double)hypergeom(a, s1, s2, N);
		prod += (double)hypergeom(a, s1, s2, N);
	}
	return prod;
}

void nrerror(char str[]) {
	printf("%s\n", str);
}

double hypergeom(int i, int s1, int s2, int N) {
    double factln(int n);
	
    return exp(factln(s1) + factln(N - s1) + factln(s2) + factln(N - s2) - 
			   factln(i) - factln(N) - factln(s1 - i) - factln(s2 - i) 
			   - factln(N - s1 - s2 + i));
    
}

double factln(int n)
//Returns ln(n!).
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	static double a[101];                                        
	//A static array is automatically initialized to zero.
	
	if (n < 0) nrerror("Negative factorial in routine factln");
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));                                    
	//In range of table.
	else return gammln(n+1.0);                                  
	//Out of range of table.
}


double gammln(double xx)
//Returns the value ln[ (xx)] for xx > 0.
{
	//Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
	//accuracy is good enough.
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
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
