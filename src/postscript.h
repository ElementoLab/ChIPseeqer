#ifndef POSTSCRIPT_H
#define POSTSCRIPT_H

void colMeans2Dplot(char* psname, float** data, int nbpeaks, int nbbins, char* xlabel, char* ylabel, int binsize);
void colMeans2Ddoubleplot(char* psname, float** data1, float** data2, int nbpeaks, int nbbins, char* xlabel, char* ylabel, int binsize, char* legend1, char* legend2);
 
#endif