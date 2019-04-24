#ifdef BNS
#include <pcre/pcre.h>
#else
#include <pcre.h>
#endif

void findSites_micombine(char* r, char* seq, unsigned short* positions, int* np, unsigned short* orientations, int* no, unsigned short maxnbsites, int singlestrand, int store_pos, int store_ori); 

void findSites(char* r, char* seq, int* positions, int* np, int* orientations, int* no, int maxnbsites, int singlestrand, int store_pos_ori); 
int re_matches(char* r, char* seq, int singlestrand); 
int raw_re_matches(pcre* re, char* seq, char* seq_c, int singlestrand); 
