//
// CS
//
#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)
#define _GNU_SOURCE

#ifndef SAM_H
#include "lib/third_party/samtools/sam.h"
#include "lib/third_party/samtools/bam.h"
#define SAM_H
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <search.h>
#include <ctype.h>


#include "hashtable.h"
#include "sequences.h"
#include "statistics.h"
#include "dataio.h"
#include "set.h"
#include "readio.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include "assert.h"


int formatToId(char* formattext)
{
	int format = -1;
	int i;
	for (i=0; i<strlen(formattext); i++)
		formattext[i] = tolower(formattext[i]);
	if (strcmp(formattext, "eland") == 0)
		format = ELAND;
	else if (strcmp(formattext, "mit") == 0)
		format = MIT;
	else if (strcmp(formattext, "bed") == 0)
		format = BED;
	else if (strcmp(formattext, "sam") == 0)
		format = SAM; 
	else if (strcmp(formattext, "bowtiesam") == 0)
		format = BOWTIESAM; 
	else if (strcmp(formattext, "bam") == 0)
		format = BAM;
	else {
		die("Format unrecognized\n");
	}
	return format;
}


int readMapData(char* file, unsigned char** map, int t)
{
	
	FILE* f2;
	char* buff = 0;
	
	int   mynmax = 100000;
	int   score = 0;
	int   i;
	int mapb = 0;
	f2 = fopen(file, "r");
	if (!f2) {
		printf("Cannot open file %s\n", file);
		exit(1);
	}
	
	buff  = (char*)malloc(mynmax * sizeof(char));
	
	i = 0;
	while (fgets(buff, mynmax, f2) != 0) {
		chomp(buff);
		score = atoi(buff);	  
		if (score >= t) {
			set_entry_in_binarized_array(*map, i);
			mapb++;
		} 
		i++;
	}
	
	return mapb;
}



long int CountAlignedNucleotides(char* file, int format)
{
	
  //int   numreads = 0;
	MappedRead     r;
	ReadI          ri;
	long int   numnt = 0;
	//int fraglen = 0;
	
	if (!ReadIopen(&ri, file, format)) {
		printf("Cannot open file %s\n", file);
	}
	
	while (nextRead(&ri, &r)) {
		if (r.uniqmap == 1) {
			numnt += r.lenread;
		}
	}
	ReadIclose(&ri);
	return numnt;
}

//
// count all reads ... meaning does not correct for clonal reads
//
int CountReads(char* file, int format)
{
	int   numreads = 0;
	MappedRead     r;
	ReadI          ri;
	
	if (!ReadIopen(&ri, file, format)) {
		printf("Cannot open file %s\n", file);
	}
	
	while (nextRead(&ri, &r)) {
		if (r.uniqmap == 1) 
			numreads++;
	}
	ReadIclose(&ri);
	
	
	return numreads;
}



//
// count all reads ... meaning does not correct for clonal reads
//
int CountNonClonalReads(char* file, int format, int chrlen)
{
  int   numreads = 0;
  MappedRead     r;
  ReadI          ri;
  unsigned char* a_unqreads_R   = 0;
  unsigned char* a_unqreads_F   = 0;
  int   pos     = 0;
  int   st	= 0;

  a_unqreads_R = create_binarized_array(chrlen);
  a_unqreads_F = create_binarized_array(chrlen);
	
  if (!ReadIopen(&ri, file, format)) {
    printf("Cannot open file %s\n", file);
  }
  
  while (nextRead(&ri, &r)) {
    
    if (r.uniqmap == 1) {
    
      st      = r.st;
      pos     = r.pos;
  
      if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 0)) ||
	  ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 0))) {
	// no read at that position !
	
	numreads++;  // add 1
	
	// now there is a read
	if (st == 1)
	  set_entry_in_binarized_array(a_unqreads_F, pos);
	else
	  set_entry_in_binarized_array(a_unqreads_R, pos);
	
      } // if there is no read 
      
    } // if unique
  }
  ReadIclose(&ri);
						
  free(a_unqreads_R);
  free(a_unqreads_F);

  return numreads;
}


int ChrData_chrlen(char* c, char** chrnames, int* chrlens, int numchroms) 
{  
	int i;
	for (i=0; i<numchroms; i++) {
		if (strcmp(chrnames[i], c) == 0) 
			return chrlens[i];
	}
	return -1;   
}

void readChrLen(char* file, int** chrlens)
{
	char* buff = 0;
	int m;
	FILE* fp;
	char** p;
	int numlines = 0;
	int i;
	int mynmax = 100000;
	
	buff  = (char*)malloc(mynmax * sizeof(char));
	fp    = fopen(file, "r");
	if (!fp) {
		printf("Please enter a valid filename (%s invalid)\n", file);
		exit(0);
	}
	
	numlines	= nbLinesInFile(file); 
	*chrlens	= (int*)calloc(numlines, sizeof(int));

	// read lines
	i = 0;
	while (fgets(buff, mynmax, fp) != 0) {
		chomp(buff);
		split_line_delim(buff, "\t", &p, &m);
		(*chrlens) [i] = atoi(p[1]);
		free(p);
		i++;
	}
	free(buff);
	fclose(fp);
}

void readChrData(char* file, char*** chrnames, int** chrlens, float** chrmap, int* numchroms)
{
	char* buff = 0;
	int m;
	FILE* fp;
	char** p;
	int numlines = 0;
	int i;
	int mynmax = 100000;
	
	buff  = (char*)malloc(mynmax * sizeof(char));
	fp    = fopen(file, "r");
	if (!fp) {
		printf("Please enter a valid filename (%s invalid)\n", file);
		exit(0);
	}
	
	numlines	= nbLinesInFile(file); 
	*chrnames	= (char**)calloc(numlines, sizeof(char*));
	*chrlens	= (int*)calloc(numlines, sizeof(int));
	*chrmap		= (float*)calloc(numlines, sizeof(float));
	for (i=0; i<numlines; i++) {
		(*chrmap)[i] = -1.0;
	}
	
	*numchroms = numlines;
	// read lines
	i = 0;
	while (fgets(buff, mynmax, fp) != 0) {
		chomp(buff);
		split_line_delim(buff, "\t", &p, &m);
		(*chrnames)[i] = strdup(p[0]);
		(*chrlens) [i] = atoi(p[1]);
		if ((m == 3) && (chrmap != 0)) {
			(*chrmap) [i] = atof(p[2]);
		}
		free(p);
		i++;
	}
	free(buff);
	fclose(fp);
}

void readChrData2(char* file1, char* file2, int hasid1, int hasid2, char*** chrnames, int* numchroms)
{
	char* buff1	= 0;
	char* buff2	= 0;
	int m;
	FILE* fp1;
	FILE* fp2;
	char** p1;
	char** p2;
	//int numlines = 0;
	int i;
	int j;
	int mynmax	= 100000;
	char *tmp	= NULL;
	
	//printf("hasid1: %d\t hasid2:%d\n", hasid1, hasid2);
	
	//open file1
	buff1  = (char*)malloc(mynmax * sizeof(char));
	fp1    = fopen(file1, "r");
	if (!fp1) {
		printf("Please enter a valid filename (%s invalid)\n", file1);
		exit(0);
	}
	
	//init set
	set_t s1;
	set_init(&s1); 
	
	//read lines of file1
	while (fgets(buff1, mynmax, fp1) != 0) {
		chomp(buff1);
		split_line_delim(buff1, "\t", &p1, &m);
		
		//add in set
		if (hasid1 == 1) {
			/* string copy */
			tmp = strdup(p1[1]);
			if (set_add_item(&s1, tmp) == 1)
				free(tmp);
		}
		else {
			/* string copy */
			tmp = strdup(p1[0]);
			if (set_add_item(&s1, tmp) == 1)
				free(tmp);
		}
		
		free(p1);
	}
	free(buff1);
	
	//open file2
	buff2  = (char*)malloc(mynmax * sizeof(char));
	fp2    = fopen(file2, "r");
	if (!fp2) {
		printf("Please enter a valid filename (%s invalid)\n", file2);
		exit(0);
	}
	
	//read lines of file2
	j=0;
	while (fgets(buff2, mynmax, fp2) != 0) {
		chomp(buff2);
		split_line_delim(buff2, "\t", &p2, &m);
		
		//add in set
		if (hasid2 == 1) {
			/* string copy */
			tmp = strdup(p2[1]);
			if (set_add_item(&s1, tmp) == 1)
				free(tmp);
		}
		else {
			/* string copy */
			tmp =  strdup(p2[0]);
			if (set_add_item(&s1, tmp) == 1)
				free(tmp);
		}
		
		free(p2);
		j++;
	}
	
	free(buff2);
	
	*numchroms	= s1.size;
	*chrnames	= (char**)calloc(s1.size, sizeof(char*));
	
	//fill in the chrnames table
	for (i = 0; i < s1.size; i++)
		(*chrnames)[i] = s1.arr[i];
	
	//delete set
	set_del(&s1);

	/* cleanup */
	(void)fclose(fp1);
	(void)fclose(fp2);
}


int SAMbitflagNum(int bitflag, int num)
{
	return (bitflag & (1 << num));
}

int SAMbitflag_ismapped(int bitflag) 
{
	return !(bitflag & ( 1 << 2 )); 
}


int SAMbitflag_strand(int bitflag)
{
	if ((bitflag & (1 << 4)) ||
		(bitflag & (1 << 5))) {
		return -1;
	} else {
		return 1;
	}
}

int SAMbitflag_isPairedProper(int bitflag)
{
	return SAMbitflagNum(bitflag, 1);
}


void geneateManyRandomSequences(float** wm, int w, int m, char*** seqs)
{
	int i=0;
	//char* s;
	(*seqs) = (char**)malloc(m * sizeof(char*));
	
	for (i=0; i<m; i++) {
		generateRandomSequence(wm, w, &((*seqs)[i]));
		//printf("%s\n", (*seqs)[i]);
	}
	
}

//
// sample a sequence from the motif distribution
//
void generateRandomSequence(float** wm, int w, char** s) 
{
	
	int i, j = 0, k;
	char* nt = "ACTG";
	float a_cum[5];
	float d = 0;
	char c = '-';
	
	*s = (char*)calloc(w+1, sizeof(char));
	
	for (j=0; j<w; j++) {
		
		i = 1;
		a_cum[0] = 0.0;
		for (k=0; k<4; k++) {
			a_cum[i] = a_cum[i-1] + wm[j][k];
			
			i++;
		}
		a_cum[i-1] = 1.0;
		
		d = default_rand();
		/*
		 printf("\n");
		 for (i=0; i<5; i++)
		 printf("%f\n", a_cum[i]);
		 printf("d=%f\n", d);
		 */
		
		for (i=1; i<=4; i++) {
			if ((d > a_cum[i-1]) && (d <= a_cum[i])) {
				c = nt[i-1];
				break;
			}
		}
		
		
		(*s)[j] = c;
	}
    
	(*s)[w] = '\0';
}




void sw(char* ss1, char* ss2, int R, int P, int O, int G, float* id, int* alnlen, int* qst, int* qen, int* sst, int* sen, int* sco, char** aln1, char** aln2) 
{
	
	int**  A;  // scores
	int*** L; // links 
	char*  s1;
	char*  s2;
	int    newsc;
	int    n1 = strlen(ss1);
	int    n2 = strlen(ss2);
	int    j_idx, i_idx;
	int    i, j, k;
	int    myi, myj;
	
	int    total_max_j = 0;
	int    total_max_i = 0;
	int    total_maxsc = 0;
	int    verbose = 0;
	int    maxsc;
	int    sc;
	int    maxil;
	int    maxjl;
	char*  as1;
	char*  as2;
	int    ii, jj, li, lj;
	
	s1 = (char*)calloc(n1+2, sizeof(char));
	s2 = (char*)calloc(n2+2, sizeof(char));
	
	s1[0] = 'N';
	memcpy(s1+1, ss1, n1);
	s1[n1+1] = '\0';
	
	s2[0] = 'N';
	memcpy(s2+1, ss2, n2);
	s2[n2+1] = '\0';
	
	n1++;
	n2++;
	
	// create score array
	A = (int**)calloc(n1, sizeof(int*));  
	for (i=0; i<n1; i++) {    
		A[i] = (int*)calloc(n2, sizeof(int));
	}
	
	for (i=0; i<n1; i++) {
		A[i][0] = 0;
	}      
	
	for (j=0; j<n2; j++) {
		A[0][j] = 0;
	}      
	
	// create link array
	L = (int***)calloc(n1, sizeof(int**));  
	for (i=0; i<n1; i++) {    
		L[i] = (int**)calloc(n2, sizeof(int*));
		for (j=0; j<n2; j++)
			L[i][j] = (int*)calloc(2, sizeof(int));
	}
	
	L[0][0][0] = -1;
	L[0][0][0] = -1;
	
	for (i=0; i<n1; i++) {    
		for (j=0; j<n2; j++) {
			
			if ((i == 0) && (j == 0))
				continue;
			
			if (s1[i] == s2[j])
				sc = R;
			else 
				sc = G;
			
			maxsc =  0;
			maxil = -1;
			maxjl = -1;
			
			// diag score
			if ((i>0) && (j>0)) {
				newsc = A[i-1][j-1] + sc;
				if (newsc > maxsc) {
					maxsc = newsc;
					maxil = i-1;
					maxjl = j-1;
				}
			} // if
			
			// up
			for (ii = i-1; ii>=0; ii--) {
				newsc = A[ii][j] + (i-ii) * G + O;
				if (newsc > maxsc) {
					maxsc = newsc;
					maxil = ii;
					maxjl = j;
				}
			}
			
			// left
			for (jj = j-1; jj>=0; jj--) {
				newsc = A[i][jj] + (j-jj) * G + O;
				if (newsc > maxsc) {
					maxsc = newsc;
					maxil = i;
					maxjl = jj;
				}
			}
			
			// store links
			L[i][j][0] = maxil;
			L[i][j][1] = maxjl;
			
			// store max score
			A[i][j]    = maxsc;
			
			// keep track of max score
			if (maxsc > total_maxsc) {
				total_maxsc = maxsc;
				total_max_i = i;
				total_max_j = j;
			}
			
		}
	}
	
	
	
	if (verbose == 1) {
		printf("    ");
		for (j=0; j<n2; j++) {
			printf("%4d", j);
		}
		printf("\n");
		
		for (i=0; i<n1; i++) { 
			printf("%4d", i);
			for (j=0; j<n2; j++) {
				printf("%4d", A[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
	
	// done with filling array
	i = total_max_i;
	j = total_max_j;
	
	*qen =  i-1;  // because added N to sequence ...
	*sen =  j-1;
	*qst = -1;
	*sst = -1;
	
	i_idx = 0;
	j_idx = 0;
	
	// string
	int maxal = max(n1,n2)*2;
	as1 = (char*)calloc(maxal, sizeof(char));
	as2 = (char*)calloc(maxal, sizeof(char));
	
	while (A[i][j] > 0) {
		
		*qst = i;
		*sst = j;
        
		myi = L[i][j][0];
		myj = L[i][j][1];
		
		//printf("Go to myi=%d, myj=%d\n", myi, myj);
		
		//if ((myi == -1) && (myj == -1))
		//  break;
		
		li = i - myi;
		lj = j - myj;
		
		
		if (li == 0) {
			for (k=j; k>myj; k--) {
				as1[i_idx++] = '-';
				as2[j_idx++] = s2[k];
			}
		} else if (lj == 0) {
			for (k=i; k>myi; k--) {
				as1[i_idx++] = s1[k];
				as2[j_idx++] = '-';
			}
		} else {
			as1[i_idx++] = s1[i];
			as2[j_idx++] = s2[j];
		}
		
		i = myi;
		j = myj;
		
		
	}
	
	*qst = *qst - 1;
	*sst = *sst - 1;
	
	// reverse strings  
	char* rev_as1 = (char*)calloc(i_idx+1, 1);
	char* rev_as2 = (char*)calloc(j_idx+1, 1);
	
	for (i=i_idx-1,j=0; i>=0; i--,j++)
		rev_as1[j] = as1[i];
	
	for (i=j_idx-1,j=0; i>=0; i--,j++)
		rev_as2[j] = as2[i];
	
	rev_as1[i_idx] = '\0';
	rev_as2[j_idx] = '\0';
	
	//printf("rev_sa1=%s\n", rev_as1);
	//printf("rev_sa2=%s\n", rev_as2);
	
	free(s1);
	free(s2);
	free(as1);
	free(as2);
	
	
	for (i=0; i<n1; i++)
		free(A[i]);
	free(A);
	
	
	for (i=0; i<n1; i++) {    
		for (j=0; j<n2; j++)
			free(L[i][j]);
		free(L[i]);
	}
	
	free(L);
	
	// sim
	for (i=0; i<i_idx; i++) {    
		if (rev_as1[i] == rev_as2[i]) 
			(*id) += 1;
	}
	*id     = *id / i_idx;
	
	*aln1   = rev_as1;
	*aln2   = rev_as2;
	*alnlen = i_idx;
	*sco    = total_maxsc;
	
	
}



int rand_chr_interval_startpoint_around_interval(char* c, int st, int en, int d)
{
	
	int strand;
	int enrand; 
	int drand;
	int len = en - st + 1;  // region length
	
	int chrlen = hg18_chrlen(c);
	if (chrlen == -1) {
		printf("Not recognizing %s\n", c);
		exit(0);
	}
	
	// pick region 
	if (default_rand() < 0.5) {
		strand = max(0, st-d - len);  // extend because we will look for start points
		enrand = st-len;              // so that the random peak does not overlap with actual peak
	} else {
		strand = en;
		enrand = min(en+d, chrlen - 1);
	}
    
	drand  = enrand - strand + 1;
	
	int p    = strand + (int)(0.5 + default_rand() * drand);
	
	return p;
}


int rand_chr_interval(char* c, int l)
{
	int chrl = hg18_chrlen(c);
	//printf("chrl=%d\n", chrl);
	if (chrl == -1) {
		printf("Not recognizing %s\n", c);
		exit(0);
	}
	
	int i    = (int)(0.5 + default_rand() * chrl - l);
	return i;
}

int rand_chr_interval2(char* c, int l, int m)
{
	//int chrl = hg18_chrlen(c);
	//printf("chrl=%d\n", chrl);
	if (m == -1) {
		printf("Not recognizing %s\n", c);
		exit(0);
	}  
	
	int i    = (int)(0.5 + default_rand() * m - l);
	return i;
}



/*
 CS
 FUNCTION: tally up read counts at all genomic positios in a vector
 file = read file
 format = 0=eland, 1=mit, 2=bed, 3=sam 
 char = useless
 chrlen = chrom len
 fraglen = actual fragment length (needed for extension); must be > 0 for extension to occur
 counts = counts (output)
 numreads = num reads (output) 
 */ 
void getCountFromReadFile(char* file, int format, char* chr, long chrlen, int readlen, int fraglen, int uniquereads,  int uniqmap, unsigned short** counts, int* numreads, int* numclonalreads)
{
	
	char** p		= 0;
	int   i		 	= 0;
	int   pos		= 0;
	int   lenread	= 0;
	int            st			  = 0;
	unsigned char* a_unqreads_R   = 0;
	unsigned char* a_unqreads_F   = 0;
	MappedRead     r;
	ReadI          ri;
	int            map;
	
	if (uniquereads == 1) {
		a_unqreads_R = create_binarized_array(chrlen);
		a_unqreads_F = create_binarized_array(chrlen);
	}
	
	*numreads		= 0;
	*numclonalreads       = 0;
	
	if (!ReadIopen(&ri, file, format)) {
		printf("Cannot open file %s\n", file);
	}
	
	while (nextRead(&ri, &r)) {
		
		//printf("%s\t%s\t%d\t%d\t%d\n", r.readname, r.seqname, r.pos, r.st, r.uniqmap);
		
		st      = r.st;
		map     = r.uniqmap;
		pos     = r.pos;
		lenread = r.lenread;
		
		if (map == 0) // next if read not uniquely mappable
			continue;
		
		// eliminate clonal read if required
		if (uniquereads == 1) {
			if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
				((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
				(*numclonalreads)++;
				continue;	      
			} else {
				if (st == 1)
					set_entry_in_binarized_array(a_unqreads_F, pos);
				else
					set_entry_in_binarized_array(a_unqreads_R, pos);
			}
		} // end eliminate
		
		
		// finally increase counts
		
		// count start and end
		int cnt_st = 0;
		int cnt_en = 0;
		
		if (st == 1) {
			
			cnt_st = pos; // count will start at beginning of read
			
			if (fraglen > 0) { // if extension required				      
				cnt_en = pos + fraglen;
			} else {	       // no extension required, just use read length
				cnt_en = pos + lenread;
			}
			
		} else if (st == -1) {
			
			cnt_en = pos+lenread; // count will start at beginning of read
			
			if (fraglen > 0) { // if extension required				      
				cnt_st = pos - (fraglen - lenread);
			} else {	       // no extension required, just use read length
				cnt_st = pos ;
			}
			
		} else {
			printf("Problem interpreting strand (%s)\n", p[8]);
			exit(0);
		}
		
		// truncate in case we reached the end of chr
		cnt_en = min(cnt_en, chrlen);
		cnt_st = max(0, cnt_st);
		
		// update counts
		for (i=cnt_st; i<cnt_en; i++) {
			if ((*counts)[i] < 65535) 
				(*counts)[i]++;					
		}
		
		(*numreads) ++;
	} // while nextRead
	
	//LineIclose(&li);
	ReadIclose(&ri);
	
	if (uniquereads == 1) {
		free(a_unqreads_F);
		free(a_unqreads_R);
	}
	
}

void getCountFromReadFileOLD(char* file, int format, char* chr, long chrlen, int readlen, int fraglen, unsigned short** counts, int* numreads)
{
	
	FILE* fp		= 0;
	int   bufflen	= 10000000;
	char* buff;
	int   cnt;
	//char* ptr;
	char* s;
	char** p;
	int   m;
	char* prevs = 0;
	int   i			= 0;
	char* mys;
	int   j			= 0;
	//int   cntc		= 0;
	int   pos		= 0;
	int   lenread	= 0;
	int   st		= 0;
	//int   mindist	= 0;
	char* tmpbuff   = 0;
	
	fp = fopen(file, "r");
	if (fp == 0) {
		printf("Cannot open read file %s.\n", file);
		exit(1);
	}
	
	*numreads		= 0;
	
	j = 0;
	while (!feof(fp)) {
		
		// allocate buffer for new chunk
		buff   = (char*)calloc(bufflen, sizeof(char)); 
		if (!buff) {
			die("Cannot allocate buff in seqI_open_fast\n");
		}    
		tmpbuff = buff; // added by OE, Dec 21, 2009
		
		cnt = fread(buff, 1, bufflen, fp);    
		
		i = 0;
		
		s = strsep(&buff, "\n");
		if (j != 0) { 
			mys = (char*)calloc(1000, sizeof(char));
			//printf("Adding %s\n", prevs);
			strcat(mys, prevs);
			//printf("Adding %s\n", s);
			if (s != 0)
				strcat(mys, s);
			s = mys;
			free(prevs);
		}
		
		while (1) {      
			
			//split_line_delim(s, "\t", &p, &m);	
			
			prevs = s;
			i++;
			
			s = strsep(&buff, "\n");
			if (s == 0) {
				prevs = strdup(prevs);
				break;
			} else {
				//printf("%s\n", prevs);
				
				
				split_line_delim(prevs, "\t", &p, &m);
				
				// ELAND
				int process = 0;
				if (format == 0) { 
					pos = atoi(p[7]);      
					if (p[8][0] == 'R')
						st = -1;
					else 
						st = 1;
					if (readlen <= 0)
						lenread = strlen(p[1]);
					else
						lenread = readlen;
					process = 1;
					
					// MIT
				} else if (format == 1) {
					pos     = atoi(p[1]);   
					lenread = atoi(p[2]) - pos + 1;
					if (p[3][0] == '-')
						st = -1;
					else 
						st = 1;
					//printf("pos=%d, leread=%d, st=%d\n", pos, lenread, st);
					process = 1;
					
					// BED
				} else if (format == 2) { // && (strcmp(p[3], "U0") == 0)) {
					pos     = atoi(p[1]);
					lenread = atoi(p[2]) - pos + 1;
					if (p[5][0] == '-')
						st = -1;
					else 
						st = 1;
					process = 1;
				}
				
				if (process == 1) {
					// increase counts
					
					// count start and end
					int cnt_st = 0;
					int cnt_en = 0;
					
					if (st == 1) {
						
						cnt_st = pos; // count will start at beginning of read
						
						if (fraglen > 0) { // if extension required				      
							cnt_en = pos + fraglen;
						} else {	       // no extension required, just use read length
							cnt_en = pos + lenread;
						}
						
					} else if (st == -1) {
						
						cnt_en = pos+lenread; // count will start at beginning of read
						
						if (fraglen > 0) { // if extension required				      
							cnt_st = pos - (fraglen - lenread);
						} else {	       // no extension required, just use read length
							cnt_st = pos ;
						}
						
					} else {
						printf("Problem interpreting strand (%s)\n", p[8]);
						exit(0);
					}
					
					// truncate in case we reached the end of chr
					cnt_en = min(cnt_en, chrlen);
					cnt_st = max(0, cnt_st);
					
					// update counts
					for (i=cnt_st; i<cnt_en; i++) {
						if ((*counts)[i] < 65535) 
							(*counts)[i]++;					
					}
					
					(*numreads) ++;
				}
				
				free(p); // p[x] are just pointers
			}
		}
		
		//printf("s=%s\n", prevs);
		
		//getchar();
		
		free(tmpbuff);
		j++;
	}
	
}

int chr2index(char* chr) 
{
	int l = strlen(chr);
	
	if (l > 5) 
		return -1;  
	else if (chr[3] == 'M')
		return 0;
	else if (chr[3] == 'X')
		return 23;
	else if (chr[3] == 'Y')
		return 24;
	else {
		return atoi(chr+3);    
	}
}


float hg18_30mer_mappability(char* c) 
{
	int i;
	char* chrnames[] = 
    {"chr18",
		"chr2",
		"chr3",
		"chr4",
		"chr6",
		"chr8",
		"chr5",
		"chr12",
		"chr11",
		"chr20",
		"chr7",
		"chr10",
		"chr17",
		"chrX",
		"chr1",
		"chr13",
		"chr16",
		"chr14",
		"chrM",
		"chr9",
		"chr19",
		"chr15",
		"chr21",
		"chr22",
		"chrY"};
	
	float chrmap[] = {
		0.123087459668913,
		0.139325292921335,
		0.139622164963933,
		0.141630052737745,
		0.142365618132972,
		0.143106107677065,
		0.143924633059643,
		0.14645199279659,
		0.149044906485258,
		0.160449000194824,
		0.164855517225434,
		0.169669404417753,
		0.184021980040252,
		0.201849540099583,
		0.223202247602959,
		0.251159230291692,
		0.265147676410215,
		0.279122120502026,
		0.279283084907368,
		0.307441635415995,
		0.317359834491667,
		0.317702957023205,
		0.359312392256674,
		0.426161556382597,
		0.747921665906161,
	};
	
	for (i=0; i<25; i++) {
		if (strcmp(chrnames[i], c) == 0) {
			return 1.0 - chrmap[i];
		}
	}	
	return -1;
}

//test values
float hg19_30mer_mappability(char* c) 
{
	int i;
	char* chrnames[] = 
    {"chr18",
		"chr2",
		"chr3",
		"chr4",
		"chr6",
		"chr8",
		"chr5",
		"chr12",
		"chr11",
		"chr20",
		"chr7",
		"chr10",
		"chr17",
		"chrX",
		"chr1",
		"chr13",
		"chr16",
		"chr14",
		"chrM",
		"chr9",
		"chr19",
		"chr15",
		"chr21",
		"chr22",
		"chrY"};
	
	float chrmap[] = {
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
		0.8,
	};
	
	for (i=0; i<25; i++) {
		if (strcmp(chrnames[i], c) == 0) {
			return 1.0 - chrmap[i];
		}
	}	
	return -1;
}



int hg18_chrlen(char* c) 
{
	
	int i;
	char* chrnames[] = 
    {"chr1",
		"chr2",
		"chr3",
		"chr4",
		"chr5",
		"chr6",
		"chr7",
		"chrX",
		"chr8",
		"chr9",
		"chr10",
		"chr11",
		"chr12",
		"chr13",
		"chr14",
		"chr15",
		"chr16",
		"chr17",
		"chr18",
		"chr19",
		"chr20",
		"chrY",
		"chr22",
		"chr21",
		"chrM"};
	
	int chrlens[] = {247249719,
		242951149,
		199501827,
		191273063,
		180857866,
		170899992,
		158821424,
		154913754,
		146274826,
		140273252,
		135374737,
		134452384,
		132349534,
		114142980,
		106368585,
		100338915,
		88827254,
		78774742,
		76117153,
		63811651,
		62435964,
		57772954,
		49691432,
		46944323,
		16571};
	
	for (i=0; i<25; i++) {
		if (strcmp(chrnames[i], c) == 0) {
			return chrlens[i];
		}
	}	
	return -1;
}

void get_hg18_chroms(char*** chroms, int* numchroms, int** chrlens)
{
	int i;
	*numchroms   = 24;
	*chroms      = (char**)calloc(*numchroms, sizeof(char*));
	(*chroms)[0] =  strdup("chrY");		//array that holds the chromosomes names
	(*chroms)[1] =  strdup("chrX");
	(*chroms)[2] =  strdup("chr9");
	(*chroms)[3] =  strdup("chr8");
	(*chroms)[4] =  strdup("chr7");
	(*chroms)[5] =  strdup("chr6");
	(*chroms)[6] =  strdup("chr5");
	(*chroms)[7] =  strdup("chr4");
	(*chroms)[8] =  strdup("chr3");
	(*chroms)[9] =  strdup("chr22");
	(*chroms)[10] = strdup("chr21");
	(*chroms)[11] = strdup("chr20");
	(*chroms)[12] = strdup("chr2");
	(*chroms)[13] = strdup("chr19");
	(*chroms)[14] = strdup("chr18");
	(*chroms)[15] = strdup("chr17");
	(*chroms)[16] = strdup("chr16");
	(*chroms)[17] = strdup("chr15");
	(*chroms)[18] = strdup("chr14");
	(*chroms)[19] = strdup("chr13");
	(*chroms)[20] = strdup("chr12");
	(*chroms)[21] = strdup("chr11");
	(*chroms)[22] = strdup("chr10");
	(*chroms)[23] = strdup("chr1");
	
	*chrlens      = (int*)calloc(*numchroms, sizeof(int));
	for (i=0; i<*numchroms; i++) {
		(*chrlens)[i] = hg18_chrlen( (*chroms)[i] );
	}
	
}


int hg19_chrlen(char* c) 
{
	
	int i;
	char* chrnames[] = 
    {"chr1",
		"chr2",
		"chr3",
		"chr4",
		"chr5",
		"chr6",
		"chr7",
		"chr8",
		"chr9",
		"chr10",
		"chr11",
		"chr12",
		"chr13",
		"chr14",
		"chr15",
		"chr16",
		"chr17",
		"chr18",
		"chr19",
		"chr20",
		"chr21",
		"chr22",
		"chrM",
		"chrX",
		"chrY"};
	
	int chrlens[] = {249250621,
		243199373,
		198022430,
		191154276,
		180915260,
		171115067,
		159138663,
		146364022,
		141213431,
		135534747,
		135006516,
		133851895,
		115169878,
		107349540,
		102531392,
		90354753,
		81195210,
		78077248,
		59128983,
		63025520,
		48129895,
		51304566,		
		16571,
		155270560,
		59373566};
	
	for (i=0; i<25; i++) {
		if (strcmp(chrnames[i], c) == 0) {
			return chrlens[i];
		}
	}	
	return -1;
}


long sequencesOverlap(long s1, long e1, long s2, long e2) {
	
	long p1 = max(s1, s2);    
	long p2 = min(e1, e2);
	
	long ll = p2 - p1 + 1;
	
	if ( ll > 0  ) {
		return ll;
	} else {
		return 0;
	}
    
}

long sequencesDistance(long s1, long e1, long s2, long e2) {
	
	long p1 = max(s1, s2);    
	long p2 = min(e1, e2);
	
	long ll = p1 - p2;
	long l;
	
	if ( ll > 0  ) {
		//long l = e2 - s1;
		if ( e1 > e2 ) {
			l = e2 - s1;
		}
		else {
			l = s2 - e1;
		}
		return l;
	} 
	else {
		return 0;
	}
    
}

//  P(A1,A2,A3) = P(A1)P(A2)P(A3)  (all col indep)
//  P(A1,A2,A3) = P(A1,A2)P(A3|A1,A2)   = P(A1)P(A2,A3|A1) = P(A1)
//                        P(A3|A1,A2) = 
//              = P(A1)P(A2|A1)P(A3|A2)

//  comes down to P(A3|A1,A2) = P(A3|A2) if A1 and A3 are indep, ie P(A1,A3)=0 


// input a file, returns an array
void readScanACEMatches(char* file, ScanACE_match** scm, int* nmp) 
{
	
	char* buff;
	int   mynmax = 1000;
	FILE* fp;
	char* genename;
	int   position;
	int   strand;  
	float score;   
	int   maxscanacematches = 100000;
	int   nm = 0;
	//
	// read scanace matches
	//  
	*scm = (ScanACE_match*)malloc(maxscanacematches * sizeof(ScanACE_match));
    
	buff  = (char*)malloc(mynmax * sizeof(char));
	fp    = fopen(file, "r");
	if (!fp) {
		printf("Please enter a valid filename (%s invalid)\n", file);
		exit(0);
	}
	
	
	// read matches
	while (!feof(fp)) {
		
		fgets(buff, mynmax, fp);
		
		if (feof(fp))
			break; 
		
		chomp(buff);
		
		//
		// get the row name 
		//
		genename =    mystrtok(buff, '\t');
		position = atoi(mystrtok(0, '\t'));
		strand   = atoi(mystrtok(0, '\t'));
		score    = atof(mystrtok(0, '\t'));
		
		(*scm)[nm].gene     = genename;
		(*scm)[nm].strand   = strand;
		(*scm)[nm].position = position;
		(*scm)[nm].score    = score;
		
		nm ++;
		
		if (nm == maxscanacematches) {
			die("maxscanacematches reached, please modify source and recompile.\n");
		}
		
	}
	
	
	*nmp = nm;
	
	fclose(fp);
	
	
}


void initialize_nt() 
{
	int i, j;
	
	N = (int*)calloc(256, sizeof(int));
	
	N['A']  = 0;
	N['C']  = 1;
	N['T']  = 2;
	N['G']  = 3;
	N['N']  = 4;
	N['X']  = 5;
	
	C = (int*)calloc(256, sizeof(int));
	
	C['A']  = 2;
	C['C']  = 3;
	C['T']          = 0;
	C['G']          = 1;
	C['N']          = 4;
	C['X']          = 4;
	
	ntoC[0]         = 'A';
	ntoC[1]         = 'C';
	ntoC[2]         = 'T';
	ntoC[3]         = 'G';
	ntoC[4]         = 'N';
	ntoC[5]         = 'X';
	
	char_compl      = (char*)calloc(256, sizeof(char));
	char_compl['A'] = 'T';
	char_compl['T'] = 'A';
	char_compl['C'] = 'G';
	char_compl['G'] = 'C';
	
	N1 = (int**)calloc(256, sizeof(int*));
	for (i=0; i<256; i++) {  
		N1[i] = (int*)calloc(256, sizeof(int));
	}
	
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			N1[ ntoC[i] ][ ntoC[j] ] = i*4+j;
			ntoC1[ i*4+j ][0]          = ntoC[i];
			ntoC1[ i*4+j ][1]          = ntoC[j]; 
			ntoC1[ i*4+j ][2]          = '\0';
		}
	}
	
	
}


void printWM(float** wm, int w)
{
	int i,j;
	for (j=0; j<4; j++) {
		
		for (i=0; i<w; i++) {
			printf("\t%3.2f", wm[i][j]);
		}
		printf("\n");
	}
}

void printBckd(float* m)
{
	int i;
	printf("BKG");
	for (i=0; i<4; i++) {
		printf("\t%4.3f", m[i]);
	}
	printf("\n");
	
}

//
// read PBM matrix 
//
void read_BULYK_WM(char* wmfile, float*** wm, int* w)
{
	char* buff;
	int   len = 100000;
	FILE* fp = 0;
	char** p;
	int    m, i, col;
	
	buff = (char*)malloc(1000*sizeof(char));
	
	fp = fopen(wmfile, "r");
	
	// order in matrices is ACGT
	// order in outout is ACTG
	
	// rid of row 1
	fgets(buff, len, fp);
	
	//
	// Motif line 1
	//
	fgets(buff, len, fp);
	chomp(buff);
	split_line_delim(buff, "\t", &p, &m);
	
	// must count nt
	*w = m-1;
	
	//
	// allocate memory for matrix (attention matrix is inverted)
	//
	*wm = (float**)malloc((m-1) * sizeof(float*));
	for (col= 0; col < m-1; col++) {
		(*wm)[col] = (float*)malloc(4 * sizeof(float));
		for (i=0; i<4; i++) {
			(*wm)[col][i] = 0;
		}
	}
	
	
	// enter data from line 1 (A)
	for (i=1; i<m; i++) {
		//printf("'%s' ", p[i]);
		(*wm)[i-1][0] = atof(p[i]);
	}
	//printf("\n");
	free(p);
	
	
	// read data from line 2 (C)
	fgets(buff, len, fp);
	chomp(buff);
	split_line_delim(buff, "\t", &p, &m); 
	for (i=1; i<m; i++) {
		(*wm)[i-1][1] = atof(p[i]);
	}
	free(p);
	
	// read data from line 3 (G, is line 3 in 0-3)
	fgets(buff, len, fp);
	chomp(buff);
	split_line_delim(buff, "\t", &p, &m);   
	for (i=1; i<m; i++) {
		(*wm)[i-1][3] = atof(p[i]);
	}
	free(p);
	
	// read data from line 3 (G, is line 3 in 0-3)
	fgets(buff, len, fp);
	chomp(buff);
	split_line_delim(buff, "\t", &p, &m);   
	for (i=1; i<m; i++) {
		(*wm)[i-1][2] = atof(p[i]);
	}
	free(p);
	
	
}



void read_JASPAR_WM(char* wmfile, int*** wm, int* w, int* nsites)
{
	char* buff;
	int   len = 100000;
	FILE* fp = 0;
	char** p;
	int    m, i, j, col;
	
	buff = (char*)malloc(1000*sizeof(char));
	
	fp = fopen(wmfile, "r");
	
	// order in matrices is ACGT
	// order in outout is ACTG
	
	//
	// Motif line 1
	//
	fgets(buff, len, fp);
	chomp(buff);
	split_line_delim(buff, " ", &p, &m);
	
	*w = m;
	
	//
	// allocate memory for matrix (attention matrix is inverted)
	//
	*wm = (int**)malloc(m * sizeof(int*));
	for (col= 0; col < m; col++) {
		(*wm)[col] = (int*)malloc(4 * sizeof(int));
		for (i=0; i<4; i++) {
			(*wm)[col][i] = 0;
		}
	}
	
	
	// enter data from line 1 (A)
	for (i=0; i<m; i++) {
		//printf("'%s' ", p[i]);
		(*wm)[i][0] = atoi(p[i]);
	}
	//printf("\n");
	free(p);
	
	
	// read data from line 2 (C)
	fgets(buff, len, fp);
	chomp(buff);
	split_line_delim(buff, " ", &p, &m); 
	
	for (i=0; i<m; i++) {
		//printf("%s ", p[i]);
		(*wm)[i][1] = atoi(p[i]);
	}
	//printf("\n");
	free(p);
	
	// read data from line 3 (G, is line 3 in 0-3)
	fgets(buff, len, fp);
	chomp(buff);
	split_line_delim(buff, " ", &p, &m); 
	
	for (i=0; i<m; i++) {
		//printf("%s ", p[i]);
		(*wm)[i][3] = atoi(p[i]);
	}
	//printf("\n");
	free(p);
	
	// read data from line 3 (G, is line 3 in 0-3)
	fgets(buff, len, fp);
	chomp(buff);
	split_line_delim(buff, " ", &p, &m); 
	
	for (i=0; i<m; i++) {
		//printf("%s ", p[i]);
		(*wm)[i][2] = atoi(p[i]);
	}
	//printf("\n");
	free(p);
	
	
	// calculate nsites
	int oldsum = 0;
	for (i=0; i<m; i++) {
		int sum = 0;
		for (j=0; j<4; j++) {
			sum += (*wm)[i][j];
		}
		if ((i > 0) && (sum != oldsum)) {
			printf("Problem: %d != %d\n", sum, oldsum);
			exit(0);
		}
		oldsum = sum;     
	}
	
	*nsites = oldsum;
	
}


//
//  read a single WM from a file (in the form of a count matrix)
//
//  transfo = 0 -> integer, 1 -> float, 2 -> log
//
void readACEintWM(char* wmfile, int*** wm, int* w, char*** sites, int* n, int** stars)
{
	
	FILE*   fp3;
	char*   buff;
	int     len = 100000;
	int**   mywm;
	int     cnt_sites;
	int     l;
	int     col;
	int     myw;
	int     i;
	int     star = 0;
	char**  mysites;
	int*    mystars;
	
	buff = (char*)malloc(1000*sizeof(char));
	
	//
	//  read in WM
	//
	fp3 = fopen(wmfile, "r");
	
	
	// Motif 1 line
	fgets(buff, len, fp3);
	
	// line 2, to get the motif width
	fgets(buff, len, fp3);
	myw = strlen(buff)-1; 
	cnt_sites = 1;
	
	// loop, get the number of sites
	while (!feof(fp3)) {
		//getline(&buff, &len, fp3);
		fgets(buff, len, fp3);  
		if (feof(fp3)) {
			break;
		}
		cnt_sites ++;
	}
	
	// start again
	rewind(fp3);
	
	
	//
	// input wm, log weights
	//
	mywm = (int**)malloc(myw * sizeof(int*));
	for (col= 0; col < myw; col++) {
		mywm[col] = (int*)malloc(5 * sizeof(int));
		for (i=0; i<5; i++) {
			mywm[col][i] = 0;
		}
	}
	
	mystars = (int*)malloc(myw * sizeof(int));
	
	//
	//  read in each site
	//
	mysites = (char**)malloc(cnt_sites * sizeof(char*));
	
	cnt_sites = 0;
	fgets(buff, len, fp3);
	
	while (!feof(fp3)) {
		fgets(buff, len, fp3);
		if (feof(fp3)) {
			break;
		}
		chomp(buff);
		l = strlen(buff);
		
		star = 0;
		for (i=0; i<l; i++) {
			if ((buff[i] == ' ') || (buff[i] == '*')) {
				star = 1;  break; 
			}
		}
		
		if (star == 1) {
			for (i=0; i<l; i++) {
				if (buff[i] == '*') 
					mystars[i] = 1;
				else
					mystars[i] = 0;
				
			}
			
			break;
		}
		
		for (i=0; i<l; i++) {
			
			mywm[i][ N[ (int)(buff[i]) ] ] += 1;
		}
		
		mysites[ cnt_sites ] = strdup(buff);
		
		cnt_sites++;
		
	}
	
	//
	//  SKIP: normalize by the number of sites
	//
	
	/*
	 for (col= 0; col < myw; col++) {
	 for (i=0; i<4; i++) {	
	 //+= bkg[i]; //epsilon;
	 int ini = mywm[col][i];
	 mywm[col][i]  = log2( ((mywm[col][i] + bkg[i]) / (float)(cnt_sites + 1)) / bkg[i] ); 
	 
	 }
	 }
	 */
	
	//
	//  output 
	//
	*w     = myw;
	*wm    = mywm;
	*sites = mysites;
	*n     = cnt_sites;
	*stars = mystars;
	
}


//
// 
//
void ACEintWMtologWM(int** intwm, int w, int nsites, float* rawbkg, float*** logwm) {
	
	int col, i;
	float** mywm;
	
	mywm = (float**)malloc(w * sizeof(float*));
	for (col= 0; col<w; col++) {
		mywm[col] = (float*)malloc(5 * sizeof(float));
		for (i=0; i<5; i++) {
			mywm[col][i] = 0.0;
		}
	}
	
	for (col= 0; col < w; col++) {
		for (i=0; i<4; i++) {	
			mywm[col][i]  = log2( ((intwm[col][i] + rawbkg[i]) / (float)(nsites + 1)) / rawbkg[i] ); 
		}
	}
	
	*logwm = mywm;
	
}



//
// 
//
void WMtologWMbkg(float** fwm, int w, float* rawbkg, float*** logwm) 
{
	
	int col, i;
	float** mywm;
	float eps = 0.000001;
	
	mywm = (float**)malloc(w * sizeof(float*));
	for (col= 0; col<w; col++) {
		mywm[col] = (float*)malloc(5 * sizeof(float));
		for (i=0; i<5; i++) {
			mywm[col][i] = 0.0;
		}
	}
	
	for (col= 0; col < w; col++) {
		for (i=0; i<4; i++) {	
			mywm[col][i]  = log2( 
								 (fwm[col][i] + eps) / rawbkg[i] 
								 ); 
		}
	}
	
	*logwm = mywm;
	
}





//
//  log score
//
float getScoreOnSeq(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn) 
{
	float score1 = 0.0;
	float score2 = 0.0;
	
	int  j, k;
	
	*hasn = 0;
	
	for (j=0; j<w; j++) {
		if (stars[j] == 1) {
			if ((seq[i + j] != 'A') && (seq[i + j] != 'T') && (seq[i + j] != 'C') && (seq[i + j] != 'G')) {
				*hasn = 1;
				return 0.0;
			}
			
			score1 += wml[j][ N[ (int)(seq[i + j]) ] ]; // - bkg[ N[ (int)(seq[i + j]) ] ];
		}
	}
	
	if (rna == 0) {  // if not RNA, also score reverse complement
		
		for (j=w-1,k=0; j>=0; j--,k++) {
			if (stars[k] == 1) {
				//printf(" += %f\n", wml[k][ C[ (int)(seq[i + j]) ] ]);
				score2 += wml[k][ C[ (int)(seq[i + j]) ] ]; // - bkg[ C[ (int)(seq[i + j]) ] ];
			}
		}
		
	}
	
	//printf("score1 0M = %f, score2 0M = %f\n", score1, score2);
	
	
	if ((rna == 1) || (score1 >= score2)) {
		*strand = 1;
		return score1;
	} else {
		*strand = -1;
		return score2;
	}
	
	
}


void findAllWeightMatrixMatches(float** wml, int w, int* stars, float* bkgl, float t, char* seq, int rna, int** matches_pos, int* num_matches, char** matches_ori, float** matches_sco, int max_num_matches)
{
	
	int i;
	int l = strlen(seq);
	float s;
	int strand;
	int my_num_matches = 0;
	int my_max_num_matches = max_num_matches;
	void* ptr;
	int  hasn;
	
	*matches_pos = (int*) malloc(my_max_num_matches * sizeof(int));  
	if (!*matches_pos) {
		die("findAllWeightMatrixMatches: not enough memory for matches_pos\n");    
	}
	*matches_ori = (char*)malloc(my_max_num_matches * sizeof(char));
	if (!*matches_ori) {
		die("findAllWeightMatrixMatches: not enough memory for matches_ori\n");    
	}
	*matches_sco = (float*)malloc(my_max_num_matches * sizeof(float));
	if (!*matches_sco) {
		die("findAllWeightMatrixMatches: not enough memory for matches_sco\n");    
	}
	
	for (i=0; i<l; i++) {
		
		s = getScoreOnSeq(wml, w, stars, bkgl, seq, i, rna, &strand, &hasn);
		
		if (!hasn && (s >= t)) {
			
			(*matches_pos)[ my_num_matches ] = i;
			(*matches_ori)[ my_num_matches ] = (char)strand;
			(*matches_sco)[ my_num_matches ] = s;
			
			my_num_matches ++;
			
			if (my_num_matches == my_max_num_matches) {
				my_max_num_matches += max_num_matches;
				ptr = realloc( *matches_pos, my_max_num_matches * sizeof(int));	
				if (ptr == 0)
					die("findAllWeightMatrixMatches: not enough memory for matches_pos (realloc)\n");
				ptr = realloc( *matches_ori, my_max_num_matches * sizeof(char));		
				if (ptr == 0)
					die("findAllWeightMatrixMatches: not enough memory for matches_ori (realloc)\n");
				ptr = realloc( *matches_sco, my_max_num_matches * sizeof(float));		
				if (ptr == 0)
					die("findAllWeightMatrixMatches: not enough memory for matches_sco (realloc)\n");
				
			}
			
		}    
	} 
	
	*num_matches = my_num_matches;
	
}


//
// this function scans all sequences and returns the max WM scores for each
//
void findAllSeqsMaxWMScores(char** seqs, int numseqs, float** logwm, int mt, int w, int* stars, float* bkg, int rna, float** scores, int** a_idx)
{
	int   i;
	float* myscores;
	int   idx;
	
	int*  my_a_idx = 0;
	
	myscores = (float*)malloc( numseqs* sizeof(float));
	if (myscores == 0) {
		die("findAllSeqsMaxWMScores: cannot allocate memory for myscores\n");
	}
	
	if (a_idx != 0)
		my_a_idx = (int*)malloc( numseqs * sizeof( int ));
	
	
	for (i=0; i<numseqs; i++) {    
		myscores[i] = findMaxWMScore(seqs[i], logwm, mt, w, stars, bkg, rna, &idx);
		if (a_idx != 0)
			my_a_idx[i] = idx;
		//printf("%f %s\n", myscores[i], substr(seqs[i], idx, w));
	}
	
	*scores  = myscores;
	
	if (a_idx != 0)
		*a_idx   = my_a_idx;
}


float findMaxWMScore(char* seq, float** logwm, int mt, int w, int* stars, float* bkg, int rna, int* idx)
{
	
	int i, l;
	float maxs = -1000000.0;
	float s = -100000;
	int   hasn, strand;
	
	l = strlen(seq);
	*idx = -1;
	
	//printWM(logwm, w);
	
	for (i=0; i<l-w; i++) {    
		
		if (mt == 0) 
			s = getScoreOnSeq(logwm, w, stars, bkg, seq, i, rna, &strand, &hasn);	    
		
		else if (mt == 1)
			s = getScoreOnSeq_1M(logwm, w, stars, bkg, seq, i, rna, &strand, &hasn);	 
		
		//printf("new max for %s at i=%d, score = %f\n", substr(seq, i, w), i, s);
		
		if ((hasn == 0) && (s > maxs)) {      
			
			//printf("new max for %s, score = %f\n", substr(seq, i, w), s);
			
			maxs = s;
			*idx = i;
		} 
		
	}
	
	return maxs;
}



//
//  get the score threshold out of a WM (log or not)
//
float getWMThreshold(float** wm, int w, int* stars, float* bkg, int transfo, char** sites, int n, float c) {
	
	float score; 
	int   i, j;
	//char* T    = "ACTG";
	float SX   = 0.0;
	float SX2  = 0.0;
	float AVG;
	float STD  = 0.0;
	float t;
	float* scores;
	int   cnt_above;
	
	scores = (float*)malloc(n * sizeof(float));
	
	for (i=0; i<n; i++) {
		
		//printf("site = %s\n", sites[i]);
		
		if (transfo == 0) {
			score = 1.0;
		} else {
			score = 0.0;
		}
		
		for (j=0; j<w; j++) {
			if (stars[j] == 0)
				continue;
			if (transfo == 0) {
				//printf(" %c mult pby %4.3f / %4.3f\n",sites[i][j], wm[j][ N[sites[i][j]] ],
				//     bkg[ N[ sites[i][j] ] ]);
				
				score *= wm[j][ N[  (int)(sites[i][j]) ] ] / bkg[ N[ (int)(sites[i][j]) ] ];
				//printf(" score = %4.3f\n", score);
			} else {
				score += wm[j][ N[ (int)(sites[i][j]) ] ]; // - bkg[ N[ (int)(sites[i][j]) ] ];
				//printf("j=%d %c mult pby %4.3f (%c) - %4.3f\n", j, sites[i][j], wm[j][ N[sites[i][j]] ], sites[i][j], bkg[ N[ sites[i][j]] ]);
				
			}
		}
		
		///printf("%s\t%4.3f\n", sites[i], score);
		
		SX  += score;
		SX2 += score * score;
		
		scores[i] = score;
		
	}
	
	AVG =  SX/n;
	STD =  (float)sqrt(( SX2 - SX*SX / n )  / (n - 1));
	//printf("AVG = %3.2f, S = %3.2f\n", AVG, STD);
	
	t = AVG - c * STD;
	
	cnt_above = 0;
	for (i=0; i<n; i++) {
		if (scores[i] >= t) {
			cnt_above ++;
		}
	}
	
	//printf("%d/%d sites score better than t\n", cnt_above, n);
	
	return t;
	
}




char* getGappedKmer(char* kmer, int gap) 
{
	
	char*    stmp;
	int      kmersize;
	int      j;
	int      half1, half2;
	
	
	kmersize = strlen(kmer);
	half1 = kmersize/2;
	half2 = kmersize-half1;
	
	stmp = (char*)calloc((kmersize+gap+1), sizeof(char));
	
	strncpy(stmp, kmer, half1);  
	for (j=0; j<gap; j++) {
		stmp[half1+j] = '.';
	}
	strncpy(stmp+half1+gap, kmer + half1, half2);
	stmp[kmersize+gap] = '\0';
	
	
	return stmp;
}


char* getGappedMotif(char* motif, int motifsize, int gap) 
{
	
	char*    newmotif;
	int      half1;
	int      l = strlen(motif);
	int      i = 0;  // cnt in ungapped motif
	int      j = 0;  // cnt in   gapped motif
	int      k = 0;  // cnt characters
	int      h = 0;
	int      p_b;
	
	newmotif = (char*)calloc((l+gap+1), sizeof(char)); 
	
	half1 = motifsize/2;
	
	while (i < l) {
		
		if (motif[i] == '[') {
			
			p_b = i;
			while (motif[p_b] != ']') {	
				p_b++; 
			} 
			
			while (i <= p_b) {
				newmotif[j] = motif[i];
				i++; j++;
			}
			
			k++;
			
		} else {
			
			newmotif[j] = motif[i];
			j++; i++; k++;
			
		}
		
		if (k == half1) {
			for (h=0; h<gap; h++) {
				newmotif[j] = '.';
				j++;
			}
		}
		
	}
	
	newmotif[l+gap] = '\0';
	
	return newmotif;
}



void getIntegerWMfromRegexp(char* motif, int gap, int*** wm, int* w)
{
	int   myw;
	int** mywm;
	
	int      l = strlen(motif);
	int      i = 0;  // cnt in ungapped motif
	int      j = 0;  // cnt in   gapped motif
	int      k = 0;  // cnt characters
	int      p_b;
	int      nl;
	int      h;
	
	myw = getRegexpMotifLength(motif, gap);
	
	// allocate memory
	mywm = (int**)calloc( myw, sizeof(int*) );
	for (i=0; i<myw; i++) {
		mywm[i] = (int*)calloc(4, sizeof(int));
	}
	
	// scan the regexp
	i = 0;
	
	while (i < l) {
		
		if (motif[i] == '[') {
			
			p_b = i;
			while (motif[p_b] != ']') {	
				p_b++; 
			} 
			
			nl = p_b - i - 1;
			
			if ((nl != 2) && (nl != 3)) {
				printf("Problem: nl should be 2 or 3, it is %d\n", nl);
				exit(0);
			}
			
			nl = 120 / nl;
			
			for (h=i+1; h<p_b; h++) {
				
				mywm[k][ N[ (int)(motif[h]) ] ] = nl;
				
			}
			
			while (i <= p_b) {	
				i++; j++;
			}
			
			k++;
			
		} else {
			
			if (motif[i] == '.') {
				for (h=0; h<4; h++) {
					mywm[k][h] = 120 / 4;
				}		
			} else if ((motif[i] == 'A') ||
					   (motif[i] == 'T') ||
					   (motif[i] == 'C') ||
					   (motif[i] == 'G')) {
				mywm[k][ N[ (int)(motif[i]) ] ] = 120;
			}
			
			j++; i++; k++;
		}    
	}  
	
	*w  = myw;
	*wm = mywm;
}

//
//  get a P(A)P(B|A)P(C|B) ... matrix
//
void getFirstOrderMarkovIntegerWMfromRegexp(char* motif, int gap, int*** wm, int* w)
{
	int      myw;
	int**    mywm;  
	int      l = strlen(motif);
	int      i = 0;  // cnt in ungapped motif
	int      j = 0;  // cnt in   gapped motif
	int      k = 0;  // cnt characters
	int      m = 0;
	int      p_b;
	int      nl;
	int      h;
	char**   a_c;
	int      nl0, nl1;
	
	myw = getRegexpMotifLength(motif, gap);
	
	// allocate memory
	mywm = (int**)calloc( myw, sizeof(int*) );
	for (i=0; i<myw; i++) {
		mywm[i] = (int*)calloc(16, sizeof(int));
	}
	
	// array of characters
	a_c = (char**)calloc(myw, sizeof(char*));
	
	// scan the regexp
	i = 0;
	
	while (i < l) {
		
		if (motif[i] == '[') {
			
			p_b = i;
			while (motif[p_b] != ']') {	
				p_b++; 
			} 
			
			nl = p_b - i - 1;
			
			if ((nl != 2) && (nl != 3)) {
				printf("Problem: nl should be 2 or 3, it is %d\n", nl);
				exit(0);
			}
			
			a_c[k] = (char*)calloc(nl+1, sizeof(char));
			
			
			// fill up with all allowed characters at that position
			for (h=i+1,m=0; h<p_b; h++,m++) {	
				a_c[k][m] = motif[h];	
			}
			a_c[k][nl] = '\0';
			
			while (i <= p_b) {	
				i++; j++;
			}
			
			k++;
			
		} else {
			
			if (motif[i] == '.') {
				
				// . characters
				a_c[k] = (char*)calloc(5, sizeof(char));	      
				for (h=0; h<4; h++) {
					a_c[k][h] = ntoC[h];
				}		
				a_c[k][4] = '\0';
				
				
			} else if ((motif[i] == 'A') ||
					   (motif[i] == 'T') ||
					   (motif[i] == 'C') ||
					   (motif[i] == 'G')) {
				
				a_c[k] = (char*)calloc(2, sizeof(char));
				a_c[k][0] = motif[i];
				a_c[k][1] = '\0';
				
			}
			
			j++; i++; k++;
		}    
	}  
	
	//for (i=0; i<myw; i++) {
	//  printf("a_c[%d] = %s\n", i, a_c[i]);
	//}
	
	
	// now use the a_c to build the weight matrix
	
	// first column, single nt   
	nl = strlen(a_c[0]);
	
	for (j=0; j<16; j++)
		mywm[0][j] = 0;
	
	for (i=0; i<nl; i++)
		mywm[0][ N[ (int)(a_c[0][i]) ] ] = 120 / nl;
	
	for (i=1; i<myw; i++) {
		nl0 = strlen(a_c[i-1]);
		nl1 = strlen(a_c[i  ]);
		
		//printf("nl0 = %s, nl1 = %s\n", a_c[i-1], a_c[i]);
		
		// init column
		for (j=0; j<16; j++)
			mywm[i][j] = 0;
		
		// 
		for (j=0; j<nl0; j++) {  // the given .. P(.|J)
			for (k=0; k<nl1; k++) {  // P(K|.)
				//printf("%c %c\n", a_c[i-1][j], a_c[i][k]);
				//printf("N1 = %d\n", N1[ (int)(a_c[i-1][j]) ][ (int)(a_c[i][k]) ]);
				// mywm[i][ N1[ (int)(a_c[i-1][j]) ][ (int)(a_c[i][k]) ] ] = (int)(120.0 / (nl0*nl1));
				
				// records P(Xi,Xi-1)
				mywm[i][ N1[ (int)(a_c[i][k]) ][ (int)(a_c[i-1][j]) ] ] = (int)(120.0 / (nl0*nl1));
				
			}
		}
		
	}
	
	
	//printIntegerWM(mywm, myw);
	
	*w  = myw;
	*wm = mywm;
}


//
// convert a 0th-order Markov WM into a first-order Markov model
//  this is done by making all previous nt equally probable (dividing by 4 basically)
//
void convert_wm_0m_to_1m(int** intwm_0m, int w, int*** intwm_1m)
{
	
	int lig;
	int** mywm2;
	int i, j;
	int dint_idx;
	
	
	// allocate memory
	mywm2 = (int**)calloc( w, sizeof(int*) );
	for (i=0; i<w; i++) {
		mywm2[i] = (int*)calloc(16, sizeof(int));
	}
	
	// 0 is unchanged
	for(lig=0; lig < 4; lig++) {
		mywm2[0][lig] = intwm_0m[0][lig];
	}  
	
	// 1+
	for (i=1; i<w; i++) {
		for(lig=0; lig < 4; lig++) {      
			for (j=0; j<4; j++) {	
				dint_idx = N1[ ntoC[lig] ][ ntoC[j] ];
				
				//printf("i=%d, lig=%d, j=%d, dint_idx=%d, intwm_0m[i][lig]=%d\n", i, lig, j, dint_idx, intwm_0m[i][lig]);
				
				mywm2[ i ][ dint_idx ] = intwm_0m[i][lig] / 4;
			}
		}
	}
	
	*intwm_1m = mywm2;
}


// transform a dinucleotide count vector into a P(Xi|Xi-1) matrix to be used for scoring.
void dint_to_scoringWM(int** diwm, int w, int*** scwm) 
{
	
	int lig, col;
	int** mywm2;
	int i, j;
	int P[4];
	int dint_idx;
	
	
	// allocate memory
	mywm2 = (int**)calloc( w, sizeof(int*) );
	for (i=0; i<w; i++) {
		mywm2[i] = (int*)calloc(16, sizeof(int));
	}
	
	// 0 is unchanged
	for(lig=0; lig < 4; lig++) {
		mywm2[0][lig] = diwm[0][lig];
	}  
	
	
	for (col= 1; col < w; col++) {
		
		// calculate P(Xi-1) by summing P(.,Xi-1)
		for (i=0; i<4; i++) {
			P[i] = 0;
			for (j=0; j<4; j++) {
				P[i] += diwm[ col ][ N1[ ntoC[j] ][ ntoC[i] ] ]; 
			}      
			
			//printf("col=%d, P[i=%d] = %d\n", col, i, P[i]);
			
		}
		
		// divide P(Xi,Xi-1)/P(Xi-1)
		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {		
				dint_idx = N1[ ntoC[j] ][ ntoC[i] ];
				
				if (diwm[ col ][ dint_idx ] == 0) 
					mywm2[col][ dint_idx ] = 0;
				else
					mywm2[col][ dint_idx ] = (int)( 120 * diwm[ col ][ dint_idx ] /  (float)P[i] ); 
			}      
		}
		
	}
	
	*scwm = mywm2;
	
	
}



//
//  log score
//
float getScoreOnSeq_1M(float** wml, int w, int* stars, float* bkg, char* seq, int i, int rna, int* strand, int* hasn) 
{
	float score1 = 0.0;
	float score2 = 0.0;
	
	int   j, k;
	int   nt_i, nt_im1, dint_idx;
	
	*hasn = 0;
	
	// exit if sequence to score has Ns
	for (j=0; j<w; j++) {
		//if (stars[j] == 1) {
		if ((seq[i + j] != 'A') && (seq[i + j] != 'T') && (seq[i + j] != 'C') && (seq[i + j] != 'G')) {
			*hasn = 1;
			return 0.0;
		}      
	}
	//}
	
	
	// score first position
	score1 = wml[0][ N[ (int)(seq[i]) ] ];
	//printf(" += %f\n", wml[0][ N[ (int)(seq[i]) ] ]);
	if (rna == 0) {
		
		//printf(" += %f\n", wml[0][ C[ (int)(seq[i+w-1]) ] ]);    
		score2 = wml[0][ C[ (int)(seq[i+w-1]) ] ];
		
	}
	
	// score rest
	for (j=1,k=w-2; j<w; j++,k--) {
		
		nt_i     = seq[ i+j   ];  // nt at i 
		nt_im1   = seq[ i+j-1 ];  // nt at i-1
		dint_idx = N1[nt_i][nt_im1];
		
		//printf(" += %f (j=%d) wml[%d][%d] for %c%c\n", wml[j][dint_idx], j, j, dint_idx, seq[ i+j   ], seq[ i+j-1   ]);
		score1  += wml[j][dint_idx];
		
		if (rna == 0) {
			
			nt_i     = char_compl[ (int)(seq[ i+k   ]) ];  // nt at i 
			nt_im1   = char_compl[ (int)(seq[ i+k+1 ]) ];  // nt at i-1
			dint_idx = N1[nt_i][nt_im1];
			score2  += wml[j][dint_idx];
			
			//printf(" += %f (for nti=%c|ntim1=%c) dint_idx=%d  \n", wml[j][dint_idx], seq[i+k], seq[i+k+1], dint_idx);    
			
			
		}
		
	}
	
	//printf("score1 1M = %f, score2 1M = %f\n", score1, score2);
	
	if ((rna == 1) || (score1 >= score2)) {
		*strand = 1;
		return score1;
	} else {
		*strand = -1;
		return score2;
	}
	
	
}


void printIntegerWM(int** wm, int w)
{
	int i,j;
	for (j=0; j<4; j++) {
		
		for (i=0; i<w; i++) {
			printf("%4d", wm[i][j]);
		}
		printf("\n");
	}
	
	
}

void printIntegerWM_1M(int** wm, int w)
{
	int i,j;
	for (j=0; j<16; j++) {
		
		printf("%4s\t", ntoC1[j]);
		for (i=0; i<w; i++) {
			printf("%4d", wm[i][j]);
		}
		printf("\n");
	}
	
	
}


void freeIntegerWM(int** wm, int w) {
	
	int i;
	for(i=0; i<w; i++) {
		free(wm[i]);
	}
	free(wm);
}



void integerWMtoACE(int** intwm, int w, char*** m) 
{
	int i, j, k, l;
	char** mym;
	
	mym = (char**)malloc(120 * sizeof(char*));
	for (i=0; i<120; i++)
		mym[i] = (char*)malloc(w * sizeof(char));
	
	for (i=0; i<w; i++) {
		
		l = 0;
		for (j=0; j<4; j++) {      
			for (k=0; k<intwm[i][j]; k++) {	
				mym[l][i] = ntoC[ j ];
				l++;
			}
		}
		
	}
	
	*m = mym;
}


//  wm2 = log(wm1)
void integerWMtoLog(int** wm1, float*** wm2, int h, int w) 
{
	int lig, col;
	float pseudo = 0.000000001;
	float** mywm2;
	int i;
	
	// allocate memory
	mywm2 = (float**)calloc( w, sizeof(float*) );
	for (i=0; i<w; i++) {
		mywm2[i] = (float*)calloc(h, sizeof(float));
	}
	
	
	
	for(lig=0; lig < h; lig++) {
		for(col= 0; col < w; col++) {
			mywm2[col][lig] = log( (wm1[col][lig] / 120.0) + pseudo );
		}
	}
	
	*wm2 = mywm2;
	
}

void shuffleRegexpMotif(char* motif, int gap, char** shuffled_motif)
{
	int  nc = 0;
	int* v  = 0;
	int  l,i;
	char** a_re;
	char*  newmotif;
	int  lg = 0, rg = 0;
	int  lnc;
	
	// motif length
	nc  = getRegexpMotifLength(motif, gap); 
	lnc = nc;
	
	//
	// printf("Motif has %d characters.\n", nc);
	// getchar();
	//
	
	l  = strlen(motif);
	
	getCharArrayFromRegexpMotif(motif, gap, &a_re); 
	
	while (a_re[lg][0] == '.') {
		lg++;
	}
	
	i = nc-1;
	while (a_re[i][0] == '.') {
		i--;
		rg++;
	}
	
	
	//printf("%d gaps to the left, %d to the right.\n", lg, rg);
	
	//for (i=0; i<nc; i++) {
	//  printf(" %s\n", a_re[i]);
	//}
	
	nc = nc - lg - rg;
	
	v  = randPermOfIndices(nc); 
	
	newmotif = (char*)calloc(l+1, sizeof(char));
	for (i=0; i<lg; i++) 
		strcat(newmotif, ".");
	
	for (i=0; i<nc; i++) {
		strcat(newmotif, a_re[ lg + v[i] ]);
	}
	
	for (i=0; i<rg; i++) 
		strcat(newmotif, ".");
	
	newmotif[l] = '\0';
	free(v);
	
	for (i=0; i<lnc; i++)
		free(a_re[i]);
	
	free(a_re);
	
	*shuffled_motif = newmotif;
}


void getCharArrayFromRegexpMotif(char* motif, int gap, char*** a_re) 
{
	
	
	int      l = strlen(motif);
	int      i = 0;  // cnt in ungapped motif
	int      j = 0;  // cnt in   gapped motif
	int      k = 0;  // cnt characters
	int      p_b;
	int      nc;
	char**   my_a_re = 0;
	char*    str;
	int      m = 0;
	
	nc = getRegexpMotifLength(motif, gap); 
	
	my_a_re = (char**)malloc(nc * sizeof(char*));
	
	while (i < l) {
		
		str = (char*)calloc(5, sizeof(char));
		m   = 0;
		
		if (motif[i] == '[') {
			
			p_b = i;
			while (motif[p_b] != ']') {	
				p_b++; 
			} 
			
			
			while (i <= p_b) {	
				str[m] = motif[i];
				i++; j++; m++;
			}
			
			
			
		} else {
			str[m] = motif[i];
			j++; i++; m++;
		}  
		
		str[m] = '\0';
		
		//printf("str=%s\n", str);
		my_a_re[k] = str;
		
		k++;
	}
	
	*a_re = my_a_re;
}


int getRegexpMotifLength(char* motif, int gap) 
{
	
	
	int      l = strlen(motif);
	int      i = 0;  // cnt in ungapped motif
	int      j = 0;  // cnt in   gapped motif
	int      k = 0;  // cnt characters
	//int      h = 0;
	int      p_b;
	
	
	while (i < l) {
		
		if (motif[i] == '[') {
			
			p_b = i;
			while (motif[p_b] != ']') {	
				p_b++; 
			} 
			
			while (i <= p_b) {	
				i++; j++;
			}
			
			k++;
		} else {
			j++; i++; k++;
		}    
	}
	
	return k + gap;
}


