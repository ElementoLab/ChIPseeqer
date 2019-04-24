#ifndef SAM_H
#include "lib/third_party/samtools/sam.h"
#include "lib/third_party/samtools/bam.h"
#include "lib/third_party/samtools/faidx.h"
#define SAM_H
#endif

#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)
#define _GNU_SOURCE


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <search.h>
#include <math.h>

#include <limits.h>

#include "hashtable.h"
#include "statistics.h"
#include "dataio.h"
#include "uthash.h"
#include "sequences.h"

#include "readio.h"

#define FREE(p) do { free(p); p = NULL; } while(0)

#define bam1_unmapped(b) (((b)->core.flag & BAM_FUNMAP) != 0)
#define bam1_mate_unmapped(b) (((b)->core.flag & BAM_FMUNMAP) != 0)
#define bam1_paired(b) (((b)->core.flag & BAM_FPAIRED) != 0)



char* my_fai_fetch(faidx_t* idx, char* c, int i, int j)
{
  int len;
  char id[1000];
  sprintf(id, "%s:%d-%d", c, i, j);
  return fai_fetch(idx, id, &len);
}



void parseTopHatCIGR(char* cigr, int** blocksizes, int** blockoffsets, int* numblocks)
{
  int maxnumblocks = 10;
  int i, blocklen, numdigits;
  int cigrlen      = strlen(cigr);
  char* snum       = (char*)calloc(100, sizeof(char));
  *blocksizes      = (int*)calloc(maxnumblocks, sizeof(int));
  *blockoffsets    = (int*)calloc(maxnumblocks, sizeof(int));
    
  int offset = 0;
  numdigits  = 0;
  for (i=0; i<cigrlen; i++) {
    if (cigr[i] == 'M') {
      blocklen  = atoi(snum);
      numdigits = 0;

      // add block
      (*blocksizes)[*numblocks]   = blocklen;
      (*blockoffsets)[*numblocks] = offset;
      offset += blocklen;
      (*numblocks)++;

    } else if (cigr[i] == 'N') {
      blocklen  = atoi(snum);
      numdigits = 0;	      

      // do not add block but update offset
      offset += blocklen;
      
    } else if (isalpha(cigr[i])) {
      printf("Char unrecognized: %c\n", cigr[i]);
      exit(1);

    } else {      
      snum[numdigits] = cigr[i]; // add digits
      numdigits++;
      snum[numdigits] = '\0'; // add end char
    }
  }
  
  free(snum);
}


int guessReadLength(char* readfile, char* format, int head)
{
  FILE* fp   = 0;
  char* buff = 0;
  int   mynmax = 10000;
  int maxl = -1;
  int minl = 100000;
  int l    = -1;
  char** p;
  int m;

  fp = fopen(readfile, "r");
  if (fp == 0) {
    printf("Cannot open %s\n", readfile);
    exit(0);
  }
  
  buff = (char*)calloc(mynmax, sizeof(char));
  
  
  int i = 0;
  while (fgets(buff, mynmax, fp) != 0) {
    
    chomp(buff);
    split_line_delim(buff, "\t", &p, &m);
  
    if (strcmp(format, "sam") == 0) {
      if (m <= 3) {
	free(p);
	continue;
      } else {
	
	l = strlen(p[9]);

	if (l > maxl) 
	  maxl = l;
	if (l < minl)
	  minl = l;
	
      }
    } else {
      die("Not implemented\n");
    }

    

    free(p);
    i++;
    if (i == head) {
      break;
    }
    
  }
  free(buff);
  fclose(fp);

  if (minl != maxl) {
    die("Problem: reads with unequal lengths detected\n");
  }

  return minl;

}

void updateMismatchCountAtReadPos(char* file, int chrlen, int** poscounts, long* numnt, long* nummm, int uniquereads)
{


	FILE* fp = 0;
	int   bufflen = 10000000;
	char* buff;
	int   cnt;

	char* s;
	char** p;
	int   m;
	char* prevs = 0;
	//int   i =0;
	int   k = 0;
	char* mys;
	int   j = 0;

	int   pos = 0;
	int   lenread = 0;
	int   st = 0;
	//char  nt;
	//int   absposmut;
	int   posmut;
	unsigned char* a_unqreads_F = 0;
	unsigned char* a_unqreads_R = 0;
	//unsigned char* co = 0;
	int     posinread = 0;
	//int     rst = 0;
	//int     ren = 0;

	fp = fopen(file, "r");
	if (fp == 0) {
		die("Cannot open read file.\n");
	}

	if (uniquereads == 1) {
		a_unqreads_F = create_binarized_array(chrlen);
		a_unqreads_R = create_binarized_array(chrlen);
	}



	//co = (unsigned char*)malloc(1000 * sizeof(unsigned char));

	j = 0;
	char * buff_st;

	while (1){		// old: !feof(fp)) 

		// allocate data for new chunk
		buff   = (char*)calloc(bufflen, sizeof(char));
		if (!buff) {
			die("Cannot allocate buff in seqI_open_fast\n");
		}
		
		// DEBUG
		// buff start -> so it's not lost after strsep
		buff_st= buff;

		// read data from fp
		cnt = fread(buff, 1, bufflen, fp);
		if (cnt == 0){
			FREE(buff_st);	
			FREE(prevs);
			break;
		}

		//i = 0;

		// split by line, get the first line in s
		s = strsep(&buff, "\n");

		// if this is not the first chunk we read
		if (j != 0) {

			// alloc new s
			mys = (char*)calloc(1000, sizeof(char));

			// concatenate what is left of prev line
			strcat(mys, prevs);

			// concate end of same line
			strcat(mys, s);

			// update s
			s = mys;

			// DEBUG
			FREE(prevs);
		}


		while (1) {

			// store current line
			prevs = s;
			//i++;

			// get next line
			s = strsep(&buff, "\n");
			if (s == 0){
				// DEBUG
				prevs= strdup(prevs);
				break;
			}else {

				// split
				split_line_delim(prevs, "\t", &p, &m);

				//
				// parse ELAND format
				//
				pos = atoi(p[7]);
				if (p[8][0] == 'R')
					st = -1;
				else
					st = 1;

				//
				// if flag set, and read already there, skip that read
				//
				if (uniquereads == 1) {
					if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
					        ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
						FREE(p);
						continue;
					} else {
						if (st == 1)
							set_entry_in_binarized_array(a_unqreads_F, pos);
						else
							set_entry_in_binarized_array(a_unqreads_R, pos);
					}
				}

				lenread = strlen(p[1]);
				(*numnt) += lenread;


				if ((strcmp(p[2], "U1") == 0) || (strcmp(p[2], "U2") == 0)) {

					for (k=10; k<m; k++) {

						p[k][ strlen(p[k]) - 1 ] = '\0';     // replace nt by \0
						posmut  = atoi(  p[ k ] );           // get pos mismatch
						if (st == 1)                         // calculate actual pos in read
							posinread = posmut - 1;
						else
							posinread = lenread - posmut;

						//printf("Inc pos=%d\n", posinread);

						(*poscounts)[ posinread ] ++;  // increase mm count at that pos
						(*nummm) ++;                   // increase overall mm count
					} //for

				} // if


				FREE(p); // p[x] are just pointers


			}
		}

		if (j != 0)
			FREE(mys);

		FREE(buff_st);
		// free(buff);
		j++;
	}

	if (uniquereads == 1) {
		FREE(a_unqreads_F);
		FREE(a_unqreads_R);
	}

}


void getReadCountFromReadFile(char* file, int chrlen, unsigned short** readcounts, unsigned short** mutacounts, unsigned short** mutacounts_F,
                              int* numreads, long* numnt, long* nummm, unsigned char*** ntcounts, unsigned char*** readpos, int uniquereads, int rw, int truncate)
{

	FILE*  fp = 0;
	int    bufflen = 10000000;
	char*  buff;
	int    cnt;
	char*  s;
	char** p;
	int    m;
	char*  prevs = 0;
	int    i =0;
	int    k = 0;
	char*  mys;
	int    j = 0, l;

	int    pos = 0;
	int    lenread = 0;
	int    st = 0;
	char   nt;
	int    absposmut;
	int    posmut;
	unsigned char* a_unqreads_F = 0;
	unsigned char* a_unqreads_R = 0;
	unsigned char* co = 0;
	//int     posinread = 0;
	int     rst = 0, ren = 0;

	fp = fopen(file, "r");
	if (fp == 0) {
		die("Cannot open read file.\n");
	}

	if (uniquereads == 1) {
		a_unqreads_F = create_binarized_array(chrlen);
		a_unqreads_R = create_binarized_array(chrlen);
	}

	*numreads = 0;

	co = (unsigned char*)malloc(1000 * sizeof(unsigned char));

	j = 0;
	char * buff_st;

	while (1){	// !feof(fp)) 

		// allocate data for new chunk
		buff   = (char*)calloc(bufflen, sizeof(char));
		if (!buff) {
			die("Cannot allocate buff in seqI_open_fast\n");
		}

		// read data from fp
		cnt = fread(buff, 1, bufflen, fp);
		buff_st = buff;

		if (cnt == 0){
			// DEBUG
			FREE(prevs);
			FREE(buff_st);
			break;
		}

		//i = 0;

		// split by line, get the first line in s
		s = strsep(&buff, "\n");

		// if this is not the first chunk we read
		if (j != 0) {

			// alloc new s
			mys = (char*)calloc(1000, sizeof(char));

			// concatenate what is left of prev line
			strcat(mys, prevs);

			// concate end of same line
			strcat(mys, s);

			// update s
			s = mys;


			// DEBUG
			FREE(prevs);
		}


		while (1) {

			// store current line
			prevs = s;
			//i++;

			// get next line
			s = strsep(&buff, "\n");
			if (s == 0){
				// DEBUG
				prevs= strdup(prevs);
				break;
			}else {

				// split
				split_line_delim(prevs, "\t", &p, &m);

				//
				// parse ELAND format
				//
				pos = atoi(p[7]);
				if (p[8][0] == 'R')
					st = -1;
				else
					st = 1;

				//
				// if flag set, and read already there, skip that read
				//
				if (uniquereads == 1) {
					if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
					        ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
						FREE(p);
						continue;
					} else {
						if (st == 1)
							set_entry_in_binarized_array(a_unqreads_F, pos);
						else
							set_entry_in_binarized_array(a_unqreads_R, pos);
					}
				}

				lenread   = strlen(p[1]);
				(*numnt) += lenread;

				rst = pos + truncate;
				ren = pos + lenread - truncate;
				int posinread = 0;
				for (i=pos,k=0; i<pos+lenread; i++,k++) {
					co[k] = 0;
					if ((i >= rst) && (i < ren) && ((*readcounts)[i-1] < USHRT_MAX)) {

						if ((*readpos)[i-1] != 0) {

							// calculate where in read we are
							if (st == 1)
								posinread = k;
							else
								posinread = lenread - k - 1;

							// add pos
							(*readpos)[i-1][  (*readcounts)[i-1] ] = (char)posinread;

						}

						// add read
						(*readcounts)[i-1] ++;
					}
				}

				if ((strcmp(p[2], "U1") == 0) || (strcmp(p[2], "U2") == 0)) {

					// first pass to calculate neighborhoods
					for (k=10; k<m; k++) {
						p[k][ strlen(p[k]) - 1 ] = '\0';     // replace nt by \0
						posmut  = atoi(  p[ k ] );            // get pos mismatch
						if (st == 1)
							posinread = posmut - 1;
						else
							posinread = lenread - posmut;

						for (l=max(0,posinread-rw); l<=min(lenread-1,posinread+rw); l++)
							co[l] ++;
					}

					// now count mismatches
					for (k=10; k<m; k++) {

						posmut  = atoi(  p[ k ] );            // get pos mismatch
						if (st == 1)
							posinread = posmut - 1;
						else
							posinread = lenread - posmut;

						if ((co[posinread] == 1) && (posmut <= lenread-truncate) && (posmut >= truncate+1)) {

							nt = p[1][ posinread ];     // get mismatch nt in reads
							if (st == -1)
								nt = char_complement(nt);

							if (nt != 'N') {

								absposmut =  pos - 1 + posmut - 1;
								(*mutacounts)[ absposmut ] ++;

								if (st == 1)
									(*mutacounts_F)[ absposmut ] ++;

								if ((*ntcounts) != 0) {
									if ( (*ntcounts)[absposmut] == 0) {
										(*ntcounts)[absposmut] = (unsigned char*)calloc(4, sizeof(unsigned char));
										if ( (*ntcounts)[absposmut] == 0 )
											die("problem allocaing *ntcounts[]\n");
									}
									if ( (*ntcounts)[absposmut][ data_nucl(nt) ] < UCHAR_MAX)
										(*ntcounts)[absposmut][ data_nucl(nt) ] ++;

								}
								(*nummm) ++;

							} // if (nt != N)
						}

					}

				}

				(*numreads) ++;

				FREE(p); // p[x] are just pointers


			}
		}

		if (j != 0)
			FREE(mys);

		FREE(buff_st);
		j++;
	}

	FREE(co);

	if (uniquereads == 1) {
		FREE(a_unqreads_F);
		FREE(a_unqreads_R);
	}


}



void getSimpleReadCountFromReadFile(char* file, int chrlen, unsigned short** readcounts, int* numreads, int uniquereads, int truncate)
{

	FILE* fp = 0;
	int   bufflen = 10000000;
	char* buff;
	int   cnt;

	char* s;
	char** p;
	int   m;
	char* prevs = 0;
	int   i =0;
	char* mys;
	int   j = 0;

	int   pos = 0;
	int   lenread = 0;
	int   st = 0;
	//  char  nt;
	//  int   absposmut;
	//  int   posmut;
	unsigned char* a_unqreads_R = 0;
	unsigned char* a_unqreads_F = 0;

	fp = fopen(file, "r");
	if (fp == 0) {
		die("Cannot open read file.\n");
	}

	if (uniquereads == 1) {
		a_unqreads_R = create_binarized_array(chrlen);
		a_unqreads_F = create_binarized_array(chrlen);
	}

	*numreads = 0;

	j = 0;
	char * buff_st;

	while (1){	// !feof(fp)) 

		// allocate data for new chunk
		buff   = (char*)calloc(bufflen, sizeof(char));
		if (!buff) {
			die("Cannot allocate buff in seqI_open_fast\n");
		}

		// read data from fp
		cnt = fread(buff, 1, bufflen, fp);
		buff_st= buff;

		if (cnt == 0){
			FREE(buff_st);
			FREE(prevs);
			break;
		}

		i = 0;

		// split by line, get the first line in s
		s = strsep(&buff, "\n");

		// if this is not the first chunk we read
		if (j != 0) {

			// alloc new s
			mys = (char*)calloc(1000, sizeof(char));

			// concatenate what is left of prev line
			strcat(mys, prevs);

			// concate end of same line
			strcat(mys, s);

			// update s
			s = mys;

			FREE(prevs);
		}


		while (1) {

			// store current line
			prevs = s;
			i++;

			// get next line
			s = strsep(&buff, "\n");
			if (s == 0){
				// DEBUG
				prevs= strdup(prevs);
				break;
			}else {

				// split
				split_line_delim(prevs, "\t", &p, &m);

				//
				// parse ELAND format
				//
				pos = atoi(p[7]);
				if (p[8][0] == 'R')
					st = -1;
				else
					st = 1;


				//
				// if flag set, and read already there, skip that read
				//
				if (uniquereads == 1) {
					if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
					        ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
						FREE(p);
						continue;
					} else {
						if (st == 1)
							set_entry_in_binarized_array(a_unqreads_F, pos);
						else
							set_entry_in_binarized_array(a_unqreads_R, pos);
					}
				}

				lenread = strlen(p[1]);

				for (i=pos+truncate; i<pos+lenread-truncate; i++) {
					if ((*readcounts)[i-1] < USHRT_MAX) {
						(*readcounts)[i-1] ++;
					}
				}

				(*numreads) ++;

				FREE(p); // p[x] are just pointers


			}
		}


		FREE(buff_st);

		if (j != 0)
			FREE(mys);

		j++;
	}

	if (uniquereads == 1) {
		FREE(a_unqreads_F);
		FREE(a_unqreads_R);
	}

}






void getINDELCountFromSAMFile(const char* file, const int chrlen, unsigned short** readcounts, unsigned short** mutacounts, unsigned short** mutacounts_F,
			     int* numreads, long* numnt, long* nummm, int del, unsigned short** indellen, int uniquereads, const int readlen, int uniqmap, int verbose)
{

  char** p;
  int    m;
  int    i =0;
  //  int    k = 0;
  //int    j = 0;
  int    pos = 0;
  int    st = 0;
  //char   nt;
  //int    absposmut;
  unsigned char* a_unqreads_F = 0;
  unsigned char* a_unqreads_R = 0;
  //unsigned char* co = 0;
  //int     rst = 0, ren = 0;
  //int     verbose = 0;
  char*   snum = 0;
  int     bitflag;
  //char* mmstr;
  //int   mmstrlen = 0;
  //char* readseq;
  //int   readseqlen;
  int   lb;
  char* cigr;
  //char* alignedseq;
  //char* myreadpos;
  //char* alignedpos;
  //int   idxnewread = 0;
  //int   idxread    = 0;
  int   cigrlen = 0;
  //int   inum = 0;
  //int   idxalnread = 0;
  //int   insertion = 0;
  //int   alignedseqlen = 0;
  //int   rennotr;
  //int   ntidx = 0;
  //int   l;
  //int cline = 0;
  //char * buff_st;
  char* line;
  LineI li;
  //int idxmd = -1;
  int posmut = 0;
  int numdigits = 0;


  snum = (char*)calloc(100, sizeof(char));

  if (uniquereads == 1) {
    a_unqreads_F = create_binarized_array(chrlen);
    a_unqreads_R = create_binarized_array(chrlen);
  }

  // init read counter
  *numreads = 0;

  // read in 
  li.bufflen = 100000000;
  LineIopen(&li, (char*)file);
  while ((line = nextLine(&li))) {

   	  	 	  
    // split
    split_line_delim(line, "\t", &p, &m);
	  
    // skip if not aligned
    if (strcmp(p[5], "*") == 0) {
      FREE(p);
      free(line);
      continue;
    }

    // non-unique mapper
    if ((uniqmap == 1) && (strcmp(p[11], "XT:A:U") != 0)) {
      FREE(p);
      free(line);
      continue;
    }

    // read pos (starts at 1)
    pos = atoi(p[3]);
	  
    // strand
    bitflag = atoi(p[1]);

    // unmapped: 
    if ( bitflag & ( 1 << 2 ) ){
      //fprintf(stderr, "WARNING: skipping unmapped read: %s %s %s\n", p[0], p[1], p[2]);
      free(p);
      free(line);
      continue;
    }


    // check bits 4 an 5 (from 0)
    if ((bitflag & (1 << 4)) ||
	(bitflag & (1 << 5))) {
      //if (bitflag == 16) 
      st = -1;
    } else {
      st = 1;
    }
	  
    //
    // if flag set, and read already there, skip that read
    //
    if (uniquereads == 1) {
      if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
	  ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
	FREE(p);
	FREE(line);
	continue;
      } else {
	if (st == 1)
	  set_entry_in_binarized_array(a_unqreads_F, pos);
	else
	  set_entry_in_binarized_array(a_unqreads_R, pos);
      }
    }
	  
    cigr       = strdup(p[5]);   // CIGAR
    numdigits  = 0;              // num digits
    cigrlen    = strlen(cigr);
    posmut     = pos - 1;                // indel position (NOTE THE -1; THIS IS BECAUSE SAM FORMAT IS 1-BASED)

    // parse CIGAR
    for (i=0; i<cigrlen; i++) {
	    
      if (cigr[i] == 'M') {
	
	// length block
	lb = atoi(snum);
	
	// update idx
	posmut += lb;
	      
	// reset digits
	numdigits = 0;  
	      
      } else if ( cigr[i] == 'N' ) {

	printf("N in CIGR (%s)?\n", cigr);
	exit(0);
	lb= atoi(snum);
	
      } else if (cigr[i] == 'D') {

	// deleting the next lb bases from reference genome
	lb = atoi(snum);

	// register deletion if we care
	if (del == 1) {
	  
	  //printf("DELETION in TS, pos %d !\n", posmut);
	  (*mutacounts)[ posmut ] ++;
	  if((*indellen)[posmut] < lb)
	    {
	      (*indellen)[posmut] = lb;
	    }
	  (*nummm)++;
	  
	  if (st == 1)
	    (*mutacounts_F)[ posmut ] ++; 

	}
	
	// increment pos in genome by block length
	posmut += lb;
	      
	numdigits = 0; // reset digits
	      
      } else if (cigr[i] == 'I') {

	// block length
	lb = atoi(snum);

		// register deletion if we care
	if (del == 0) {
	  
	  //printf("DELETION in TS, pos %d !\n", posmut);
	  (*mutacounts)[ posmut ] ++;
	  if((*indellen)[posmut] < lb)
	    {
	      (*indellen)[posmut] = lb;
	    }
	  (*nummm)++;
	  
	  if (st == 1)
	    (*mutacounts_F)[ posmut ] ++; 

	}
	      
	numdigits = 0; // reset digits
	      
      } else if (cigr[i] == 'S') {

	//printf("S in CIGR (%s)?\n", cigr);
	//exit(0);
	// ignore 
	lb = atoi(snum);
	numdigits = 0; // reset digits
	      
      } else {
	// num
	snum[numdigits] = cigr[i]; // add digits
	numdigits++;
	snum[numdigits] = '\0'; // add end char
      }
    }

    // aligned read length
    int lenread =  posmut - pos;
    (*numnt)    += lenread;
    for (i=pos; i<min(chrlen,pos+lenread); i++) {
      if ((*readcounts)[i-1] < USHRT_MAX) {
	(*readcounts)[i-1] ++;
      }
    }

    (*numreads)++;

    FREE(p);
    FREE(cigr);
    FREE(line);

    if ((verbose == 1) && ( ((*numreads) % 100000) == 0)) {
      printf("Read %d reads from %s\n", *numreads, file);
    }
  }

  if (verbose == 1) 
    printf("Number of reads = %d\n", *numreads);
	
  LineIclose(&li);

  if (uniquereads == 1) {
    FREE(a_unqreads_F);
    FREE(a_unqreads_R);
  }

  if (verbose == 1) 
    printf("About to exit after reading %s\n", file);

  FREE(snum);

  
}






//
//
//  SAM format
//	DEBUG: adding support for extended CIGAR strings (i.e. incorporate the N operation)
//
//  returns/changes:
//	readcounts (along the chromosome/transcript length)
//	numreads
//

// removed: char * seq -> never used
void getSimpleReadCountFromReadFile_SAM(const char* file, const int chrlen, unsigned short** readcounts, int* numreads, const int uniquereads, const int truncate, const int readlen)
{


  //char* buff;
  //int   cnt;

  //char* s;
  char** p;
  int   m;
  //char* prevs = 0;
  int   i =0;
  int   k = 0;
  //char* mys = 0;
  int   j = 0; //, l;

  int   pos = 0;
  int   st = 0;
  unsigned char* a_unqreads_F = 0;
  unsigned char* a_unqreads_R = 0;
  int     rst = 0; //, ren = 0;

  //int    verbose = 1;

  char* snum = (char*)calloc(100, sizeof(char));
  int   bitflag;
  char* mmstr;
  char* readseq;
  int   readseqlen;
  int   lb;
  char* cigr;
  //char* alignedseq;
  //char* readpos;
  //char* alignedpos;
  int   idxnewread = 0;
  int   idxread    = 0;
  int   cigrlen = 0;

  // new version of aligned seq
  my_aln *aln_seq= NULL, *aln_nt= NULL;

  //DEBUG
  // printf("%s, chrlen= %d\n", file, chrlen);
  if (uniquereads == 1) {
    // a binarized array of length chrlen (chromosome length)
    // DEBUG:
    //	maybe at some point p in create_binarized_array() is given up so a_unqreads_F/R is no longer safe
    a_unqreads_F = create_binarized_array(chrlen);
    a_unqreads_R = create_binarized_array(chrlen);
    // create_safe_binarized_array(chrlen, &a_unqreads_F);
    // create_safe_binarized_array(chrlen, &a_unqreads_R);
  }

  *numreads = 0;
  j = 0;

  char* line;
  LineI li;
  li.bufflen = 10000000;
  LineIopen(&li, (char*)file);
  while ((line = nextLine(&li))) {
	  
    char* myline = strdup(line);

    // split
    split_line_delim(line, "\t", &p, &m);
	  
    // get mm string
    mmstr = p[m - 1] + 5;
	  

    // skip if not aligned
    if (strcmp(p[5], "*") == 0) {
      FREE(line); FREE(myline); FREE(p);
      continue;
    }
	  
    // read pos (starts at 1)
    pos = atoi(p[3]);
	  
    // strand
    bitflag = atoi(p[1]);
	  
    // unmapped:  
    // 	0x0004: the query sequence itself is unmapped
    if ( bitflag & ( 1 << 2 ) ){
    //if ( bitflag & 4 != 0 ){      
      //fprintf(stderr, "WARNING: skipping unmapped read: %s %s %s\n", p[0], p[1], p[2]);
      FREE(line); FREE(myline); FREE(p);
      continue;
    }
    // strand: check bits 4 and 5 (from 0)
    if ((bitflag & (1 << 4)) ||
	(bitflag & (1 << 5))) {
      //if (bitflag == 16)
      st = -1;
      //if (verbose == 1)
      //	printf("STRAND = +++++++++++\n");
    } else {
      st = 1;
      //if (verbose == 1)
      //	printf("STRAND = -----------\n");
	    
    }
    //
    // if flag set, and read already there, skip that read
    //
    if (uniquereads == 1) {
      if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
	  ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
	// i.e. position has a read already
	//if (verbose == 1)
	//	printf("we already have %d, pos= %d\n", st, pos);
	      
	FREE(line); FREE(myline); FREE(p);
	continue;
      } else {
	if (st == 1){
	  set_entry_in_binarized_array(a_unqreads_F, pos);
	}else{
	  set_entry_in_binarized_array(a_unqreads_R, pos);
	}	
      }
    }
	  
    // actual read sequence
    readseq    = strdup(p[9]);
    readseqlen = strlen(readseq);
	  
    // parse CIGAR
	  
    cigr       = strdup(p[5]);
    // printf("readseq= %s, len= %d, cigar= %s\n", readseq, readseqlen, cigr);
	  
    // cigr       = p[5];
	  
    idxnewread = 0;
    idxread    = 0;
	  
    k          = 0; // num digits
    cigrlen    = strlen(cigr);
	  
    // printf("cigar= %s\n", cigr);
    for (i=0; i<cigrlen; i++) {
      //if (verbose)
      //	printf("snum= %s\n", snum);
	    
      if (cigr[i] == 'M') {
	lb = atoi(snum);
	      
	// DEBUG
	// readseq= what's in the SAM file
	// alignedseq= what's reconstructed
	//	-> should use a hash table
	//	-> alignedpos is NOT USED -> maybe just drop it.
	//if (verbose) {
	//printf("idxnewread= %d, idxread= %d, lb= %d\n", idxnewread, idxread, lb);
	// printf("before: %s\n", alignedseq);
	      
	// printf("adding %s to alignedseq position %d\n", strndup(readseq+idxread, lb), idxnewread);
	//}
	      
	// DEBUG: new hash based alignedseq ... add a block of nucleotides
	int tmpidx;
	for (tmpidx= 0; tmpidx< lb; tmpidx++) {
	  if ( (aln_nt = (my_aln*)malloc(sizeof(my_aln))) == NULL) {
	    fprintf(stderr, "no memory left!\n");
	    exit(-1);
	  }
	  aln_nt->tpos= idxnewread + tmpidx;
	  aln_nt->nt = (char)readseq[idxread + tmpidx];
		
	  // if (verbose)
	  //	printf("adding %c to pos %d\n",  aln_nt->nt, aln_nt->tpos);
	  HASH_ADD_INT(aln_seq, tpos, aln_nt);
	}
	      
	// reconstructing aligned sequence on the genome
	      
	// update both idx
	idxnewread += lb;
	idxread    += lb;
	      
	k = 0;  // reset digits
	      
      } else if (cigr[i] == 'N') {
	// DEBUG: new operation: N= skipped from reference= gap between exons
	lb= atoi(snum);
	// printf("we have N: %d\n", lb);
	idxnewread += lb;
	k= 0;
	      
      } else if (cigr[i] == 'D') {
	// deletion from reference detected, add one - to alignedseq
	lb = atoi(snum);
	      
	// add gaps
	// update index in new read
	idxnewread += lb;
	      
	k = 0; // reset digits
	      
      } else if (cigr[i] == 'I') {
	lb = atoi(snum);

	// update idx in read
	idxread += lb;
	      
	k = 0; // reset digits
	      
      } else if (cigr[i] == 'S') {
	lb = atoi(snum);
	k = 0; // reset digits
	      
      } else if (isalpha(cigr[i])) {
	fprintf(stderr, "WARNING: CIGAR string format (%s) not recognized!\n", cigr);

	k= 0;
      } else {
	// num
	snum[k] = cigr[i]; // add digits
	k++;
	snum[k] = '\0'; // add end char
      }
    } //end for cigrlen
	  
    // print keys of aln_seq -> check if they are in ascending order (no, but does not matter)
    // also tally the nucleotides
    my_aln *sss;
	  
    rst = pos + truncate; // truncate not used anymore (=0)
    for(sss=aln_seq; sss != NULL; sss=sss->hh.next) {
      // printf("position in hash: %d, nt: %c, pos= %d, truncate= %d\n", sss->tpos, sss->nt, pos, truncate);
      i= rst+ sss->tpos;		// location on the genome
	    
      if (i> chrlen-1){
	fprintf(stderr, "getSimpleReadCountFromReadFile_SAM: ERROR! inconsistent chromosome/sequence length in %s start= %s, position %d (length on file= %d)!!\n", 
		file, p[3], i-1, chrlen);
	continue;
	      
      }
      if (((*readcounts)[i-1] < USHRT_MAX) && ( sss->nt != '-')) {
	      
	// BUG HERE
	// 	why adding to position 2739 when mRNA length is only 2700???
	(*readcounts)[i-1] ++;
	      
	// printf("hash: adding one count to position: %d, %c\n", i-1, sss->nt);
      }
    }
	  
    // printf("\n");
	  
    // tally up nucleotides
	  
    (*numreads) ++;
	  
    FREE(p); // p[x] are just pointers
	  
    FREE(cigr);	// don't forget strdup
    FREE(readseq);	
	  
    // reset the hash table aln_seq
    my_aln *current_base;
	  
    while(aln_seq) {
      current_base = aln_seq;          /* copy pointer to first item     */
      HASH_DEL(aln_seq, current_base);  /* delete; users advances to next */
      FREE(current_base);            /* optional- if you want to free  */
    }
	  
    FREE(line); FREE(myline);
  }  
    // else (line exists)

  LineIclose(&li);


	
  if (uniquereads == 1) {
    FREE(a_unqreads_F);
    FREE(a_unqreads_R);
  }
	
  FREE(snum);
}


//
// quality scores
//   qualscores is chrlen * [A,C,G,T] (4) * maxnumreads
 
void getQualityScoresFromReadFile_SAM(const char* file, char* chrname, const int chrlen, unsigned char* snvpos, unsigned char**** qualscores, unsigned short*** numqualscores, HASH* hsnv, const int uniquereads)
{
  char** p;
  int   m;
  int   i                     = 0;
  int   k                     = 0;
  int   j                     = 0;
  int   pos                   = 0;
  int   st                    = 0;
  unsigned char* a_unqreads_F = 0;
  unsigned char* a_unqreads_R = 0;
  int     rst                 = 0; 
  char* snum                  = (char*)calloc(100, sizeof(char));
  int   bitflag;
  char* mmstr;
  char* readseq;
  char* qscores;
  int   readseqlen;
  int   lb;
  char* cigr;
  int   idxnewread = 0;
  int   idxread    = 0;
  int   cigrlen    = 0;
  int   numqs      = 0;
  int   maxnumqs   = 255; //SHRT_MAX;

  /*
  // qualscores
  *qualscores = (unsigned char***)calloc(chrlen, sizeof(unsigned char**));
  if (*qualscores == 0) {
    die("Cannot allocate *qualscores\n");
  }
  for (i=0; i<chrlen; i++) {
    (*qualscores)[i] = (unsigned char**)calloc(4, sizeof(unsigned char*));
    if ((*qualscores)[i] == 0) {
      die("Cannot allocate (*qualscores)[i]\n");
    }
  }


  // num qualscores
  *numqualscores = (unsigned short**)calloc(chrlen, sizeof(unsigned short*));
  if (*numqualscores == 0) {
    die("Cannot allocate *numqualscores\n");
  }
  for (i=0; i<chrlen; i++) {
    (*numqualscores)[i] = (unsigned short*)calloc(4, sizeof(unsigned short));
    if ((*numqualscores)[i] == 0) {
      die("Cannot allocate (*numqualscores)[i]\n");
    }
  }
  */

  // new version of aligned seq
  my_aln *aln_seq= NULL, *aln_nt= NULL;

  if (uniquereads == 1) {
    // a binarized array of length chrlen (chromosome length)
    // DEBUG:
    //	maybe at some point p in create_binarized_array() is given up so a_unqreads_F/R is no longer safe
    a_unqreads_F = create_binarized_array(chrlen);
    a_unqreads_R = create_binarized_array(chrlen);
    // create_safe_binarized_array(chrlen, &a_unqreads_F);
    // create_safe_binarized_array(chrlen, &a_unqreads_R);
  }

  //*numreads = 0;
  j         = 0;

  char* line;
  LineI li;
  li.bufflen = 10000000;
  LineIopen(&li, (char*)file);
  while ((line = nextLine(&li))) {
	  
    char* myline = strdup(line);

    // split
    split_line_delim(line, "\t", &p, &m);
	  
    // get mm string
    mmstr = p[m - 1] + 5;
	  
    // skip if not aligned
    if (strcmp(p[5], "*") == 0) {
      FREE(line); FREE(myline); FREE(p);
      continue;
    }
	  
    // read pos (starts at 1)
    pos = atoi(p[3]);
	  
    // strand
    bitflag = atoi(p[1]);
	  
    // unmapped:  
    // 	0x0004: the query sequence itself is unmapped
    if ( bitflag & ( 1 << 2 ) ){
    //if ( bitflag & 4 != 0 ){      
      //fprintf(stderr, "WARNING: skipping unmapped read: %s %s %s\n", p[0], p[1], p[2]);
      FREE(line); FREE(myline); FREE(p);
      continue;
    }
    // strand: check bits 4 and 5 (from 0)
    if ((bitflag & (1 << 4)) ||
	(bitflag & (1 << 5))) {
      //if (bitflag == 16)
      st = -1;
      //if (verbose == 1)
      //	printf("STRAND = +++++++++++\n");
    } else {
      st = 1;
      //if (verbose == 1)
      //	printf("STRAND = -----------\n");
	    
    }
    //
    // if flag set, and read already there, skip that read
    //
    if (uniquereads == 1) {
      if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
	  ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
	// i.e. position has a read already
	//if (verbose == 1)
	//	printf("we already have %d, pos= %d\n", st, pos);
	      
	FREE(line); FREE(myline); FREE(p);
	continue;
      } else {
	if (st == 1){
	  set_entry_in_binarized_array(a_unqreads_F, pos);
	}else{
	  set_entry_in_binarized_array(a_unqreads_R, pos);
	}	
      }
    }
	  
    // actual read sequence
    readseq    = strdup(p[9]);
    readseqlen = strlen(readseq);
    
    // qscores
    qscores    = strdup(p[10]);

    // parse CIGAR
	  
    cigr       = strdup(p[5]);
	  
    idxnewread = 0;
    idxread    = 0;
	  
    k          = 0; // num digits
    cigrlen    = strlen(cigr);
	  
    for (i=0; i<cigrlen; i++) {
	    
      if (cigr[i] == 'M') {
	lb = atoi(snum);
	      
	// DEBUG: new hash based alignedseq ... add a block of nucleotides
	int tmpidx;
	for (tmpidx= 0; tmpidx< lb; tmpidx++) {
	  if ( (aln_nt = (my_aln*)malloc(sizeof(my_aln))) == NULL) {
	    fprintf(stderr, "no memory left!\n");
	    exit(-1);
	  }
	  aln_nt->tpos= idxnewread + tmpidx;
	  aln_nt->nt = (char)readseq[idxread + tmpidx];
	  aln_nt->qs = (char)qscores[idxread + tmpidx];

	  // if (verbose)
	  //	printf("adding %c to pos %d\n",  aln_nt->nt, aln_nt->tpos);
	  HASH_ADD_INT(aln_seq, tpos, aln_nt);
	}
	      
	// reconstructing aligned sequence on the genome
	      
	// update both idx
	idxnewread += lb;
	idxread    += lb;
	      
	k = 0;  // reset digits
	      
      } else if (cigr[i] == 'N') {
	// DEBUG: new operation: N= skipped from reference= gap between exons
	lb= atoi(snum);
	// printf("we have N: %d\n", lb);
	idxnewread += lb;
	k= 0;
	      
      } else if (cigr[i] == 'D') {
	// deletion from reference detected, add one - to alignedseq
	lb = atoi(snum);
	      
	// add gaps
	// update index in new read
	idxnewread += lb;
	      
	k = 0; // reset digits
	      
      } else if (cigr[i] == 'I') {
	lb = atoi(snum);

	// update idx in read
	idxread += lb;
	      
	k = 0; // reset digits
	      
      } else if (cigr[i] == 'S') {
	lb = atoi(snum);
	k = 0; // reset digits
	      
      } else if (isalpha(cigr[i])) {
	fprintf(stderr, "WARNING: CIGAR string format (%s) not recognized!\n", cigr);

	k= 0;
      } else {
	// num
	snum[k] = cigr[i]; // add digits
	k++;
	snum[k] = '\0'; // add end char
      }
    } //end for cigrlen
	  
    // print keys of aln_seq -> check if they are in ascending order (no, but does not matter)
    // also tally the nucleotides
    my_aln *sss;
    char nt;
    int ntidx ;
    rst = pos; // truncate not used anymore (=0)
    for(sss=aln_seq; sss != NULL; sss=sss->hh.next) {
      // printf("position in hash: %d, nt: %c, pos= %d, truncate= %d\n", sss->tpos, sss->nt, pos, truncate);
      i= rst+ sss->tpos;		// location on the genome
	    
      // code addded to limit memory usage ! skip if no SNV at that pos
      if (get_entry_in_binarized_array(snvpos, i-1) == 0)      
	continue;

      if (i> chrlen-1){
	fprintf(stderr, "getSimpleReadCountFromReadFile_SAM: ERROR! inconsistent chromosome/sequence length in %s start= %s, position %d (length on file= %d)!!\n", 
		file, p[3], i-1, chrlen);
	continue;
	      
      }

      nt = sss->nt;
      ntidx = data_nucl(nt);
      
      if ((ntidx != -1) && ( sss->nt != '-')) {

	// get SNV idx in
	int snvidx = -1;
	char* key = (char*)calloc(1000, 1);
	sprintf(key, "%s-%d", chrname, i-1); // notice the i-1
	HASH_find(hsnv, key, &snvidx);
	if (snvidx == -1) {
	  printf("Problem ... %s not found\n", key);
	  exit(1);
	}
	free(key);

	numqs = (int)((*numqualscores)[snvidx][ntidx]);
	if (numqs == 0) {
	  (*qualscores)[snvidx][ntidx] = (unsigned char*)calloc(maxnumqs, sizeof(unsigned char));
	  if ((*qualscores)[snvidx][ntidx] == 0) {
	    die("cannot allocate mem for (*qualscores)[i-1][ntidx] !!\n");
	  }
	}
	if (numqs < maxnumqs-1) {
	  (*qualscores)[snvidx][ntidx][numqs] = sss->qs;
	  (*numqualscores)[snvidx][ntidx] ++;
	}
	      
	// printf("hash: adding one count to position: %d, %c\n", i-1, sss->nt);
      }
    } // for loop over read pos
	  
    // printf("\n");
	  
    // tally up nucleotides
	  
    //(*numreads) ++;
	  
    FREE(p); // p[x] are just pointers
    FREE(qscores);
    FREE(cigr);	// don't forget strdup
    FREE(readseq);	
	  
    // reset the hash table aln_seq
    my_aln *current_base;
	  
    while(aln_seq) {
      current_base = aln_seq;          /* copy pointer to first item     */
      HASH_DEL(aln_seq, current_base);  /* delete; users advances to next */
      FREE(current_base);            /* optional- if you want to free  */
    }
	  
    FREE(line); FREE(myline);
  }  
    // else (line exists)

  LineIclose(&li);


	
  if (uniquereads == 1) {
    FREE(a_unqreads_F);
    FREE(a_unqreads_R);
  }
	
  FREE(snum);
}




// more complicated than getSimpleReadCOunt...
// returns:
//	mutacounts
//	mutacounts_F
// removed:
//	seq (since it's not used)
void getReadCountFromReadFile_SAM(const char* file, const int chrlen, unsigned short** readcounts, unsigned short** mutacounts, unsigned short** mutacounts_F,
                                  int* numreads, long* numnt, long* nummm, unsigned char*** ntcounts, unsigned char*** readpos, int uniquereads, int rw, int truncate, const int readlen)
{

  char** p;
  int   m;
  int   i =0;
  int   k = 0;
  int   j = 0;
  int   pos = 0;
  int   st = 0;
  char  nt;
  int            absposmut;
  unsigned char* a_unqreads_F = 0;
  unsigned char* a_unqreads_R = 0;
  unsigned char* co = 0;
  int            rst = 0, ren = 0;
  int            verbose = 0;

  // printf("we are in getReadCountFromReadFile_SAM()\n");

  char* snum = (char*)calloc(100, sizeof(char));
  int   bitflag;
  char* mmstr;
  int   mmstrlen = 0;
  char* readseq;
  int   readseqlen;
  int   lb;
  char* cigr;
  char* alignedseq;
  char* myreadpos;
  char* alignedpos;
  int   idxnewread = 0;
  int   idxread    = 0;
  int   cigrlen = 0;
  int   inum = 0;
  int   idxalnread = 0;
  int   insertion = 0;
  int   alignedseqlen = 0;
  int   rennotr;
  int   ntidx = 0;
  int   l;
  //int cline = 0;
  //char * buff_st;
  char* line;
  LineI li;
  int idxmd = -1;

  // int   todebug= 0;

  // DEBUG
  // int   bugposition= 129501800; // 204703056;

  alignedseq = (char*)calloc(readlen * 2, sizeof(char));
  myreadpos  = (char*)calloc(readlen * 2, sizeof(char));
  alignedpos = (char*)calloc(readlen * 2, sizeof(char));

  // DEBUG
  my_aln *aln_seq= NULL, *aln_nt= NULL;
  my_aln *readpos2nt= NULL, *tmpread= NULL;

  //fp = fopen(file, "r");
  //if (fp == 0) {
  //  printf("Cannot open read file %s.\n", file);
  //  exit(1);
  //}

  if (uniquereads == 1) {
    a_unqreads_F = create_binarized_array(chrlen);
    a_unqreads_R = create_binarized_array(chrlen);
  }

  *numreads = 0;

  // DEBUG: 
  // what is co? length= 1000
  co = (unsigned char*)malloc(1000 * sizeof(unsigned char));

  j = 0;

  li.bufflen = 10000000;
  LineIopen(&li, (char*)file);
  while ((line = nextLine(&li))) {
	  	 	  
	  	  
    //char* line = strdup(prevs);
    char* myline = strdup(line);

    // split
    split_line_delim(line, "\t", &p, &m);
	  
    // skip if not aligned
    if (strcmp(p[5], "*") == 0) {
      FREE(p);
      free(line);
      free(myline);
      continue;
    }
    
    // get mm string
    // DEBUG: m is varying. so should check for the tag MD:Z:*
    //if (idxmd == -1) {
    idxmd      = -1;
    int tmpidx = m-1;
    while (tmpidx >= 0) {
      if (strstr(p[tmpidx], "MD:Z") != 0) {
	idxmd = tmpidx;
	break;
      }
      tmpidx --;
    }
    //}
	      
    if (idxmd == -1) {
      die("Cannot find MD:Z tag, check format.\n");
    }

    // useless
    if (strstr(p[idxmd], "MD:Z")== NULL){
      //printf("ERROR! no MD tag found!!!\n");
      printf("ERROR! no MD tag found!!! (last tag is %s, idxmd=%d)\n", p[idxmd], idxmd);
      printf("LINE: %s\n", myline);
    }
    
    mmstr = p[idxmd] + 5;

    // read pos (starts at 1)
    pos = atoi(p[3]);
	  
    // strand
    bitflag = atoi(p[1]);

    // unmapped: 
    if ( bitflag & ( 1 << 2 ) ){
	    //fprintf(stderr, "WARNING: skipping unmapped read: %s %s %s\n", p[0], p[1], p[2]);
      FREE(p);
      free(line);
      free(myline);
      continue;
    }

    // check bits 4 an 5 (from 0)
    if ((bitflag & (1 << 4)) ||
	(bitflag & (1 << 5))) {
      //if (bitflag == 16) 
      st = -1;
    } else {
      st = 1;
    }
	  
    //
    // if flag set, and read already there, skip that read
    //
    if (uniquereads == 1) {
      if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
	  ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
	FREE(p);
	FREE(line);
	FREE(myline);
	continue;
      } else {
	if (st == 1)
	  set_entry_in_binarized_array(a_unqreads_F, pos);
	else
	  set_entry_in_binarized_array(a_unqreads_R, pos);
      }
    }
	  
    //if (verbose == 1)
    //	printf("unique read starting pos= %d\n", pos);
	  
    // printf("mmstr= %s\n", mmstr);
	  
    // actual read sequence
    readseq    = strdup(p[9]);
    readseqlen = strlen(readseq);
	  
    // parse CIGAR
    cigr       = strdup(p[5]);

    for (i=0; i<readseqlen; i++)
      myreadpos[i] = (char)i;
	  
    idxnewread = 0;
    idxread    = 0;
	  
    k          = 0; // num digits
    cigrlen    = strlen(cigr);
	  
    // DEBUG
    /***
	todebug= 0;
	if ( abs( pos - bugposition ) < 5)
	todebug= 0;
	      
	//	if (pos == bugposition)
	      
	if (todebug == 1){
	verbose= 1;
	printf("hi bug here, cigar= %s\n", cigr);
	}else{
	verbose= 0;
	}
    ***/
	  
    int idxmmpos_read= 0;		// match/mismatch position on the read	-> 
    for (i=0; i<cigrlen; i++) {
	    
      if (cigr[i] == 'M') {
	lb = atoi(snum);
	      
	// DEBUG: now use hash based alignedpos and alignedseq
	// aln_nt:
	//	-> pos=       position on the chromosome relative to the start of the mapping
	//	-> nt=        nucleotide in the read
	//	-> myreadpos= position in the read (including all match/mismatch/insert bases)
	// tmpread:
	//	-> pos=       n-th match/mismatch base
	//	-> nt         corresponding nucleotide
	//	-> 
	      
	int tmpidx;
	for (tmpidx=0; tmpidx < lb; tmpidx++){
	  if ( (aln_nt = (my_aln*)malloc(sizeof(my_aln))) == NULL){
	    fprintf(stderr, "no memory left!\n");
	    exit(-1);
	  }
	  aln_nt->tpos=     idxnewread + tmpidx;
	  aln_nt->nt=      readseq[idxread + tmpidx ] ;
	  aln_nt->myreadpos= myreadpos[idxread + tmpidx ] ;
	  //if (verbose)
	  //	printf("cigar loop: adding %c to relative chr pos %d, readpos= %d\n", aln_nt->nt, aln_nt->tpos, aln_nt->myreadpos);
	  HASH_ADD_INT(aln_seq, tpos, aln_nt);
		
	  if ( (tmpread = (my_aln*)malloc(sizeof(my_aln))) == NULL){
	    fprintf(stderr, "no memory left!\n");
	    exit(-1);
	  }
		
	  // read position (i.e. n-th match/mismatch base) -> nt
		
	  // tmpread->tpos= myreadpos[idxread + tmpidx ] ;
	  tmpread->tpos= idxmmpos_read;
	  tmpread->nt= aln_nt->nt ;
	  tmpread->chrpos= idxnewread + tmpidx ;	// position in chromosome 
		
							//if (verbose)
							//	printf("readpos2nt: adding %c to read mm position %d (idxmmpos_read= %d), chrpos= %d\n", 
							// tmpread->nt, tmpread->tpos, idxmmpos_read, tmpread->chrpos);
	  HASH_ADD_INT(readpos2nt, tpos, tmpread);
	  idxmmpos_read++;
	}
	      
	/*
	  memcpy(alignedseq + idxnewread,   readseq + idxread, lb);
	  memcpy(alignedpos + idxnewread, myreadpos + idxread, lb);
	*/
	      
	// update both idx
	idxnewread += lb;
	idxread    += lb;
	      
	k = 0;  // reset digits
	      
      } else if ( cigr[i] == 'N' ){
	lb= atoi(snum);
	      
	idxnewread+= lb;
	k= 0;
	      
      } else if (cigr[i] == 'D') {
	// deleting the next lb bases from reference genome
	lb = atoi(snum);
	      
	my_aln *s;
	      
	// add gaps
	// printf("DEPRECATED!\n");
	for (l=idxnewread; l<idxnewread+lb; l++) {
	  HASH_FIND_INT( aln_seq, &l, s );  /* s: output pointer */
	  if (s != NULL){
	    fprintf(stderr, "WARNING! found pos= %d that's not supposed to be deleted\n", l);
	    // s->nt= '-';
	    // s->myreadpos= -1;
	  }
		
	  /***
	   // not adding this deletion to aln_seq 
	   //	-> if seg fault then it means we erroneously entered a deleted position
	   printf("warning: adding chr position %d (deletion) to aln_seq \n", l);
		 
	   if ( (aln_nt = (my_aln*)malloc(sizeof(my_aln))) == NULL){
	   fprintf(stderr, "no memory left!\n");
	   exit(-1);
	   }
	   aln_nt->tpos=     l;
	   aln_nt->nt=      '-';
	   aln_nt->myreadpos= -1;
		 
	   // DEBUG:
	   // 	-> do we add this position?
	   HASH_ADD_INT(aln_seq, tpos, aln_nt);
		 
	  ****/
		
		
	  /* 
	     alignedseq[l] = '-';
	     alignedpos[l] = -1;
	  */
	}
	      
	// update index in new read
	idxnewread += lb;
	      
	k = 0; // reset digits
	      
      } else if (cigr[i] == 'I') {
	lb = atoi(snum);
	      
	// update idx in read
	idxread += lb;
	      
	k = 0; // reset digits
	      
      } else if (cigr[i] == 'S') {
	lb = atoi(snum);
	k = 0; // reset digits
	      
      } else {
	// num
	snum[k] = cigr[i]; // add digits
	k++;
	snum[k] = '\0'; // add end char
      }
    }
	  
    if (0){
      printf("DEPRECATED\n");
      alignedseq[idxnewread] = '\0';
    }
	  
    if (st == -1) {
      my_aln *s;
      for(s= aln_seq; s != NULL; s=s->hh.next) {
	//if (verbose)
	//	printf("st= -1 => reversing indices: position %d-> nt %c-> readpos %d\n", s->tpos, s->nt, s->myreadpos);
	if (s->myreadpos != -1)
	  // turn around others
	  s->myreadpos = (char)( readseqlen - s->myreadpos - 1 );
      }
	    
      /*
	for (i=0; i<strlen(alignedseq); i++) {
	// leaves -1 as they are
	if (alignedpos[i] != -1)
	// turn around others
	alignedpos[i] = (char)( readseqlen - alignedpos[i] - 1 );
	}
      */
    }
	  
    // printf("hi there cigar= %s\n", cigr);
	  
    alignedseqlen = HASH_COUNT(aln_seq);
	  
    // printf("alignedseqlen= %d\n", alignedseqlen);
	  
    // alignedseqlen = strlen(alignedseq);
	  
    // count aligned nucleotides
    rst     = pos + truncate;
    ren     = min(pos + alignedseqlen - truncate, chrlen );
    rennotr = min(pos + alignedseqlen, chrlen);

    //if (verbose)
    //	printf("pos= %d, , mmstr= %s, cigar= %s, rst= %d, ren= %d, rennotr= %d, truncate= %d\n", pos, mmstr, cigr, rst, ren, rennotr, truncate);
	  
    int numalignednts = 0;
    int posinread     = 0;
	  
    // DEBUG: continue here
    // go through all keys of aln_seq
    // i= pos+ truncate + key
    my_aln *s;
	  
    // working on readpos 
    // go through all the reads and register on the chromosome where the reads come from (i.e. n-th base in a read)
    for(s= aln_seq; s != NULL; s=s->hh.next) {
      //if (verbose)
      //	printf("traverse: position %d -> nt %c -> readpos %d\n", s->tpos, s->nt, s->myreadpos);

      // assuming...
      k= s->myreadpos;
      co[k]= 0;	// what is co??
      // printf("setting co[%d]= 0\n", k);
      i= s->tpos + pos;		// position on chromosome
	    
      //if (verbose)
      //	printf("trying to add 'pos= %d' to readpos %d [ %d ], nt= %c, readcounts= %d, current readpos[%d][]= %p\n", 
      //		posinread, i-1, (*readcounts)[i-1], s->nt, (*readcounts)[i-1], i-1, (*readpos)[i-1]  );
	    
      // if ((i >= rst) && (i < ren) && (s->nt != '-') && ((*readcounts)[i-1] < USHRT_MAX)) 
	    
      // removed boundary checking 
      if ( (s->nt != '-') && ((*readcounts)[i-1] < USHRT_MAX)) {
	      
	// if space allocated already, register read positions (so we can use position-specific err rate later)
	if ((*readpos)[i-1] != 0) {
		
	  // calculate where in read we are
	  // posinread = alignedpos[k];
	  posinread = s->myreadpos;
		
	  // add pos
	  //if (verbose)
	  //	printf("\tadding 'posinread= %d' to readpos %d [ %d th ]\n", posinread, i-1, (*readcounts)[i-1]);
	  (*readpos)[i-1][  (*readcounts)[i-1] ] = (char)posinread;
		
	}
	      
	// add read
	(*readcounts)[i-1] ++;
	numalignednts ++;
	      
      }
	    
    }
	  
	  
    /*
      for (i=pos, k=0; i<rennotr; i++,k++) {
      co[k] = 0;
	    
      if ((i >= rst) && (i < ren) && (alignedseq[k] != '-') && ((*readcounts)[i-1] < USHRT_MAX)){ 
	    
      if ((*readpos)[i-1] != 0) {
	    
      // calculate where in read we are
      posinread = alignedpos[k];
	    
      // add pos
      (*readpos)[i-1][  (*readcounts)[i-1] ] = (char)posinread;
	    
      }
	    
      // add read
      (*readcounts)[i-1] ++;
      numalignednts ++;
	    
      }
      }
    */
	  
    (*numnt) += numalignednts;
	  
    // mmstr -> to uppercase
    for (i= 0; mmstr[i]; i++)
      mmstr[i]= toupper(mmstr[i]);
	  
    //
    // parse MUTATIONS
    //
    i          = 0; // iterate thru mmstr
    k          = 0; // num digits
    idxalnread = 0; // index in aligned read
    insertion  = 0; // are we in an insertion
    mmstrlen   = strlen(mmstr);
	  
    // DEBUG: beware of lowercase mmstr !!!
    // printf("entering MMSTR loop, MMSTR= %s, pos= %d\n", mmstr, pos);

    while (i<mmstrlen) {
      // printf("mmstr [%d] = %c\n", i, mmstr[i]);
	    
      if ((mmstr[i] == 'A') ||
	  (mmstr[i] == 'C') ||
	  (mmstr[i] == 'G') ||
	  (mmstr[i] == 'T')) {

	// jump by snum
	//	note that we are advancing according to the MD tag
	//		-> i.e. probabily the MD tag itself does not tell you the absolute chromosome position !!

	if (k > 0) {
	  snum[k] = '\0';
	  inum = atoi(snum);
	  idxalnread += inum;
	  k = 0; // reset
	  // printf("advancing by %d bases; idxalnread= %d\n", inum, idxalnread);
	}

	// if we are not in an insertion (or a deletion in the read)
	if (insertion == 0) {

	  // make sure the readpos is not outside sequence
	  if (pos + idxalnread <= chrlen) {

	    /*
	      if (verbose == 1)
	      printf("Variation %c (REF=%c) -> %c at pos = %d (chrlen=%d)\n",
	      mmstr[i], toupper(seq[ pos + idxalnread - 1]),
	      alignedseq[ idxalnread ], pos + idxalnread, chrlen);
	    */	       

	    // increase count for mismatch at that position in the read
	    // mismatch vector
	    // idxalnread should be <= alignedseqlen

	    co[ idxalnread ] ++;
	    // printf("co[%d]= %d\n", idxalnread, co[idxalnread]);

	    // if mutation in central part, then l
	    //	-> central part means: idxalnread is >= truncate from either end of the read
	    //	-> truncate: by default = 0; can be set in command line
	    if ((co[idxalnread] == 1)    &&
		(idxalnread >= truncate) &&
		(idxalnread < alignedseqlen-truncate)) {

	      // DEBUG: important! 
	      //	want to find out which nucleotide the base is mutated to...
	      //	should match positions after taking into acocunt insertions !!!

	      HASH_FIND_INT( readpos2nt, &idxalnread, s );  // s: output pointer 
	      // printf("looking for readpos2nt [ %d ] in readpos2nt\n", idxalnread);
	      // printf("\t-> position= %d, nt= %c, chr pos= %d\n", s->tpos, s->nt, s->chrpos);

	      // printf("looking for aln_seq [ %d ] \n", idxalnread);
	      // HASH_FIND_INT( aln_seq, &idxalnread, s );  // s: output pointer 

	      // if (alignedpos[ idxalnread ] == -1)
	      // 	die("We have a problem, position should not be -1 (in deletion)\n");

	      // printf("WARNING: idxalnread should be mapped to the chr position!!!\n");

	      // use s->chrpos (mm position -> chr position) since there may be gaps inside a read on the genome
	      // absposmut = pos + idxalnread - 1;
	      absposmut = pos + s->chrpos - 1;

	      // this is the refgen nt
	      // nt        = toupper( seq[ absposmut ] );

	      // nt        = alignedseq[ idxalnread ];

	      // get A/C/G/T index for this mutation
	      //	(mutation called based on MD tag, stored in readpos2nt)
	      nt        = s->nt;
	      ntidx     = data_nucl(nt);
	      // printf("we got %c at pos %d\n", nt, idxalnread);

	      // not an N or IUPAC code
	      if (ntidx >= 0) {
		// mutation count per position
		(*mutacounts)[ absposmut ] ++;

		if (st == 1)
		  (*mutacounts_F)[ absposmut ] ++;

		if ((*ntcounts) != 0) {
		  // allocate memory for mutated nucleotide
		  if ( (*ntcounts)[absposmut] == 0) {
		    (*ntcounts)[absposmut] = (unsigned char*)calloc(4, sizeof(unsigned char));
		    if ( (*ntcounts)[absposmut] == 0 )
		      die("problem allocaing *ntcounts[]\n");
		  }

		  // register the mutation on chromosome
		  // note: ucsc genome browser coordinate is 1 based
		  if ( (*ntcounts)[absposmut][ ntidx ] < UCHAR_MAX){
		    (*ntcounts)[absposmut][ ntidx ] ++;

		    // if (verbose == 1)
		    // 	printf("adding one to ntcounts [%d][%d] (%c)\n", absposmut, ntidx, nt);

		  }	
		}
		(*nummm) ++;

	      } // if (ntidx

	    } // if mutation in central part

	  } // if mut pos within chrom range

	} // if not in insertion

	// in any case, move by one nt
	//	DEBUG: 
	//		should consider the case of deletion
	//		e.g. ^AC means deleting 2 bases from the reference genome
	idxalnread ++;
	// printf("\tadvancing by 1 nt (ref= %c), current idxalnread= %d\n", mmstr[i], idxalnread);

      } else if (mmstr[i] == '^') {
	// printf("\twe have an insertion here, k= %d\n", k);

	// jump by snum
	//	i.e. leading matching nucleotides.
	if (k > 0) {
	  snum[k] = '\0';
	  inum = atoi(snum);
	  idxalnread += inum;
	  // printf("\tadvancing idxalnread by %d, idxalnread= %d\n", inum, idxalnread);
	  k = 0; // reset

	}

	// new:
	//	look for trailing nucleotides (those deleted from the reference genome
	int tmpidx=0;
	for (tmpidx= i+1; mmstr[tmpidx]; tmpidx++){
	  // are there trailing nucleotide(s)?
	  // printf("\tmmstr[%d]= %c => numeric= %d\n", tmpidx, mmstr[tmpidx], isdigit(mmstr[tmpidx]));
	  if (isdigit(mmstr[tmpidx])){
	    // found the end of trailing, deleted nucleotides	
	    break;
	  }
	}

	// printf("\tlast index= %d, should skip i by %d\n", tmpidx, (tmpidx-i-1));
	i+= (tmpidx-i-1);

	insertion = 1; // we entered an insert

      } else if ( mmstr[i] == 'N' ){
	if (k > 0){
	  snum[k]= '\0';
	  inum = atoi(snum);
	}

	if (mmstr[i+1] == '0' ){
	  //
	  inum++;
	}

	//if (verbose == 1)
	//	printf("MMSTR: advanced by %d nt\n", inum);

	idxalnread+= inum;
	k=0;

      } else if ( isalpha(mmstr[i])){
	fprintf(stderr, "mmstr format not supported %c\n", mmstr[i]);
					
      } else {
	// must be number
	snum[k] = mmstr[i];
	k++;
	// leave insertion if we are in one
	if (insertion == 1) {
	  insertion = 0;
	}
      }
    
      i++;
    }
  
    // printf("end of MMSTR loop\n");

    // final advance
    if (k > 0) {
      snum[k] = '\0';
      inum = atoi(snum);
      idxalnread += inum;
      //if (verbose == 1)
      //	printf("MMSTR: advance by %d nt\n", inum);
      k = 0;

    }
  
    // get length of aligned sequence
     alignedseqlen = HASH_COUNT(aln_seq);
  
  
    //if (verbose)
    //	printf("[hash] alignedseqlen= %d\n", alignedseqlen);
  
    // if (idxalnread != strlen(alignedseq)) 
    if (idxalnread != alignedseqlen) {
      printf("MMSTR: Houston, we have a pb in getReadCountFromReadFile_SAM(), idxalnread= %d, len= %d, pos= %d, mmstr= %s, cigar= %s, file= %s\n", idxalnread, alignedseqlen, pos, mmstr, cigr, file);
    }
  
    FREE(p); // p[x] are just pointers
    FREE(line);
    FREE(myline);

    // clean up and reclaim memory 
    my_aln *current_base;
    while (aln_seq){
      current_base= aln_seq;
      HASH_DEL(aln_seq, current_base);
      FREE(current_base);
    }
    while (readpos2nt){
      current_base= readpos2nt;
      HASH_DEL(readpos2nt, current_base);
      FREE(current_base);
    }

    FREE(cigr);
    FREE(readseq);

    /**
       if ( todebug == 1){
       // die("let's stop and check!\n");
       }
    **/
					 			
    (*numreads) ++;
				
    if ((verbose == 1) && ( ((*numreads) % 100000) == 0)) {
      printf("Read %d reads from %s\n", *numreads, file);
    }
  }

  FREE(co);

  if (verbose == 1) 
    printf("Number of reads = %d\n", *numreads);
	
  LineIclose(&li);

  if (uniquereads == 1) {
    FREE(a_unqreads_F);
    FREE(a_unqreads_R);
  }

  if (verbose == 1) 
    printf("About to exit after reading %s\n", file);

  //fclose(fp);
  FREE(snum);
  FREE(alignedseq);
  FREE(myreadpos);
  FREE(alignedpos);

  if (verbose == 1) 
    printf("About to exit after reading %s (2)\n", file);
  
}

// goal: count mismatches
// DEBUG: working on extended CIGAR format
// to return:
// 	poscounts
//	numnt
//	nummm
void updateMismatchCountAtReadPos_SAM(const char* file, const int chrlen, int** poscounts, long* numnt, long* nummm, const int uniquereads, const int readlen)
{

  FILE* fp = 0;
  int   bufflen = 10000000;
  char* buff;
	int   cnt;

	char* s;
	char** p;
	int   m;
	char* prevs = 0;
	int   i =0;
	int   k = 0;
	char* mys;
	int   j = 0, l;
	int   pos = 0;
	int   st = 0;
	unsigned char* a_unqreads_F = 0;
	unsigned char* a_unqreads_R = 0;
	int   verbose = 0;
	char* snum = (char*)calloc(100, sizeof(char));
	int   bitflag;
	char* mmstr;
	int   mmstrlen = 0;
	char* readseq;
	int   readseqlen;
	int   lb;
	char* cigr;
	char* alignedseq;
	char* readpos;
	char* alignedpos;
	int   idxnewread = 0;
	int   idxread    = 0;
	int   cigrlen = 0;
	int   inum = 0;
	int   idxalnread = 0;
	int   insertion = 0;
	int   idxmd = -1;

	// DEBUG:
	//	-> we just want to know the number of mismatches along the read length
	alignedseq = (char*)calloc(readlen * 2, sizeof(char));
	readpos    = (char*)calloc(readlen * 2, sizeof(char));
	alignedpos = (char*)calloc(readlen * 2, sizeof(char));


	fp = fopen(file, "r");
	if (fp == 0) {
		printf("Cannot open read file %s.\n", file);
		exit(1);
	}

	if (uniquereads == 1) {
		a_unqreads_F = create_binarized_array(chrlen);
		a_unqreads_R = create_binarized_array(chrlen);
	}



	j = 0;
	int cline = 0;

	// printf("we are in updateMismatchCountAtReadPos_SAM()\n");

	// DEBUG
	char * buff_st;

	while (1){	// !feof(fp)) 

		// allocate data for new chunk
		buff   = (char*)calloc(bufflen, sizeof(char));
		if (!buff) {
			die("Cannot allocate buff in seqI_open_fast\n");
		}

		// read data from fp
		cnt = fread(buff, 1, bufflen, fp);
		buff_st= buff;

		if (cnt == 0) {
			FREE(buff_st);
			FREE(prevs);
			break;
		}

		// split by line, get the first line in s
		// 	strsep problem
		s = strsep(&buff, "\n");

		// if this is not the first chunk we read
		if (j != 0) {

			// alloc new s
			mys = (char*)calloc(1000, sizeof(char));

			// concatenate what is left of prev line
			strcat(mys, prevs);

			if (s != 0) {

				// concate end of same line
				strcat(mys, s);
			}

			// update s
			s = mys;

			// needed to tell next block to free mys
			cline = 1;

			// needed because the only way to have gotten here is if (s == 0) in next block
			FREE(prevs);

		}


		while (1) {

			// store current line
			prevs = s;


			// get next line ahead of time
			s = strsep(&buff, "\n");

			if (s == 0) {

				prevs = strdup(prevs);
				break;

			} else {

				// split
				// printf("working on: %s\n", prevs);
				split_line_delim(prevs, "\t", &p, &m);

				// skip if not aligned
				if (strcmp(p[5], "*") == 0) {
					FREE(p);
					//free(line);
					continue;
				}

				// get mm string
				// DEBUG: m is varying. so should check for the tag MD:Z:*
				//if (idxmd == -1) {
				idxmd      = -1;				
				int tmpidx = m-1;
				while (tmpidx >= 0) {
				  if (strstr(p[tmpidx], "MD:Z") != 0) {
				    idxmd = tmpidx;
				    break;
				  }
				  tmpidx --;
				}
				//}

				if (idxmd == -1) {
				  die("Cannot find MD:Z tag, check format.\n");
				}

				if (strstr(p[idxmd], "MD:Z")== NULL){
					//printf("ERROR! no MD tag found!!!\n");
					printf("ERROR! no MD tag found!!! (last tag is %s)\n", p[idxmd]);
      
				}

				mmstr = p[idxmd] + 5;

				

				//if (verbose == 1)
				//printf("line=%s\n", line);
				//free(line);
				//


				if (verbose == 1)
					printf("TOTOTTOO=%s\n", p[5]);


				//
				// parse SAM format
				//
				if (verbose == 1)
					printf("MMSTR = %s (%s)\n", mmstr, p[idxmd]);


				// read pos (starts at 1)
				pos = atoi(p[3]);

				if (verbose == 1)
					printf("pos = %d\n", pos);

				// strand
				//	+1 => minus strand
				//	-1 => plus strand
				bitflag = atoi(p[1]);

				// unmapped
				if ( bitflag & ( 1 << 2 ) ){
					//fprintf(stderr, "WARNING: skipping unmapped read: %s %s %s\n", p[0], p[1], p[2]);
					FREE(p);
					continue;
				}

				// check bits 4 anb 5 (from 0)
				if ((bitflag & (1 << 4)) ||
				        (bitflag & (1 << 5))) {

					//if (bitflag == 16)
					st = -1;
					if (verbose == 1)
						printf("STRAND = +++++++++++\n");
				} else {
					st = 1;
					if (verbose == 1)
						printf("STRAND = -----------\n");

				}

				//
				// if flag set, and read already there, skip that read
				//
				if (uniquereads == 1) {
					if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
					        ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
						FREE(p);
						continue;
					} else {
						if (st == 1)
							set_entry_in_binarized_array(a_unqreads_F, pos);
						else
							set_entry_in_binarized_array(a_unqreads_R, pos);
					}
				}

				// actual read sequence
				readseq    = strdup(p[9]);
				readseqlen = strlen(readseq);

				// parse CIGAR
				cigr       = strdup(p[5]);

				if (verbose == 1) {
					printf("ALLOC: %ld\n", (long)(alignedseq));
					printf("ALLOC: %ld\n", (long)(readpos));
					printf("ALLOC: %ld\n", (long)(alignedpos));
				}

				//if (readseqlen > readlen)
				//  printf("ERROR! readseqlen = %d > readlen! = %d!!\n", readseqlen, readlen);

				// set up index along the read length
				for (i=0; i<readseqlen; i++)
					readpos[i] = (char)i;
				
				idxnewread = 0;
				idxread    = 0;

				k          = 0; // num digits
				cigrlen    = strlen(cigr);

				// parse cigar tag
				for (i=0; i<cigrlen; i++) {
					if (verbose == 1)
						printf("cigr[i=%d] = %c\n", i, cigr[i]);

					if (cigr[i] == 'M') {
						lb = atoi(snum);

						// actually don't care about the absolute position of the read in the genome, just care about the position in the read
						memcpy(alignedseq + idxnewread, readseq + idxread, lb);
						memcpy(alignedpos + idxnewread, readpos + idxread, lb);

						// update both idx
						idxnewread += lb;
						idxread    += lb;

						k = 0;  // reset digits
						if (verbose == 1)
							printf("M, lb=%d\n", lb);
					} else if (cigr[i] == 'N') {
						// long gaps 

						lb= atoi(snum);
						
						// we don't care about the absolute genomic position here, so don't advance the index when seeing a gap
						// idxnewread += lb;
						
						k= 0;

					} else if (cigr[i] == 'D') {
						// deletion from genome
						//	-> does not count as an error 
						
						lb = atoi(snum);

						// add gaps
						for (l=idxnewread; l<idxnewread+lb; l++) {
							alignedseq[l] = '-';
							alignedpos[l] = -1;
						}

						// update index in new read
						idxnewread += lb;

						k = 0; // reset digits
						if (verbose == 1)
							printf("D, lb=%d\n", lb);
					} else if (cigr[i] == 'I') {
						// insertion to genome
						lb = atoi(snum);

						// update idx in read
						// DEBUG: removed since an insertion should not count as an error
						// idxread += lb;

						k = 0; // reset digits
						if (verbose == 1)
							printf("I, lb=%d\n", lb);

					} else if (cigr[i] == 'S') {
						lb = atoi(snum);
						k = 0; // reset digits
						if (verbose == 1)
							printf("S, lb=%d\n", lb);
					} else if (isalpha(cigr[i])) {
						// DEBUG
						fprintf(stderr, "WARNING: unsupported CIGAR string (%s)\n", cigr);
						k=0;

					} else {
						// num
						snum[k] = cigr[i]; // add digits
						k++;
						snum[k] = '\0'; // add end char
					}
				} // end for cigarlen

				alignedseq[idxnewread] = '\0';

				if (verbose == 1)
					printf("OLD = %s, len= %d\n", readseq, (int)strlen(readseq));
				if (verbose == 1)
					printf("NEW = %s, len= %d\n", alignedseq, (int) strlen(readseq));

				if (verbose == 1)
					printf("strand= %d\n", st);

				if (st == -1) {
					if (verbose == 1)
						printf("reversing indices\n");

					for (i=0; i<strlen(alignedseq); i++) {
						// leaves -1 as they are
						if (alignedpos[i] != -1)
							// turn around others
							alignedpos[i] = (char)( readseqlen - alignedpos[i] - 1 );
					}
				}

				if (verbose == 1) {
					// print read position indices
					for (i=0; i<strlen(alignedseq); i++) {
						printf(" %d", (int)(alignedpos[i]));
					}
					printf("\n");
				}
				
				// DEBUG
				//	make mmstr uppercase??
				// 	note: find out why mmstr is sometimes in lowercase
				
				
				// printf("OLD: mmstr= %s %p\n", mmstr, mmstr);	
			
				for (i=0; mmstr[i]; i++)
					mmstr[i]= toupper(mmstr[i]);

				// printf("NEW: mmstr= %s %p\n", mmstr, mmstr);	

				i = 0; // iterate thru mmstr
				k = 0; // num digits

				idxalnread = 0; // index in aligned read
				insertion  = 0; // are we in an insertion

				mmstrlen = strlen(mmstr);


				//
				// MD field (not necessarily the last field)
				// example MMSTR= 10A5^AC6
				//	-> 10 matches + A on reference + 5 matches + 2 deletions (AC) from referene + 6 matches
				//
				// note: sum of MD string is always <= read length
				while (i<mmstrlen) {
					if ((mmstr[i] == 'A') ||
					        (mmstr[i] == 'C') ||
					        (mmstr[i] == 'G') ||
					        (mmstr[i] == 'T')) {

						// jump by snum
						// 	i.e. inum bases before this A/C/G/T
						//	e.g. 75=> 75 matches
						if (k > 0) {
							snum[k] = '\0';
							inum = atoi(snum);
							idxalnread += inum;
							if (verbose == 1)
								printf("MMSTR: advance by %d nt\n", inum);
							k = 0; // reset

							// increase numnt
							(*numnt) += inum;

						}

						if (verbose == 1)
							printf("MMSTR: %c\n", mmstr[i]);

						// if we are not in an insertion (or a deletion in the read)
						if (insertion == 0) {

							// make sure the readpos is not outside sequence
							if (pos + idxalnread <= chrlen) {

								// should store the sequences ONLY when verbose is ON !!!

								/**
								if (verbose == 1 && seq != 0){
									// only works when seq is provided (fasta input)
									printf("Variation %c (REF=%c) -> %c at pos = %d (chrlen=%d)\n",
									       mmstr[i], toupper(seq[ pos + idxalnread - 1]),
									       alignedseq[ idxalnread ], pos + idxalnread, chrlen);
									       }
								**/

								// increase count for mismatch at that position in the read
								if (alignedpos[ idxalnread ] == -1)
									die("We have a problem, position should not be -1\n");

								// DEBUG: the only thing that needs changing
								(*poscounts)[ (int)( alignedpos[ idxalnread ] ) ] ++;
								if (verbose == 1)
									printf("adding one count to position %d, starting pos= %d\n", (int)alignedpos[ idxalnread ], pos );

								(*nummm) ++;                   // increase overall mm count

								// increase numnt
								(*numnt) ++;

								// DEBUG:
								// 	currently turned off to allow setting seq len to a large number 
								//	for cases where the user only has the alignment output and no original fasta sequence..
								//	(will figure out a better way to do this.. 
								//		all that's needed is the seq. length and not the sequence itself) 

								/**
								if ((seq != 0) && (mmstr[i] !=  toupper(seq[ pos + idxalnread - 1]))) {
									fprintf(stderr, "Attention:");
									fprintf(stderr, "Variation %c (REF=%c) -> %c at pos = %d (chrlen=%d)\n",
									        mmstr[i], toupper(seq[ pos + idxalnread - 1]),
									        alignedseq[ idxalnread ], pos + idxalnread, chrlen);
									//exit(1);
								}
								**/

							}

						}

						// move by one nt
						idxalnread ++;
						if (verbose == 1)
							printf("MMSTR: advance by 1\n");

					} else if (mmstr[i] == '^') {
						// insertion mode (== deletion from genome)

						// jump by snum
						// 	previous number means "number of matches"
						//	trailing letters are the letters that should be deleted from the reference
						if (k > 0) {
							snum[k] = '\0';
							inum = atoi(snum);
							idxalnread += inum;
							if (verbose == 1)
								printf("MMSTR: advance by %d nt\n", inum);
							k = 0; // reset

							// increase numnt
							(*numnt) += inum;

						}

						insertion = 1; // we entered an insert
						if (verbose == 1)
							printf("MMSTR: entered insert\n");
					} else if (mmstr[i] == 'N'){
						// no sequence determined on the reference genome
						// printf("k= %d\n", k);
						if (k > 0) {
							snum[k] = '\0';
							inum = atoi(snum);
							// inum is usually 0
						}

		
						if (mmstr[i+1] == '0'){
							// N followed by 0 => advance one more
							inum++;
						}

						if (verbose == 1)
							printf("MMSTR: advance by %d nt\n", inum);
	
						idxalnread+= inum;
						k= 0;
					
					} else if (isalpha(mmstr[i])){
						fprintf(stderr, "WARNING! MMSTR format not recognized!! mmstr=%s, pos= %d, cigr= %s, file= %s\n", mmstr, pos, cigr, file);
					
					} else {
						// must be number
						snum[k] = mmstr[i];
						k++;
						// leave insertion if we are in one
						if (insertion == 1) {
							if (verbose == 1)
								printf("MMSTR: left insert\n");
							insertion = 0;
						}
					}

					i++;
				}
				// if last letter is number => advance by 1
				// if (mmstr[i-1]== '0')
				//	printf("WARNING! MMSTR format not recognized!! mmstr=%s, pos= %d, cigr= %s, file= %s\n", mmstr, pos, cigr, file);
				
				// final advance
				if (k > 0) {
					snum[k] = '\0';
					inum = atoi(snum);

					/*
					if (mmstr[i-2]== 'N' && inum== 0){
						inum= 1;
						printf("fixing the last N\n");
					}
					*/
					idxalnread += inum;
					if (verbose == 1)
						printf("MMSTR: final advance by %d nt\n", inum);
					k = 0;
					// increase numnt
					(*numnt) += inum;

				}


				// check if last index equals read length 
				if (idxalnread != strlen(alignedseq)) {
					// if (idxalnread != strlen(readseq)) 
					printf("MMSTR: Houston, we have a pb in updateMismatchCountAtReadPos_SAM(), pos= %d, cigr= %s, mmstr= %s, idxalnread= %d, len of alignedseq= %d\n", 
							pos, cigr, mmstr, idxalnread, (int)strlen(alignedseq));
				}

				if (verbose == 1) {
					printf("FREE: %ld\n", (long)(alignedseq));
					printf("FREE: %ld\n", (long)(readpos));
					printf("FREE: %ld\n", (long)(alignedpos));
				}

				FREE(readseq);
				FREE(cigr);

				FREE(p); // p[x] are just pointers

				if (verbose == 1)
					printf("\n");

			} // else (line exists)
		} // while (1)

		if (j != 0)
			FREE(mys);

		FREE(buff_st);
		j++;
	}

	if (uniquereads == 1) {
		FREE(a_unqreads_F);
		FREE(a_unqreads_R);
	}

	fclose(fp);
	FREE(snum);
	FREE(alignedseq);
	FREE(readpos);
	FREE(alignedpos);

}



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
  else if (strcmp(formattext, "bamunq") == 0)
    format = BAMUNQ;
  else if (strcmp(formattext, "bam") == 0)
    format = BAM;
  else {
    die("Format unrecognized\n");
  }
  return format;
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



int CountReadsInChrInterval(char* file, int format, char* chr, int i1, int i2)
{
  
  LineI li;
  char* l;
  char** p;
  int   m;
  int   numreads = 0;
  //int   pos = 0;
  
  li.bufflen = 1000000;
  LineIopen(&li, file);
  while ((l = nextLine(&li))) {

    if (format == 3) {
      chomp(l);
      split_line_delim(l, "\t", &p, &m);

      if ((m >= 12) && (strcmp(p[11], "XT:A:U") == 0) && (strcmp(p[2], chr) == 0) && (atoi(p[3]) >= i1) && (atoi(p[3]) <= i2)) {	
	numreads ++;
      }

      free(p);
    } else {
      die("Unsupported format.\n");
    }
    free(l);
  }

  LineIclose(&li);
  
  return numreads;
}


//
// count !!! unique mappers !!!
//
int CountReads(char* file, int format)
{
  
  //LineI li;
  //char* l;
  // char** p;
  //int   m;
  int   numreads = 0;

  MappedRead     r;
  ReadI          ri;
    
  if (!ReadIopen(&ri, file, format)) {
    printf("Cannot open file %s\n", file);
  }
  
  while (nextRead(&ri, &r)) {
    
    if (r.uniqmap != 0) {
      numreads++;
    } 

  }

  ReadIclose(&ri);

  
  return numreads;
}




/*
  SN
 FUNCTION: tally up read counts at all genomic positios in a vector
 file = read file
 format = 0=eland, 1=mit, 2=bed, 3=sam 
 char = useless 
 chrlen = chrom len
 fraglen = actual fragment length (needed for extension); must be > 0 for extension to occur; BUT if fraglen eq -10, then don't even count the whole read but pos
 uniqmap
 counts = counts (output)
 numreads = num reads (output) 
 */ 
void getCountFromReadFile(char* file, int format, char* chr, long chrlen, int readlen, int fraglen, int uniquereads,  int uniqmap, unsigned short** counts, int* numreads, int* numclonalreads)
{
	  
  char** p;
  //int   m;
  int   i		 	= 0;
  int   pos		        = 0;
  int   lenread	                = 0;
  int            st		= 0;
  //int            bitflag        = 0;
  unsigned char* a_unqreads_R   = 0;
  unsigned char* a_unqreads_F   = 0;
  //char*          cigr           = 0;
  //LineI          li;
  //char*          l;
  //int*           blocksizes;
  //int*           blockoffsets;
  //int            numblocks;
  //int            mapq;
  MappedRead     r;
  ReadI          ri;
  int            map;
 
  
  
  if (uniquereads == 1) {
    a_unqreads_R = create_binarized_array(chrlen);
    a_unqreads_F = create_binarized_array(chrlen);
  }
  
  *numreads		= 0;
  *numclonalreads       = 0;

  //li.bufflen = 10000000;
  //LineIopen(&li, file);
  //while ((l = nextLine(&li))) {
  //chomp(l);
  //split_line_delim(l, "\t", &p, &m);
  
  
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
      } else if (fraglen == STARTONLY) {
	cnt_en = cnt_st+1;
      } else {	       // no extension required, just use read length
	cnt_en = pos + lenread;
      }
      
    } else if (st == -1) {
      
	cnt_en = pos+lenread; // count will start at beginning of read
	
	if (fraglen > 0) { // if extension required				      
	  cnt_st = pos - (fraglen - lenread);
	} else if (fraglen == STARTONLY) {
	  cnt_st = cnt_en-1;
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
  
  //free(p); // p[x] are just pointers
  //free(l);
  
  //LineIclose(&li);
  ReadIclose(&ri);
  
  if (uniquereads == 1) {
    free(a_unqreads_F);
    free(a_unqreads_R);
  }
  
}




/*
  FUNCTION: tally up read counts at all genomic positios in a vector
  file = read file
  format = 0=eland, 1=mit, 2=bed, 3=sam 
  char = useless
  chrlen = chrom len
  fraglen = actual fragment length (needed for extension)
  counts = counts (output)
  numreads = num reads (output)

*/
void getCountFromReadFileOLDBEFCLONAL(char* file, int format, char* chr, long chrlen, int readlen, int fraglen, int uniquereads, int uniqmap, unsigned short** counts, int* numreads)
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
	int   bitflag = 0;
	unsigned char* a_unqreads_R = 0;
	unsigned char* a_unqreads_F = 0;
	
	
	fp = fopen(file, "r");
	if (fp == 0) {
		printf("Cannot open read file %s.\n", file);
		exit(1);
	}
	
	if (uniquereads == 1) {
		a_unqreads_R = create_binarized_array(chrlen);
		a_unqreads_F = create_binarized_array(chrlen);
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
					
					// SAM
				} else if (format == 3) {
					
				  // must be sure that XT:A:U is always here though !!! (even for SOLID)
				  if ( (m >= 12) && ( (uniqmap == 0) || ((uniqmap == 1) && (strcmp(p[11], "XT:A:U") == 0)) ) ) {	    
				  		// TO DO:
						// 	put in code to get ride of unmapped reads if needed. 
						bitflag = atoi(p[1]);
						pos     = atoi(p[3]);
						st      = SAMbitflag_strand(bitflag);
						lenread = strlen(p[9]);
						process = 1;
					}
					
				}
				
				
				if (process == 1) {
					
					// eliminate clonal read if required
					if (uniquereads == 1) {
						if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F, pos) == 1)) ||
							((st == -1) && (get_entry_in_binarized_array(a_unqreads_R, pos) == 1))) {
							free(p);
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

					// special case, where we just want to count the number of reads mapping at beginning
					if (fraglen == -10) {
					  cnt_en = cnt_st+1;
					}

					// update counts
					for (i=cnt_st; i<cnt_en; i++) {
						if ((*counts)[i] < 65535) 
							(*counts)[i]++;					
						else 
						  fprintf(stderr, "Warning: max count (short) reached\n"); 
					}
					
					(*numreads) ++;
				} // end if (process == 1)
				
				free(p); // p[x] are just pointers
			}
		}
		
		//printf("s=%s\n", prevs);
		
		//getchar();
 		
		free(tmpbuff);
		j++;
	}
	
	if (uniquereads == 1) {
		free(a_unqreads_F);
		free(a_unqreads_R);
	}
	
}

 
void getCountFromReadFileOLD(char* file, int format, char* chr, int readlen, int fraglen, int** counts, int* numreads)
{

  FILE* fp = 0;
  int   bufflen = 10000000;
  char* buff;
  int   cnt;
  //char* ptr;
  char* s;
  char** p;
  int   m;
  char* prevs = 0;
  int   i =0;
  char* mys;
  int   j = 0;
  //int   cntc = 0;
  int   pos = 0;
  int   lenread = 0;
  int   st = 0;

  fp = fopen(file, "r");
  if (fp == 0) {
    die("Cannot open read file.\n");
  }
  
  *numreads = 0;
  
  j = 0;
  while (!feof(fp)) {

    // allocate buffer for new chunk
    buff   = (char*)calloc(bufflen, sizeof(char)); 
    if (!buff) {
      die("Cannot allocate buff in seqI_open_fast\n");
    }    
    cnt = fread(buff, 1, bufflen, fp);    
    
    i = 0;

    s = strsep(&buff, "\n");
    if (j != 0) { 
      mys = (char*)calloc(1000, sizeof(char));
      //printf("Adding %s\n", prevs);
      strcat(mys, prevs);
      //printf("Adding %s\n", s);
      strcat(mys, s);
	s = mys;
    }

    while (1) {      
      
      //split_line_delim(s, "\t", &p, &m);	
      
      prevs = s;
      i++;

      s = strsep(&buff, "\n");
      if (s == 0)
	break;
      else {
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
	  if (st == 1) {

	    for (i=pos; i<pos+fraglen; i++) {
	      (*counts)[i]++;
	    }

	  } else if (st == -1) {
	    	    
	    for (i=pos+lenread; i>pos-fraglen-lenread; i--) {
	      (*counts)[i]++;
	    }

	  } else {
	    printf("Problem interpreting strand (%s)\n", p[8]);
	    exit(0);
	  }
	  
	  (*numreads) ++;
	}
	
	free(p); // p[x] are just pointers
	
	
      }
    }

    //printf("s=%s\n", prevs);

    //getchar();
    
    free(buff);
    j++;
  }

}




int ReadIopen(ReadI* ri, char* file, int format) {

  ri->formatid = format;
  ri->verbose  = 0;

  //if (strcmp(format, "bam") == 0) {
  //  ri->formatid = BAM;
  //}
  
  if ((ri->formatid == BAM) || (ri->formatid == BAMUNQ)) {    
    ri->in = samopen(file, "rb", 0);
    if (ri->in == 0) {
      fprintf(stderr, "Fail to open BAM file %s\n", file);
      return 0;
    }
    ri->b = bam_init1();
    ri->inread = 0; // not in a read (need to get one)
    return 1;
  } else { // txt format
    //ri->li.verbose = 0;
    return LineIopen(&(ri->li), file);
  }



}

void ReadIclose(ReadI* ri) {
  if ((ri->formatid == BAM) || (ri->formatid == BAMUNQ)) {
    samclose(ri->in);
  } else {
    LineIclose(&(ri->li));
  }
      
}


int nextRead(ReadI* ri, MappedRead* r) {
  
  char* line;
  int ret;
  char* xt;
  char** p;
  int    m;
  int    bitflag;
  int    st = 0;

  if ((ri->formatid == BAM) || (ri->formatid == BAMUNQ)) {
    
    // find the next mapped read
    while (1) { 
      // read
      
      if (ri->verbose) {	
	printf("while loop iteration\n");	
      }
            
      if (ri->verbose)	  
	printf("ri->inread == 0, get new read\n");
      ret = samread(ri->in, ri->b);
      if (ret <= 0) // done
	return 0;
      
      // if read  mapped
      if (!bam1_unmapped(ri->b)) {
	if (ri->verbose)	    
	  printf("Read mapped\n");

	if (ri->formatid == BAM) {
	  xt = (char*)(bam_aux_get(ri->b, "XT"));
	  if ((xt == 0) || (xt[0] != 'A')) {
	    die("Problem interpreting XT field\n");
	  }
	  if (xt[1] == 'U')
	    r->uniqmap = 1;
	  else 
	    r->uniqmap = 0; // R or M
	} else {  // BAMUNQ
	  r->uniqmap = 1;
	}

	r->pos     = ri->b->core.pos;  // 0-based
	r->seqname = ri->in->header->target_name[ri->b->core.tid];
	st         =  bam1_strand(ri->b);
	if (st == 1)
	  r->st      = 1;
	else
	  r->st      = -1;

	r->readname = bam1_qname(ri->b);  // will be erased upon free ?
	r->lenread  = ri->b->core.l_qseq;
	//printf("XT:%s\n", (char*)(bam_aux_get(ri->b, "XT")));


	r->n_cigar =  ri->b->core.n_cigar;
	r->cigar  = bam1_cigar(ri->b);
	/*@discussion In the CIGAR array, each element is a 32-bit integer. The
  lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
  length of a CIGAR.*/
	/*
	int j;
	for (j=0; j<n_cigar; j++) {
	  //uint32_t x = cigar[j];
	  int op = cigar[j] & BAM_CIGAR_MASK; // operation
	  int l = cigar[j] >> BAM_CIGAR_SHIFT; // length
	  
	  if (op == BAM_CMATCH) {
	    printf("M");
	  }
	  if (op == BAM_CMATCH) {
	    printf("M");
	  }
	  

	}
	*/
	return 1;
      } // else ? do nothing and wait for next while (1) iter
      
    } // while (1)    

    // other format
  } else {

    while (1) { // iterate for as long as we haven't found a good read (ie a mapped one)
      
      line = nextLine(&(ri->li));
      if (line == 0) {
	return 0;
      } else {
	
	// parse line      
	chomp(line);
	split_line_delim(line, "\t", &p, &m);
	
	// format parsing
	// ELAND 
	int process = 0;
	if (ri->formatid == ELAND) { 
	  r->pos = atoi(p[7]);      
	  if (p[8][0] == 'R')
	    r->st = -1;
	  else 
	    r->st = 1;	
	  r->lenread = strlen(p[1]);
	  if (p[2][0] == 'U') {
	    r->uniqmap = 1;
	  } else {
	    r->uniqmap = 0;
	  }
	  process = 1;
	  // MIT
	} else if (ri->formatid == MIT) {
	  r->pos     = atoi(p[1]);   
	  r->lenread = atoi(p[2]) - r->pos + 1;
	  if (p[3][0] == '-')
	    r->st = -1;
	  else 
	    r->st = 1;
	  r->uniqmap = 1; // assumed
	  process = 1;
	  
	  // BED
	} else if (ri->formatid == BED) { // && (strcmp(p[3], "U0") == 0)) {      
	  r->pos     = atoi(p[1]);
	  r->lenread = atoi(p[2]) - r->pos + 1;
	  if (p[5][0] == '-')
	    r->st = -1;
	  else 
	    r->st = 1;
	  if (p[3][0] == 'U') {
	    r->uniqmap = 1;
	  } else {
	    r->uniqmap = 0;
	  }
	  process = 1;
	  
	  // SAM
	} else if (ri->formatid == SAM) {

	  if (p[0][0] == '@') {
	    process = 0; // header
	  } else {
	    
	    bitflag = atoi(p[1]);	  
	    if (SAMbitflag_ismapped(bitflag)) { // mapped
	      r->pos     = atoi(p[3]);
	      r->st      = SAMbitflag_strand(bitflag);
	      r->lenread = strlen(p[9]);
	      if (strcmp(p[11], "XT:A:U") == 0) 
		r->uniqmap = 1;
	      else 
		r->uniqmap = 0;
	      process = 1;
	    } else {
	      process = 0; // unmapped skip
	    }

	  }
	  
	  // BOWTIE SAM
	} else if (ri->formatid == BOWTIESAM) { // assumes that reads are mapped unambigously !
	  
	  bitflag = atoi(p[1]);
	  r->pos     = atoi(p[3]);
	  r->st      = SAMbitflag_strand(bitflag);
	  r->lenread = strlen(p[9]);
	  process = 1;	    	    	  
	}
	
	free(p);
	free(line);
	
	if (process == 0) {
	  continue;
	} else {
	  return 1;
	}
      
      } // end parse line
    } // while (1)
  
  } // else other format

  return 0;
}




//
// simply count mutations
// 
void CountMismatchesInSAMFile(const char* file, int readlen, int uniquereads, 
			      char** seqnames, int* seqlens, int numseqs, 
			      int* numreads, long* numnt, long* nummm, int verbose)
{

  char** p;
  int   m;
  int   i =0;
  int   numdigits = 0;
  //int   j = 0;
  int   pos = 0;
  int   st = 0;
  char  nt;
  int            absposmut;
  
  int            rst = 0, ren = 0;
  //int            verbose = 0;

  // printf("we are in getReadCountFromReadFile_SAM()\n");

  char* snum = (char*)calloc(100, sizeof(char));
  //int   bitflag;
  char* mmstr;
  int   mmstrlen = 0;
  char* readseq;
  int   readseqlen;
  int   lb;
  char* cigr;
  char* alignedseq;
  char* myreadpos;
  char* alignedpos;
  int   idxnewread = 0;
  int   idxread    = 0;
  int   cigrlen = 0;
  int   inum = 0;
  int   idxalnread = 0;
  int   insertion = 0;
  int   alignedseqlen = 0;
  int   rennotr;
  int   ntidx = 0;
  int   l;
  char* line;
  LineI li;
  int idxmd = -1;
  unsigned char** a_unqreads_F = 0;
  unsigned char** a_unqreads_R = 0;
  char* seqname = 0;
  int   idx = -1;
  int   chrlen = -1;

  // alloc memory 
  alignedseq = (char*)calloc(readlen * 2, sizeof(char));
  myreadpos  = (char*)calloc(readlen * 2, sizeof(char));
  alignedpos = (char*)calloc(readlen * 2, sizeof(char));

  // DEBUG
  my_aln *aln_seq= NULL, *aln_nt= NULL;
  my_aln *readpos2nt= NULL, *tmpread= NULL;

  // prep work
  if (uniquereads == 1) {
    a_unqreads_F = (unsigned char**)calloc(numseqs, sizeof(unsigned char*));
    a_unqreads_R = (unsigned char**)calloc(numseqs, sizeof(unsigned char*));
  }
  HASH hc;
  HASH_init(&hc, numseqs);
  for (i=0; i<numseqs; i++) {
    HASH_enter(&hc, seqnames[i], i);
    if (uniquereads == 1) {
      a_unqreads_R[i] = create_binarized_array(seqlens[i]);
      a_unqreads_F[i] = create_binarized_array(seqlens[i]);
    }
  }


  *numreads = 0;

  li.bufflen = 10000000;
  LineIopen(&li, (char*)file);
  while ((line = nextLine(&li))) {


    if (line[0] == '@') {
      FREE(line);
      continue;
    }

    //char* line = strdup(prevs);
    char* myline = strdup(line);

    // split
    split_line_delim(line, "\t", &p, &m);
	  
    int bitflag = atoi(p[1]);

    if (!SAMbitflag_ismapped(bitflag)) {
      FREE(p); free(line); free(myline);
      continue;
    }
    
    // look for MD:Z string
    int tmpidx = m-1;
    idxmd      = -1;
    while (tmpidx >= 0) {
      if (strstr(p[tmpidx], "MD:Z") != 0) {
	idxmd = tmpidx;
	break;
      }
      tmpidx --;
    }   	      
    if (idxmd == -1) {
      die("Cannot find MD:Z tag, check format.\n");
    }

    // mmstr
    mmstr = p[idxmd] + 5;

    // identify seq
    seqname = p[2];
    idx     = -1;
    if (!HASH_find(&hc, seqname, &idx)) {
      FREE(p); free(line); free(myline);
      fprintf(stderr, "# odd, could not identify chr %s\n", seqname);
      continue;
    }
    if (idx == -1) {
      die("Problem, idx should not be -1\n");
    }

    // chrlen
    chrlen = seqlens[idx];

    // read pos (starts at 1)
    pos = atoi(p[3]);
	  
    // strand
    st  = SAMbitflag_strand(bitflag);
    //
    // if flag set, and read already there, skip that read
    //
    if (uniquereads == 1) {
      if (((st ==  1) && (get_entry_in_binarized_array(a_unqreads_F[idx], pos) == 1)) ||
	  ((st == -1) && (get_entry_in_binarized_array(a_unqreads_R[idx], pos) == 1))) {
	FREE(p); FREE(line); FREE(myline); continue;
      } else {
	if (st == 1)
	  set_entry_in_binarized_array(a_unqreads_F[idx], pos);
	else
	  set_entry_in_binarized_array(a_unqreads_R[idx], pos);
      }
    }
	  
    // actual read sequence
    readseq    = strdup(p[9]);
    readseqlen = strlen(readseq);
	  
    // parse CIGAR
    cigr       = strdup(p[5]);

    for (i=0; i<readseqlen; i++)
      myreadpos[i] = (char)i;
	  
    idxnewread = 0;
    idxread    = 0;
	  
    numdigits  = 0; // num digits (was k)
    cigrlen    = strlen(cigr);
	  
	  
    int idxmmpos_read= 0;		// match/mismatch position on the read	-> 
    for (i=0; i<cigrlen; i++) {
	    
      if (cigr[i] == 'M') {
	lb = atoi(snum);
	int tmpidx;
	for (tmpidx=0; tmpidx < lb; tmpidx++){
	  if ( (aln_nt = (my_aln*)malloc(sizeof(my_aln))) == NULL){
	    fprintf(stderr, "no memory left!\n");
	    exit(-1);
	  }
	  aln_nt->tpos=     idxnewread + tmpidx;
	  aln_nt->nt=      readseq[idxread + tmpidx ] ;
	  aln_nt->myreadpos= myreadpos[idxread + tmpidx ] ;
	  HASH_ADD_INT(aln_seq, tpos, aln_nt);
	  if ( (tmpread = (my_aln*)malloc(sizeof(my_aln))) == NULL){
	    fprintf(stderr, "no memory left!\n");
	    exit(-1);
	  }
	  // read position (i.e. n-th match/mismatch base) -> nt
	  tmpread->tpos= idxmmpos_read;
	  tmpread->nt= aln_nt->nt ;
	  tmpread->chrpos= idxnewread + tmpidx ;	// position in chromosome 
	  HASH_ADD_INT(readpos2nt, tpos, tmpread);
	  idxmmpos_read++;
	}
	      
	// update both idx
	idxnewread += lb;
	idxread    += lb;
	
	numdigits = 0;  // reset digits
	      
      } else if ( cigr[i] == 'N' ){
	lb= atoi(snum);
	      
	idxnewread+= lb;
	numdigits  = 0;
	      
      } else if (cigr[i] == 'D') {
	// deleting the next lb bases from reference genome
	lb = atoi(snum);
	      
	my_aln *s;
	      
	// add gaps
	// printf("DEPRECATED!\n");
	for (l=idxnewread; l<idxnewread+lb; l++) {
	  HASH_FIND_INT( aln_seq, &l, s );  /* s: output pointer */
	  if (s != NULL){
	    fprintf(stderr, "WARNING! found pos= %d that's not supposed to be deleted\n", l);
	  }
		
	}
	      
	// update index in new read
	idxnewread += lb;
	      
	numdigits = 0; // reset digits
	      
      } else if (cigr[i] == 'I') {
	lb = atoi(snum);
	      
	// update idx in read
	idxread  += lb;
	      
	numdigits = 0; // reset digits
	      
      } else if (cigr[i] == 'S') {
	lb = atoi(snum);
	numdigits = 0; // reset digits
	      
      } else {
	// num
	snum[numdigits] = cigr[i]; // add digits
	numdigits++;
	snum[numdigits] = '\0'; // add end char
      }
    }
	  
	  
    if (st == -1) {
      my_aln *s;
      for(s= aln_seq; s != NULL; s=s->hh.next) {
	if (s->myreadpos != -1)
	  // turn around others
	  s->myreadpos = (char)( readseqlen - s->myreadpos - 1 );
      }
	    
    }
	  
    alignedseqlen = HASH_COUNT(aln_seq);
	  
    // count aligned nucleotides
    rst     = pos;
    ren     = min(pos + alignedseqlen, chrlen );
    rennotr = min(pos + alignedseqlen, chrlen);

    int numalignednts = 0;
    //int posinread     = 0;
	  
    // DEBUG: continue here
    // go through all keys of aln_seq
    // i= pos+ truncate + key
    my_aln *s;
	  
    // go through all the reads and register on the chromosome where the reads come from (i.e. n-th base in a read)
    for(s= aln_seq; s != NULL; s=s->hh.next) {
      if  (s->nt != '-') { 
	numalignednts ++;	      
      }      
    }
	  
    (*numnt) += numalignednts;
	  
    // mmstr -> to uppercase
    for (i= 0; mmstr[i]; i++)
      mmstr[i]= toupper(mmstr[i]);
	  
    //
    // parse MUTATIONS
    //
    i          = 0; // iterate thru mmstr
    numdigits  = 0; // num digits
    idxalnread = 0; // index in aligned read
    insertion  = 0; // are we in an insertion
    mmstrlen   = strlen(mmstr);
	  

    while (i<mmstrlen) {
	    
      if ((mmstr[i] == 'A') ||
	  (mmstr[i] == 'C') ||
	  (mmstr[i] == 'G') ||
	  (mmstr[i] == 'T')) {

	if (numdigits > 0) {
	  snum[numdigits] = '\0';
	  inum = atoi(snum);
	  idxalnread += inum;
	  numdigits = 0; // reset
	}

	// if we are not in an insertion (or a deletion in the read)
	if (insertion == 0) {

	  // make sure the readpos is not outside sequence
	  if (pos + idxalnread <= chrlen) {

	      //	want to find out which nucleotide the base is mutated to...
	      //	should match positions after taking into acocunt insertions !!!

	      HASH_FIND_INT( readpos2nt, &idxalnread, s );  // s: output pointer 
	      absposmut = pos + s->chrpos - 1;
	      nt        = s->nt;
	      ntidx     = data_nucl(nt);
	      // not an N or IUPAC code
	      if (ntidx >= 0) {
		(*nummm) ++;
		
	      } // if (ntidx

	  } // if mut pos within chrom range
	  
	} // if not in insertion

	// in any case, move by one nt
	//	DEBUG: 
	//		should consider the case of deletion
	//		e.g. ^AC means deleting 2 bases from the reference genome
	idxalnread ++;
	// printf("\tadvancing by 1 nt (ref= %c), current idxalnread= %d\n", mmstr[i], idxalnread);

      } else if (mmstr[i] == '^') {
	if (numdigits > 0) {
	  snum[numdigits] = '\0';
	  inum = atoi(snum);
	  idxalnread += inum;
	  numdigits = 0; // reset
	}
	
	// new:
	//	look for trailing nucleotides (those deleted from the reference genome
	int tmpidx=0;
	for (tmpidx= i+1; mmstr[tmpidx]; tmpidx++){
	  // are there trailing nucleotide(s)?
	  // printf("\tmmstr[%d]= %c => numeric= %d\n", tmpidx, mmstr[tmpidx], isdigit(mmstr[tmpidx]));
	  if (isdigit(mmstr[tmpidx])){
	    // found the end of trailing, deleted nucleotides	
	    break;
	  }
	}

	// printf("\tlast index= %d, should skip i by %d\n", tmpidx, (tmpidx-i-1));
	i+= (tmpidx-i-1);

	insertion = 1; // we entered an insert
	
      } else if ( mmstr[i] == 'N' ){
	if (numdigits > 0){
	  snum[numdigits]= '\0';
	  inum = atoi(snum);
	}

	if (mmstr[i+1] == '0' ){
	  //
	  inum++;
	}

	//if (verbose == 1)
	//	printf("MMSTR: advanced by %d nt\n", inum);

	idxalnread+= inum;
	numdigits  = 0;

      } else if ( isalpha(mmstr[i])){
	fprintf(stderr, "mmstr format not supported %c\n", mmstr[i]);
					
      } else {
	// must be number
	snum[numdigits] = mmstr[i];
	numdigits++;
	// leave insertion if we are in one
	if (insertion == 1) {
	  insertion = 0;
	}
      }
    
      i++;
    }
  
    // printf("end of MMSTR loop\n");

    // final advance
    if (numdigits > 0) {
      snum[numdigits] = '\0';
      inum = atoi(snum);
      idxalnread += inum;
      //if (verbose == 1)
      //	printf("MMSTR: advance by %d nt\n", inum);
      numdigits = 0;

    }
  
    // get length of aligned sequence
     alignedseqlen = HASH_COUNT(aln_seq);
  
  
    //if (verbose)
    //	printf("[hash] alignedseqlen= %d\n", alignedseqlen);
  
    // if (idxalnread != strlen(alignedseq)) 
    if (idxalnread != alignedseqlen) {
      printf("MMSTR: Houston, we have a pb in getReadCountFromReadFile_SAM(), idxalnread= %d, len= %d, pos= %d, mmstr= %s, cigar= %s, file= %s\n", idxalnread, alignedseqlen, pos, mmstr, cigr, file);
    }
  
    FREE(p); // p[x] are just pointers
    FREE(line);
    FREE(myline);

    // clean up and reclaim memory 
    my_aln *current_base;
    while (aln_seq){
      current_base= aln_seq;
      HASH_DEL(aln_seq, current_base);
      FREE(current_base);
    }
    while (readpos2nt){
      current_base= readpos2nt;
      HASH_DEL(readpos2nt, current_base);
      FREE(current_base);
    }

    FREE(cigr);
    FREE(readseq);

					 			
    (*numreads) ++;
				
    if ((verbose == 1) && ( ((*numreads) % 100000) == 0)) {
      printf("Read %d reads from %s                 \r", *numreads, file);
      fflush(stdout);
    }
  }



  if (verbose == 1) 
    printf("Number of reads = %d\n", *numreads);
  
  LineIclose(&li);

  if (uniquereads == 1) {
    for (i=0; i<numseqs; i++) {
      FREE(a_unqreads_F[i]);
      FREE(a_unqreads_R[i]);
    }
  }

  if (verbose == 1) 
    printf("About to exit after reading %s\n", file);

  //fclose(fp);
  FREE(snum);
  FREE(alignedseq);
  FREE(myreadpos);
  FREE(alignedpos);

  if (verbose == 1) 
    printf("About to exit after reading %s (2)\n", file);
  
}
