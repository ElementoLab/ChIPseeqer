#include <search.h>
#include "hashtable.h"
#include "uthash.h"

#ifndef SAM_H
#include "lib/third_party/samtools/faidx.h"
#include "lib/third_party/samtools/sam.h"
#include "lib/third_party/samtools/bam.h"
#define SAM_H
#endif

#define ELAND 0
#define MIT 1
#define BED 2
#define SAM 3
#define BOWTIESAM 4
#define BAM 5
#define BAMUNQ 6
#define STARTONLY -10
 

// data structure for aligned read
typedef struct my_aln {
  int tpos;                      /* key */
  char nt;			// nucleotide at the position
  int myreadpos; 		// position on the read?
  int chrpos; 			// position on the read?
  char qs;                      // quality score
  UT_hash_handle hh;            /* makes this structure hashable */
} my_aln;


typedef struct _ReadI {
  
  int formatid;
  
  /* for all txt formats */
  LineI li;
  
  /* BAM */
  samfile_t* in;
  bam1_t* b;
  int inread;
  int verbose;


} ReadI;


typedef struct _MappedRead {
  char* seqname;
  int   pos;
  int   st;
  char* readname;
  int   uniqmap;
  int   lenread;
  uint32_t* cigar;
  int n_cigar;
} MappedRead;


char* my_fai_fetch(faidx_t* idx, char* c, int i, int j);

void getQualityScoresFromReadFile_SAM(const char* file, char* chrname, const int chrlen, unsigned char* snvpos, unsigned char**** qualscores, unsigned short*** numqualscores, HASH* hsnv, const int uniquereads);

int ReadIopen(ReadI* ri, char* file, int format);
void ReadIclose(ReadI* ri);
int nextRead(ReadI* ri, MappedRead* r);

int formatToId(char* formattext);

void parseTopHatCIGR(char* cigr, int** blocksizes, int** blockoffsets, int* numblocks);

int guessReadLength(char* readfile, char* format, int head);

//
// ELAND
//
void getSimpleReadCountFromReadFile(char* file, int chrlen, unsigned short** readcounts, int* numreads, int uniquereads, int truncate);
void getReadCountFromReadFile(char* file, int chrlen, unsigned short** readcounts, unsigned short** mutacounts, unsigned short** mutacounts_F,
                              int* numreads, long* numnt, long* nummm, unsigned char*** ntcounts, unsigned char*** readpos, int uniquereads, int rw, int truncate);
void updateMismatchCountAtReadPos(char* file, int chrlen, int** poscounts, long* numnt, long* nummm, int uniquereads);


//
// SAM
//

void getINDELCountFromSAMFile(const char* file, const int chrlen, unsigned short** readcounts, unsigned short** mutacounts, unsigned short** mutacounts_F,
			     int* numreads, long* numnt, long* nummm, int del, unsigned short** indellen, int uniquereads, const int readlen, int uniqmap, int verbose);

void getSimpleReadCountFromReadFile_SAM(const char* file, const int chrlen, unsigned short** readcounts, int* numreads, const int uniquereads, const int truncate, const int readlen);
void getReadCountFromReadFile_SAM(const char* file, const int chrlen, unsigned short** readcounts, unsigned short** mutacounts, unsigned short** mutacounts_F,
                                  int* numreads, long* numnt, long* nummm, unsigned char*** ntcounts, unsigned char*** readpos, int uniquereads, int rw, int truncate, const int readlen);
int CountReadsInChrInterval(char* file, int format, char* chr, int i1, int i2);

void updateMismatchCountAtReadPos_SAM(const char* file, const int chrlen, int** poscounts, long* numnt, long* nummm, const int uniquereads, const int readlen);

int SAMbitflagNum(int bitflag, int num);
int SAMbitflag_strand(int bitflag);
int SAMbitflag_isPairedProper(int bitflag);
int SAMbitflag_ismapped(int bitflag); 
//int SAMbitflag_isPairedProper(int bitflag)

void getCountFromReadFile(char* file, int format, char* chr, long chrlen, int readlen, int fraglen, int uniquereads, int uniqmap, unsigned short** counts, int* numreads, int* numclonalreads);


int CountReads(char* file, int format);

void CountMismatchesInSAMFile(const char* file, int readlen, int uniquereads, 
			      char** seqnames, int* seqlens, int numseqs, 
			      int* numreads, long* numnt, long* nummm, int verbose);
