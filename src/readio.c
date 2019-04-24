#include <stdio.h>
#ifndef SAM_H
#include "lib/third_party/samtools/sam.h"
#include "lib/third_party/samtools/bam.h"
#define SAM_H
#endif
#include "sequences.h"
#include "dataio.h"
#include "readio.h"

#define bam1_unmapped(b) (((b)->core.flag & BAM_FUNMAP) != 0)
#define bam1_mate_unmapped(b) (((b)->core.flag & BAM_FMUNMAP) != 0)
#define bam1_paired(b) (((b)->core.flag & BAM_FPAIRED) != 0)


int ReadIopen(ReadI* ri, char* file, int format) {

  ri->formatid = format;
  ri->verbose  = 0;

  //if (strcmp(format, "bam") == 0) {
  //  ri->formatid = BAM;
  //}
  
  if (ri->formatid == BAM) {    
    ri->in = samopen(file, "rb", 0);
    if (ri->in == 0) {
      fprintf(stderr, "Fail to open BAM file %s\n", file);
      return 0;
    }
    ri->b = bam_init1();
    ri->inread = 0; // not in a read (need to get one)
    return 1;
  } else { // txt format
    //ri->li.bufflen = 10000000;
    return LineIopen(&(ri->li), file);
  }



}

void ReadIclose(ReadI* ri) {
  if (ri->formatid == BAM) {
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

  if (ri->formatid == BAM) {
    
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

	xt = (char*)(bam_aux_get(ri->b, "XT"));
	if ((xt == 0) || (xt[0] != 'A')) {
	  die("Problem interpreting XT field\n");
	}
	if (xt[1] == 'U')
	  r->uniqmap = 1;
	else 
	  r->uniqmap = 0; // R or M

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
	  
	  r->uniqmap = 1;

	  //if (p[3][0] == 'U') {
	  //  r->uniqmap = 1;
	  //} else {
	  //  r->uniqmap = 0;
	  //}

	  process = 1;
	  
	  // SAM
	} else if (ri->formatid == SAM) {

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
