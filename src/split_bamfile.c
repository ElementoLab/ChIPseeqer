#include <stdio.h>
#include <search.h>
#include <sys/stat.h>
#include "lib/third_party/samtools/sam.h"
#include "lib/third_party/samtools/bam.h"

#include "dataio.h"
#include "hashtable.h"

#define bam1_unmapped(b) (((b)->core.flag & BAM_FUNMAP) != 0)
#define bam1_mate_unmapped(b) (((b)->core.flag & BAM_FMUNMAP) != 0)
#define bam1_paired(b) (((b)->core.flag & BAM_FPAIRED) != 0)
 

int main(int argc, char** argv) {

  if (argc < 2) {
    printf("Usage: split_bamfile FILE [ -outdir DIR -seqlist FILE ] \n");
    exit(0);
  }

  fprintf(stderr, "Warning: this tool only supports splitting 1 BAM file\n");

  char* outdir = 0;
  if (exist_parameter(argc, argv, "-outdir"))
    outdir  = get_parameter(argc, argv, "-outdir");
  
  char* file = argv[1];

  if (outdir == 0) {
    outdir = (char*)calloc(1000, 1);
    sprintf(outdir, "%s_SPLIT", file);
    mkdir(outdir, S_IRWXU | S_IRWXO | S_IRWXG);
  }
  
  int i;
  samfile_t* in;
//int prevpos = -1;
//int prevst  = -1;
  samfile_t** a_out;
  samfile_t* out;

  int maxoutfile = 100;
  int numoutfile = 0;
  // reserve mem for outfile
  a_out = (samfile_t**)calloc(maxoutfile, sizeof(samfile_t*));
  HASH hc;
  int cidx = -1;
  HASH_init(&hc, maxoutfile);
  
  in = samopen(file, "rb", 0);
  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", file);
    return 1;    
  }
        
  bam1_t* b = bam_init1();
  int ret;
  char* seq = 0;
  char tmpfile[1000];
  int numlines = 0;
  while ((ret = samread(in, b)) >= 0) {
    numlines ++;
    if ((numlines % 10000) == 0) {
      printf("Read %d lines         \r", numlines);
      fflush(stdout);
    }
    if (bam1_unmapped(b)) // no need to register these
      continue;

    // get seq
    seq = in->header->target_name[b->core.tid];

    if (HASH_find(&hc, seq, &cidx)) {
      out = a_out[cidx];
    } else {
      if (numoutfile == maxoutfile) {
	die("maximum number of chrom reached\n");
      }
      HASH_enter(&hc, seq, numoutfile);
      sprintf(tmpfile, "%s/read.%s.bam", outdir, seq);
      out = samopen(tmpfile, "wb", in->header);
      if (out == 0) {
	fprintf(stderr, "Fail to open BAM file %s for writing\n", tmpfile);
	return 1;    
      }
      a_out[numoutfile] = out;      
      numoutfile ++;
    }
    
    samwrite(out, b);
    
    
  }

  samclose(in);
  for (i=0; i<numoutfile; i++) {
    samclose(a_out[i]);
  }
  return 0;
}
  
    
