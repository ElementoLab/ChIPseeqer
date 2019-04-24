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
} MappedRead;

int ReadIopen(ReadI* ri, char* file, int format);
void ReadIclose(ReadI* ri);
int nextRead(ReadI* ri, MappedRead* r);
