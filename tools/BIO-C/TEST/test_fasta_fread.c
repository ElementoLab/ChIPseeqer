#include <stdio.h>


typedef struct _seqI {
  
  FILE* fp;
  int   nextSequence_started;  // boolean
  int   nextSequence_ended;    // boolean
  char* nextSequence_currentLine;
  int   seqlen_inc;
  int   verbose;

  // fast 
  char* seq;
  int   pos;
  int   max_seqlen;

} seqI;


char* seqI_nextSequence_fast(seqI* s, char** name, int* size);
int seqI_open_fast(seqI* s, char* file);
void seqI_set_max_seqlen(seqI* s, int n);
void seqI_set_seqlen_inc(seqI* s, int n);

int main(int argc, char** argv) {



  char* name;
  int   size;
  seqI  si;
  char* myseq;

  seqI_set_max_seqlen(&si, 150000000);
  seqI_set_seqlen_inc(&si, 30000000);
  seqI_open_fast(&si, argv[1]);
  
  while ( myseq = seqI_nextSequence_fast(&si, &name, &size) ) {
   
    printf("%d\n", strlen(myseq));
    
  }
  


}


void seqI_set_seqlen_inc(seqI* s, int n)
{
  s->seqlen_inc = n;

}


int seqI_open_fast(seqI* s, char* file)
{
  
  int   cnt       = 0;
  int   total_cnt = 0;
  char* buff;
  char* seq;
  int   bufflen   = 1000000;

  s->fp = fopen(file, "r");
  if (s->fp == 0) {
    return 0;
  }  

  seq    = (char*)calloc(s->max_seqlen, sizeof(char));
  if (seq == 0) {
    printf("cannot allocate seq with n=%d\n", s->max_seqlen);
    exit(0);
  }
  buff   = (char*)calloc(bufflen, sizeof(char)); 

  while (!feof(s->fp)) {
    cnt = fread(buff, 1, bufflen, s->fp);    
    //printf("read %d char.\n", cnt);
    total_cnt += cnt;
    strcat(seq, buff);    
  }
  
  seq[total_cnt] = '\0';

  printf("len = %d\n", strlen(seq));
    
  s->seq = seq;
  s->pos = 0;
  
  s->verbose    = 0;
  return 1;
}





void seqI_set_max_seqlen(seqI* s, int n)
{
  s->max_seqlen = n;

}

char* seqI_nextSequence_fast(seqI* s, char** name, int* size)
{

  int l = strlen(s->seq);
  int i = 0, j = 0;
  int inseq = 0;
  char* seq = 0;

  
  

  j = 0;
  while (s->pos < l) {

    if (s->seq[ s->pos ] == '>') {
  
      
    
      // two cases, if not yet started, enter seq
      if (inseq == 0) {
	
	printf("allocating %d\n", s->seqlen_inc);
	
	seq = (char*)calloc(s->seqlen_inc, sizeof(char));
	if (seq == 0) {
	  printf("Problem allocatin seq, n=%d\n", s->seqlen_inc);
	  exit(0);
	}

	*name = (char*)calloc(100, sizeof(char));
	i     = 0;
	s->pos++;
	while (s->seq[ s->pos ] != '\n') {
	  (*name)[i] = s->seq[ s->pos ];
	  s->pos++;
	  i ++;
	}

	(*name)[i] = '\0';
	
	printf("name = %s\n", *name);
	fflush(stdout);

	inseq = 1;

      } else {
	
	// means we are reaching the begining of a new seq
	// stay at the same position in file
	
	return seq;
      }

    } else {

      // add to seq
      if ((s->seq[ s->pos ] != '\n') && (s->seq[ s->pos ] != '\r')) {	
	seq[j] = s->seq[ s->pos ];
	j ++;
      }
      
      s->pos ++;
    }
    
    
  }

  return seq;
  
}


