//
// Author : Olivier Elemento, Princeton University
//
// Permission is granted to do whatever you like with this code, provided you accept that
// I cannot be help responsible for any misuse of it or of the program it codes for.
//
//
#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)
#define MAXMOTIFS 10000
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct _ACE_motif {
  float* v;
  int    w;
  int*   s;

} ACE_motif;

float compare_motifs(float *m1, int l1, float* m2, int l2, int ss); 


void  readMotif(char* file, float** m, int** st, int* l);
int   nucl(char c); 
float pearson(float* data1, float* data2, int m); 
void  chomp(char* s);
int*  bubbleSortIndex(float *ia, int n);
int   exist_flag(int argc, char** argv, char* param);
void  readMotifs(char* file, ACE_motif** a_am, int* n);


int main(int argc, char** argv) { 

  //int    l1, l2;
  //float* m1, *m2, *m3;
  //int*   st1, *st2, *st3;
  
  int    verbose = 0;
  int    ss = 0;

  ACE_motif* a_m1;
  ACE_motif* a_m2;
  int nm1, nm2;
  float pc = -1.0;
  int   i,j;
  int   simple = 0;

  ss      = exist_flag(argc, argv, "-ss");
  verbose = exist_flag(argc, argv, "-verbose");
  simple  = exist_flag(argc, argv, "-simple");

  //readMotif(argv[1], &m1, &st1, &l1);
  //readMotif(argv[2], &m2, &st2, &l2); 
 
  readMotifs(argv[1], &a_m1, &nm1);
  readMotifs(argv[2], &a_m2, &nm2); 
  
  if (verbose == 1) {
    printf("Read %d and %d motifs ... \n", nm1, nm2);
  }
  
  if (simple == 1) {
    pc = compare_motifs(a_m1[0].v, a_m1[0].w, a_m2[0].v, a_m2[0].w, ss);
    printf("%4.3f\n", pc);

  } else {
    
    for (i=0; i<nm1; i++) {
      for (j=0; j<nm2; j++) {
	pc = compare_motifs(a_m1[i].v, a_m1[i].w, a_m2[j].v, a_m2[j].w, ss);
	printf("%d\t%d\t%4.3f\n", i, j, pc);
      }
    }
  }
  
}

float compare_motifs(float *m1, int l1, float* m2, int l2, int ss) 
{

  int    min_inf_cols = 6;
  int    verbose      = 0;
  int   *st1, *st2, *st3;
  float* in1;
  float* in2;
  int*   in1_idx;
  int*   in2_idx;
  int    cnt_st1, cnt_st2;
  float  lg2;
  int    i, j;
  float  *m3;
  int    i1, i2;
  float  pe, pe_max;
  int    cnt;
  float* V1, *V2;
  int    start_i1, start_i2;
  float  min_inf_t    = 0.8;

  float  last_i;

  lg2 = log(2.0);

  // minimum number of informative columns is of course bounded by motif width
  min_inf_cols = min(min_inf_cols, l1);
  min_inf_cols = min(min_inf_cols, l2);
  
  st1 = (int*)calloc( l1, sizeof(int));  
  st2 = (int*)calloc( l2, sizeof(int));

  if (verbose == 1)
    printf("min_inf_cols = %d\n", min_inf_cols);

  // calculate information content of each column of motif 1
  in1 = (float*)calloc(l1, sizeof(float));
  for (i=0; i<l1; i++) {
    in1[i] = 2.0;
    for (j=0; j<4; j++) {
      in1[i] += m1[i*4+j] * log( m1[i*4+j] ) / lg2;
    }
    //printf("pos = %d, I=%4.1f\n", i, in1[i]);
  }


  // flag the 6+ most informative columns
  in1_idx = bubbleSortIndex(in1, l1);
  for (i=0; i<min_inf_cols; i++) {
    //printf("I = %4.3f (%d, 1)\n", in1[ in1_idx[i] ], in1_idx[i]);

    st1[ in1_idx[i] ] = 1;
    last_i            = in1[ in1_idx[i] ];

  }
  i=min_inf_cols;
  //while ((i<l1) && ( ( in1[ in1_idx[min_inf_cols-1] ] - in1[ in1_idx[i] ] ) < 0.001 )) {
  while ((i<l1) && ( in1[ in1_idx[i] ] >= min_inf_t)) {

    //printf("5 - %d = %4.3f\n", i, ( in1[ in1_idx[5] ] - in1[ in1_idx[i] ] ));    
    //printf("I = %4.3f (%d, 1)\n", in1[ in1_idx[i] ], in1_idx[i]);
    st1[ in1_idx[i] ] = 1;
    last_i            = in1[ in1_idx[i] ];

    i++;
  }
  for (j=i; j<l1; j++) {
    //printf("I = %4.3f (%d, 0)\n", in1[ in1_idx[j] ], in1_idx[j]);
    if ( (last_i - in1[ in1_idx[i] ]) < 0.001)
      st1[ in1_idx[j] ] = 1;
    else
      st1[ in1_idx[j] ] = 0;
  }
  
  /*
  for (i=0; i<l1; i++) {
    printf(" %4d", st1[i]);
  }
  printf("\n\n");
  */

  // same for motif 2
  in2 = (float*)calloc(l2, sizeof(float));
  for (i=0; i<l2; i++) {
    in2[i] = 2.0;
    for (j=0; j<4; j++) {
      in2[i] += m2[i*4+j] * log( m2[i*4+j] ) / lg2;
    }
    //printf("pos = %d, I=%4.1f\n", i, in2[i]);
  }

  // flag the 6+ most informative columns
  in2_idx = bubbleSortIndex(in2, l2);
  for (i=0; i<min_inf_cols; i++) {
    //printf("I = %4.3f (%d, 1)\n", in2[ in2_idx[i] ], in2_idx[i]);    
    st2[ in2_idx[i] ] = 1;
    last_i            = in2[ in2_idx[i] ];

  }

  // if (m
  i=min_inf_cols;
  while ((i<l2) && ( in2[ in2_idx[i] ] >= min_inf_t)) {

  //while ((i<l2) && ( ( in2[ in2_idx[min_inf_cols-1] ] - in2[ in2_idx[i] ] ) < 0.001 )) {    
    //printf("5 - %d = %4.3f\n", i, ( in2[ in2_idx[5] ] - in2[ in2_idx[i] ] ));    
    //printf("I = %4.3f (%d, 1)\n", in2[ in2_idx[i] ], in2_idx[i]);

    st2[ in2_idx[i] ] = 1;
    last_i            = in2[ in2_idx[i] ];

    i++;
  }

  for (j=i; j<l2; j++) {
    //printf("I = %4.3f (%d, 0)\n", in2[ in2_idx[j] ], in2_idx[j]);	
    //st2[ in2_idx[j] ] = 0;
    
    if ( (last_i - in2[ in2_idx[i] ]) < 0.001)
      st2[ in2_idx[j] ] = 1;
    else
      st2[ in2_idx[j] ] = 0;

  }
  
  /*
  for (i=0; i<l2; i++) {
    printf(" %4d", st2[i]);
  }
  printf("\n\n");
  */

  //for (i=0; i<l2; i++) {
  //  printf(" %3.2f", in2[ in2_idx[i] ]);
  //}
    
  if (verbose == 1) {

    printf("l1=%d\n", l1);
    for (i=0; i<l1; i++) {
      for (j=0; j<4; j++) {
	printf(" %4.1f", m1[i*4+j]);
      }
      printf("\t%3.2f\t%4d\n", in1[i], st1[i]);
    }
    
    printf("l2=%d\n", l2);
    for (i=0; i<l2; i++) {
      for (j=0; j<4; j++) {
	printf(" %5.2f", m2[i*4+j]);
      }
      printf(" %5d\n", st2[i]);
    }

  }

  
  V1 = (float*)malloc(4*l1 * sizeof(float));
  V2 = (float*)malloc(4*l2 * sizeof(float));

  //
  // m1 |    x    m2 -
  //
  // diago starting on m1
  pe_max = 0.0;
  for (start_i1=0; start_i1<l1; start_i1++) {

    i1      = start_i1;
    i2      = 0;
    cnt     = 0;						
    cnt_st1 = 0;
    cnt_st2 = 0;

    if (verbose == 1)	
      printf("Aln starting at i1=%d, i2=%d\n", i1, i2);

    while ((i1 < l1) && (i2 < l2)) {
      
      //printf("%d-%d ", i1, i2);
      if (st1[i1] == 1)
	cnt_st1 ++;
      if (st2[i2] == 1)
	cnt_st2 ++;
	
      /*
      printf("%d-%d\t", i1, i2);      
      for (j=0; j<4; j++) {
	printf(" %4.2f", m1[i1*4+j]);
      }
      printf(" and ");
      for (j=0; j<4; j++) {
	printf(" %4.2f", m2[i2*4+j]);
      }   
      printf("\n");
      */

      for (j=0; j<4; j++) {
	V1[cnt] = m1[i1*4+j];
	V2[cnt] = m2[i2*4+j];      
	cnt ++;      
      }
      
      //}

      i1 ++;
      i2 ++;
      
    }
    //printf(", cnt=%d\n", cnt/4);

    if (verbose == 1)
      printf("cnt_st1=%d and cnt_st2=%d\n", cnt_st1, cnt_st2);
	
    
    if ((cnt_st1 >= min_inf_cols) && (cnt_st2 >= min_inf_cols)) {
      
      pe = pearson(V1, V2, cnt);
      
      /*
      for (i=0; i<cnt; i++) {
	printf("%f\t%f\n", V1[i], V2[i]);
      }
      */

      if (verbose == 1)
	printf("pe = %4.3f\n", pe);
      
      if (pe > pe_max) {
	pe_max = pe;
      }
      
    }	
  }

  
  // diago starting on m2
  for (start_i2=0; start_i2<l2; start_i2++) {
    
    i1  = 0;
    i2  = start_i2;
    cnt = 0;
    cnt_st1 = 0;
    cnt_st2 = 0;

    if (verbose == 1)
      printf("Aln starting at i1=%d, i2=%d\n", i1, i2);
  
    while ((i1 < l1) && (i2 < l2)) {
      

      if (st1[i1] == 1)
	cnt_st1 ++;
      if (st2[i2] == 1)
	cnt_st2 ++;

      /*
      printf("%d-%d\t", i1, i2);      
      for (j=0; j<4; j++) {
	printf(" %4.2f", m1[i1*4+j]);
      }
      printf(" and ");
      for (j=0; j<4; j++) {
	printf(" %4.2f", m2[i2*4+j]);
      }   
      printf("\n");
      */
 
      for (j=0; j<4; j++) {
	V1[cnt] = m1[i1*4+j];
	V2[cnt] = m2[i2*4+j];      
	cnt ++;      
      }
      
      i1 ++;
      i2 ++;
      

    }

    //}
    
    //printf(", cnt=%d\n", cnt/4);

    //if ((cnt / 4) >= 6) {

    if (verbose == 1)
      printf("cnt_st1=%d and cnt_st2=%d\n", cnt_st1, cnt_st2);


    if ((cnt_st1 >= min_inf_cols) && (cnt_st2 >= min_inf_cols)) {

      pe = pearson(V1, V2, cnt);

      if (verbose == 1)
	printf("pe = %4.3f\n", pe);
      
      if (pe > pe_max) {
	pe_max = pe;
      }
    
    }	
  }

  if (ss == 0) {

    // take rev comp of motif 2
    m3  = (float*)malloc(l2*4 * sizeof(float));
    st3 = (int*)  malloc(l2 * sizeof(int  ));
    
    for (i=l2-1,j=0; i>=0; i--,j++) {
      
      //printf("i=%d, j=%d\n", i, j);
      
      // A == T
      m3[ j*4 + 0 ] = m2[ i*4 + 3 ];
    
      // T == A
      m3[ j*4 + 3 ] = m2[ i*4 + 0 ];
      
      // C == G
      m3[ j*4 + 1 ] = m2[ i*4 + 2 ];
      
      // G == C
      m3[ j*4 + 2 ] = m2[ i*4 + 1 ];
      
      st3[ j ] = st2[ i ];
      
    }


      if (verbose == 1) {

	printf("l3=%d\n", l2);
	for (i=0; i<l2; i++) {
	  for (j=0; j<4; j++) {
	    printf(" %4.1f", m3[i*4+j]);
	  }
	  printf(" %4d\n", st3[i]);
	}
      }
    
    
  
    
    for (start_i1=0; start_i1<l1; start_i1++) {
      
      i1  = start_i1;
      i2  = 0;
      cnt = 0;
      cnt_st1 = 0;
      cnt_st2 = 0;

      if (verbose == 1)
	printf("Aln starting at i1=%d, i2=%d\n", i1, i2);
      
      
      while ((i1 < l1) && (i2 < l2)) {
      
	//printf("%d-%d ", i1, i2);            
	//if ((st1[i1] == 1) && (st3[i2] == 1)) {

	if (st1[i1] == 1)
	  cnt_st1 ++;
	if (st3[i2] == 1)
	  cnt_st2 ++;
	
	for (j=0; j<4; j++) {
	  V1[cnt] = m1[i1*4+j];
	  V2[cnt] = m3[i2*4+j];      
	  cnt ++;      
	}
      
	//}

	i1 ++;
	i2 ++;
      
      }

      // printf(", cnt=%d\n", cnt/4);
    
      // if ((cnt / 4) >= 6) {

      if (verbose == 1)
	printf("cnt_st1=%d and cnt_st2=%d\n", cnt_st1, cnt_st2);


      if ((cnt_st1 >= min_inf_cols) && (cnt_st2 >= min_inf_cols)) {

	pe = pearson(V1, V2, cnt);
      
	if (verbose == 1)	
	  printf("pe = %4.3f\n", pe);
      
	if (pe > pe_max) {
	  pe_max = pe;
	}
      
      }	
    }

  
    // diago starting on m2
    for (start_i2=0; start_i2<l2; start_i2++) {
    
      i1  = 0;
      i2  = start_i2;
      cnt = 0;
      cnt_st1 = 0;
      cnt_st2 = 0;


      if (verbose == 1)
	printf("Aln starting at i1=%d, i2=%d\n", i1, i2);

  
      while ((i1 < l1) && (i2 < l2)) {
      
	// printf("%d-%d ", i1, i2);      
	// if ((st1[i1] == 1) && (st3[i2] == 1)) {

	if (st1[i1] == 1)
	  cnt_st1 ++;
	if (st3[i2] == 1)
	  cnt_st2 ++;

	for (j=0; j<4; j++) {
	  V1[cnt] = m1[i1*4+j];
	  V2[cnt] = m3[i2*4+j];      
	  cnt ++;      
	}
      
	//}
      
	i1 ++;
	i2 ++;
      
      }
    
      // printf(", cnt=%d\n", cnt/4);

      //if ((cnt / 4) >= 6) {

      if (verbose == 1)
	printf("cnt_st1=%d and cnt_st2=%d\n", cnt_st1, cnt_st2);
    
      if ((cnt_st1 >= min_inf_cols) && (cnt_st2 >= min_inf_cols)) {

	pe = pearson(V1, V2, cnt);

      
	if (verbose == 1)
	  printf("pe = %4.3f\n", pe);

	if (pe > pe_max) {
	  pe_max = pe;
	}

      }  
      //}	
    }

  }

  //printf("%f\n", pe_max);
  
  return pe_max;
}

void chomp(char* s)
{

  int j;

  j = strlen(s);
  while (j >= 0) {
    if (s[j] == '\n') {
       s[j] = '\0';
       return;
    }
    j--;
  }


}

//
// this function reads a motif file (motifs separated by empty lines)
//
void readMotifs(char* file, ACE_motif** a_am, int* n)
{

  char* buff = 0;
  FILE* fp;
  int   maxlen = 1000;
  int   i, len;
  int   cnt_stars, cnt_lines = 0;
  int   mw = 0;
  float* v;
  int*  myst = 0;
  int   c;
  int   k;
  int   nummotifs = 0;

  buff = (char*)calloc(maxlen, sizeof(char));
  
  fp    = fopen(file, "r");
  if (!fp) {
    printf("Cannot open motif file %s. Exiting.\n", file);
    exit(0);
  }
  
  
  *a_am = (ACE_motif*)malloc(MAXMOTIFS * sizeof(ACE_motif));
  
  while (!feof(fp)) {

    fgets(buff, maxlen, fp);


    len = -1;
    
    if (!feof(fp)) {
      chomp(buff);
      len = strlen(buff);
    }


    if (feof(fp) || (len == 0)) {
      
      //printf("store motif .. \n");
      
      // calculate frequencies here, using a pseudocount-like approach
      for (i=0; i<mw*4; i++) {
	v[i] = ( v[i] + 0.5 ) / (cnt_lines + 0.5);
      }  
 
      // add the motif to the list 
      (*a_am)[nummotifs].v = v;
      (*a_am)[nummotifs].w = mw;
      (*a_am)[nummotifs].s = myst;
      
      nummotifs ++;

    }

    if (feof(fp)) {
      break;
    }


  
    if (buff[0] == 'M') {
      // motif line, initialize
      cnt_lines = 0;
      

    } else {
      
      cnt_stars = 0;
      for (i=0; i<len; i++) {
	if (buff[i] == '*') {
	  cnt_stars ++;
	}
      }


      // star line or nucleotide line ?
      if (cnt_stars > 0) {
	
	//printf(" Star line\n");
	if (mw <= 0) {
	  printf("pb .. mw=%d should be > 0 \n", mw);
	}
	
	myst = (int*)calloc(mw, sizeof(int));
	for (i=0; i<mw; i++) {
	  if (buff[i] == '*') {
	    myst[i] = 1;	  
	  } else {
	    myst[i] = 0;
	  }	  
	}
	
	//*st = myst;

	// nt line
      } else {
	
	if (cnt_lines == 0) {
	  // estimate the width of the motif (mw)
	  mw = 0;
	  while ((mw < len) && (buff[mw] != '\t') && (buff[mw] != ' ')) {
	    mw++; 
	  }
	  //*l = mw;
	  v = (float*)calloc(mw * 4, sizeof(float));
	}
	
	i = 0;
	while ((i < len) && (buff[i] != '\t') && (buff[i] != ' ')) {
	  c = nucl( buff[i] );
	  if (c == -1) {
	    //printf("c should not be -1 ..\n"); exit(1);
	    if (buff[i] == 'N') {
	      for (k=0; k<4; k++) {
		v[ i*4+k ] += 0.25;
	      }
	      
	    }
	  } else {
	    v[ i * 4 + c ] += 1.0;
	  }
	  i ++;
	}		 
	
	cnt_lines ++;
      }
    }

    
  }

  
  //printf("Motif def has %d lines.\n", cnt_lines);
  
  
  //*m = v;
  
  
  *n = nummotifs;



}


void readMotif(char* file, float** m, int** st, int* l) 
{
  char* buff = 0;
  FILE* fp;
  int   maxlen = 1000;
  int   i, len;
  int   cnt_stars, cnt_lines;
  int   mw = 0;
  float* v;
  int*  myst = 0;
  int   c;
  int   k;

  buff = (char*)calloc(maxlen, sizeof(char));
  
  fp    = fopen(file, "r");
  if (!fp) {
    printf("Cannot open motif file %s. Exiting.\n", file);
    exit(0);
  }
  

  cnt_lines = 0;
  while (!feof(fp)) {

    fgets(buff, maxlen, fp);
    if (feof(fp)) {
      //printf("EOF, exit loop\n");
      break;
    }


    chomp(buff);
    len = strlen(buff);

    //printf("Reading %s\n", buff);
    
    if (buff[0] == 'M') {
      //printf(" Motif line\n");

    } else {
      
      cnt_stars = 0;
      for (i=0; i<len; i++) {
	if (buff[i] == '*') {
	  cnt_stars ++;
	}
      }
      


      // star line or nucleotide line ?
      if (cnt_stars > 0) {
	
	//printf(" Star line\n");
  
	
	if (mw <= 0) {
	  printf("pb .. mw=%d should be > 0 \n", mw);
	}
	

	myst = (int*)calloc(mw, sizeof(int));
	
	for (i=0; i<mw; i++) {

	  if (buff[i] == '*') {
	    myst[i] = 1;	  
	  } else {
	    myst[i] = 0;
	  }
	  
	}
	
	*st = myst;
	
      } else {
	
	if (cnt_lines == 0) {
	  // estimate the width of the motif
	  mw = 0;
	  while ((mw < len) && (buff[mw] != '\t') && (buff[mw] != ' ')) {
	    mw++; 
	  }
	  *l = mw;
	  v = (float*)calloc(mw * 4, sizeof(float));
	}
	
	i = 0;
	while ((i < len) && (buff[i] != '\t') && (buff[i] != ' ')) {
	  c = nucl( buff[i] );
	  if (c == -1) {
	    //printf("c should not be -1 ..\n"); exit(1);
	    if (buff[i] == 'N') {
	      for (k=0; k<4; k++) {
		v[ i*4+k ] += 0.25;
	      }
	      
	    }
	  } else {
	    v[ i * 4 + c ] += 1.0;
	  }
	  i ++;
	}		 
	
	cnt_lines ++;
      }
    }

    
  }

    
  //printf("Motif def has %d lines.\n", cnt_lines);
  
  // calculate frequencies here, using a pseudocount-like approach
  for (i=0; i<mw*4; i++) {
    v[i] = ( v[i] + 0.5 ) / (cnt_lines + 0.5);
    //v[i] = v[i] / cnt_lines;
  }  
  
  *m = v;
  
  
  
}
  

int nucl(char c) 
{
  
  if (c == 'A')
    return 0;
  else if (c == 'C')
    return 1;
  else if (c == 'G')
    return 2;
  else if (c == 'T')
    return 3;
  else return -1;
  
  
}


float pearson(float* data1, float* data2, int m) 
{
  
  float sum1   = 0.0;
  float sum2   = 0.0;
    
  float sum1_2 = 0.0;
  float sum2_2 = 0.0;

  float avg1   = 0.0;
  float avg2   = 0.0;
    
  float a, b, c;

  int i;

  for (i=0; i<m; i++) {

    sum1   += data1[i];
    sum2   += data2[i];
    
    sum1_2 += data1[i] * data1[i];
    sum2_2 += data2[i] * data2[i];
    
  }
    
  avg1 = sum1 / m;
  avg2 = sum2 / m;

  // calc the Pearson correlation
    
  a = 0.0;
  b = 0.0;
  c = 0.0;

  for (i=0; i<m; i++) {
    a += (data1[i] - avg1) * (data2[i] - avg2);
    b += (data1[i] - avg1) * (data1[i] - avg1);
    c += (data2[i] - avg2) * (data2[i] - avg2);
  }
  
    
  return a / sqrt ( b * c );

}

int* bubbleSortIndex(float *ia, int n)
{
  int i, j;

  float* a;
  int* b;


  // for storing temporary values
  int   itmp;
  float ftmp;

  a = (float*)malloc(n * sizeof(float));
  b = (int*)malloc(n * sizeof(int));

  for (i=0; i<n; i++) {
    a[i] = ia[i];
    b[i] = i;
  }

  for (i=0; i<n-1; i++) {
    for (j=0; j<n-1-i; j++)
      if (a[j+1] > a[j]) {  /* compare the two neighbors */

        ftmp = a[j];         /* swap a[j] and a[j+1]      */
        a[j] = a[j+1];
        a[j+1] = ftmp;

        itmp = b[j];         /* swap a[j] and a[j+1]      */
        b[j] = b[j+1];
        b[j+1] = itmp;

      }
  }

  free(a);

  return b;

}



int exist_flag(int argc, char** argv, char* param)
{
  int i = 0;
  while ((i < argc) && (strcmp(param, argv[i])))
    i++;
  if (i<argc)
    return 1;
  else
    return 0;
}
