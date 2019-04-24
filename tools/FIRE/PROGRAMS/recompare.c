// input: set of sequences containing the motif
//        regexp les plus "fit"
//
//

#include <stdio.h>
#ifdef BNS
#include <pcre/pcre.h>
#else
#include <pcre.h>
#endif

#include <string.h>
#include <math.h>
#include <search.h>

#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)

char* get_parameter(int argc, char** argv, char* param);
void get_overlap_array(int* site1, int* site2, int size, int nbcopies, int* i, int* s1, int* s2, char* site3); 

//int findSites(char* r, char* name, char* seq, int* positions, int* np); 
int findSites(char* r, char* name, char* seq, int* positions, int* np, int* orientations, int* no, int singlestrand); 


double lcumhyper(int i, int s1, int s2, int N); 
double cumhyper(int i, int s1, int s2, int N); 
double bico(int n, int k);
double factln(int n);
double gammln(double xx);
double factrl(int n) ;
void   nrerror(char str[]);
double hypergeom(int i, int s1, int s2, int N);

void chomp(char* s); 
char* complement(char* s); 
char *substr(char *string, int start, int length); 

int verbose = 1;


char* nextSequence(FILE* fp, char** name, int* size); 
char* nextSequence_currentLine; 
int   nextSequence_started = 0;
int   nextSequence_ended = 0;



int main(int argc, char** argv) {

  int   i, len, l=0, j;
  FILE* fp1, *fp2, *fp;
  char* buff;
  
  char* seqa = 0;
  char* seq = 0;
  char* tmp = 0;
  char* name;
  int size;
  int maxdist;

  int* positions;
  int* orientations;

  int np;
  int no;

  char* fasta1;
  char* fasta2;
  int nbgenes;
  char* re;

  int* sites_species1;
  int* sites_species2;
  char* sites_species1_species2;
  int s;
  int overlap_i;
  int overlap_s1;
  int overlap_s2;

  double score;


  char* outfile;
  char* posfile;
  char* orifile;

  char** names;

  int** sites_species1_positions;
  int** sites_species2_positions;

  int* sites_species1_nbpositions;
  int* sites_species2_nbpositions;


  int** sites_species1_orientations;
  int** sites_species2_orientations;

  int* sites_species1_nborientations;
  int* sites_species2_nborientations;

  

  int mu=-1, ts=-1, p1, l1;

  int nbcopies = 1;

  int add = 0;

  int singlestrand = 0;

  int   exclude = 0;
  char* excludefile;
  ENTRY e, *ep;
 

  if (argc < 4) {
    printf("Usage : recompare -fasta1 -fasta2 -nbgenes -re -out -pos -ori -mu -2s -nbcopies -maxdist\n");
    exit(0);
  }

  fasta1 = get_parameter(argc, argv, "-fasta1");
  fasta2 = get_parameter(argc, argv, "-fasta2");
  nbgenes      = atoi(get_parameter(argc, argv, "-nbgenes"));
  re         = get_parameter(argc, argv, "-re");

  outfile     = get_parameter(argc, argv, "-out");

  posfile     = get_parameter(argc, argv, "-pos");

  orifile     = get_parameter(argc, argv, "-ori");


  if (strcmp(get_parameter(argc, argv, "-mu"), "") != 0)
    mu          = atoi(get_parameter(argc, argv, "-mu"));

  if (strcmp(get_parameter(argc, argv, "-2s"), "") != 0)
    ts          = atoi(get_parameter(argc, argv, "-2s"));

  if (strcmp(get_parameter(argc, argv, "-nbcopies"), "") != 0)
    nbcopies         = atoi(get_parameter(argc, argv, "-nbcopies"));

 

  if (strcmp(get_parameter(argc, argv, "-maxdist"), "") != 0)
    maxdist      = atoi(get_parameter(argc, argv, "-maxdist"));
  else 
    maxdist      = 0;



  if (strcmp(get_parameter(argc, argv, "-rna"), "") != 0)
    singlestrand      = atoi(get_parameter(argc, argv, "-rna"));
  else 
    singlestrand      = 0;
  
  

  if (strcmp(get_parameter(argc, argv, "-exclude"), "") != 0) {
    exclude          = 1;
    excludefile      = get_parameter(argc, argv, "-exclude");
  } else 
    exclude = 0;


  positions = (int*)malloc(5000 * sizeof(int));

  orientations = (int*)malloc(5000 * sizeof(int));


  
  /** simple array of 0/1 representing  **/
  sites_species1 = (int*)calloc(nbgenes, sizeof(int));
  sites_species2 = (int*)calloc(nbgenes, sizeof(int));
  
  sites_species1_species2 = (char*)calloc(nbgenes, sizeof(char));


  sites_species1_positions   = (int**)malloc(nbgenes * sizeof(int*));
  sites_species2_positions   = (int**)malloc(nbgenes * sizeof(int*));
  sites_species1_nbpositions = (int*)malloc(nbgenes * sizeof(int));
  sites_species2_nbpositions = (int*)malloc(nbgenes * sizeof(int));
  
  
  sites_species1_orientations   = (int**)malloc(nbgenes * sizeof(int*));
  sites_species2_orientations   = (int**)malloc(nbgenes * sizeof(int*));
  sites_species1_nborientations = (int*)malloc(nbgenes * sizeof(int));
  sites_species2_nborientations = (int*)malloc(nbgenes * sizeof(int));

  
 

  //
  //  load the exclusion file
  //
  if (exclude == 1) {
    
    fp = fopen(excludefile, "r");
    if (!fp1) {
      printf("cannot open excludefile %s ..\n", excludefile);
      exit(0);
    }
    
    hcreate(1000000);
    buff = (char*)malloc(10000 * sizeof(char));

    while (!feof(fp)) {
      fscanf(fp, "%s\n", buff);
      e.key  = strdup(buff);
      e.data = (char *)(1);
      hsearch(e, ENTER);
    }

    fclose(fp);
  }

  names          = (char**)calloc(nbgenes, sizeof(char*));


  nextSequence_currentLine = (char*)malloc(50000 * sizeof(char));
  
  fp1 = fopen(fasta1, "r");
  if (!fp1) {
    printf("cannot open %s ..\n", fasta1);
    exit(0);
  }


  s = 0;
  while (seq = nextSequence(fp1, &name, &size)) { 
    //printf("needle %s in stack %s\n", re, seq);
    

    if (exclude == 1) {
      e.key  = name;
      ep = hsearch(e, FIND);
      if (ep) {
	continue;
	free(seq);
      }
    }

    

    names[s] = strdup(name);

    if (seq && (strlen(seq) > 0)) {

      np = 0;
      no = 0;

      if (maxdist > 0)

	seqa = seq + (strlen(seq) - maxdist);

      else if ((ts >= 0) && (mu >= 0)) {

	p1   = strlen(seq) - (mu+ts);
	p1   = max(p1, 0);
	l1   = min(strlen(seq) - p1,      ts*2);
	seqa = substr(seq, p1, l1);

      }
      else
	seqa = seq;

      findSites(re, name, seqa, positions, &np, orientations, &no, singlestrand);
      

      //
      //  store positions
      //
      sites_species1_nbpositions[s] = np;
      if (np > 0) {
	sites_species1[s] = np;
	sites_species1_positions[s] = (int*)malloc(np * sizeof(int));
	memcpy(sites_species1_positions[s], positions, np * sizeof(int));
      }

      //
      //  store orientations
      //
      sites_species1_nborientations[s] = no;
      if (no > 0) {
	sites_species1_orientations[s] = (int*)malloc(no * sizeof(int));
	memcpy(sites_species1_orientations[s], orientations, no * sizeof(int));
      }

      
    }
    free(seq);
    s++;
  }


  fp2 = fopen(fasta2, "r");
  if (!fp2) {
    printf("cannot open %s ..\n", fasta2);
    exit(0);
  }
  
  nextSequence_started = 0;
  nextSequence_ended = 0;
  
  s = 0;
  while (seq = nextSequence(fp2, &name, &size)) {   

    if (exclude == 1) {
      e.key  = name;
      ep = hsearch(e, FIND);
      if (ep) {
	continue;
	free(seq);
      }
    }

    if (seq && (strlen(seq) > 0)) {
      
      // reset the numbers of positions / orientation
      np = 0;
      no = 0;

      if (maxdist > 0)
	seqa = seq + (strlen(seq) - maxdist);
      else if ((ts >= 0) && (mu >= 0)) {

	p1   = strlen(seq) - (mu+ts);
	p1   = max(p1, 0);

	l1   = min(strlen(seq) - p1,      ts*2);
	
	

	seqa = substr(seq, p1, l1);
      }
      else
	seqa = seq;
      
      findSites(re, name, seqa, positions, &np, orientations, &no, singlestrand);

      sites_species2_nbpositions[s] = np;

      if (np > 0) {
	sites_species2[s] = np;

	sites_species2_positions[s] = (int*)malloc(np * sizeof(int));
	memcpy(sites_species2_positions[s], positions, np * sizeof(int));


	

	//printf("%s\t", name);
	//for (i=0; i<np; i++) {
	//  printf("%d ", positions[i]);
	//}
	//printf("\n");
      }
    }
    free(seq);
    s++;
  }

  if (exclude == 1) {
    nbgenes = s;
 
  }

  /** now look at the results ! **/
  get_overlap_array(sites_species1, sites_species2, nbgenes, nbcopies, &overlap_i, &overlap_s1, &overlap_s2, sites_species1_species2);
  
  

  score = lcumhyper(overlap_i, overlap_s1, overlap_s2, nbgenes);
  printf("species1=%d, species2=%d, overlap=%d, nb. genes=%d\n",  overlap_s1, overlap_s2, overlap_i, nbgenes);
  printf("Conservation score for pattern %s = %5.4f\n", re, -score); 

  if (strcmp(outfile, "") != 0) {
    
    fp = fopen(outfile, "w");
    if (!fp) {
      printf("Could not open %s ..\n" ,outfile);
    }
    for (i=0; i<nbgenes; i++) {
      
      if (sites_species1_species2[i] == 1)
	fprintf(fp, "%s\n", names[i]);
      
    }

    fclose(fp);
    
  } 


  //
  //  save the positions
  //
  if (strcmp(posfile, "") != 0) {
    
    fp = fopen(posfile, "w");
    if (!fp) {
      printf("Could not open %s ..\n" ,posfile);
    }
    for (i=0; i<nbgenes; i++) {
      
      if (sites_species1_species2[i] == 1) {
	for (j=0; j<sites_species1_nbpositions[i]; j++) {
	  if ((mu >= 0) && (ts >= 0)) 
	    add = mu - ts;
	  else 
	    add = 0;
	  fprintf(fp, "%d\n", add + sites_species1_positions[i][j]);
	}
      }
      
    }
    fclose(fp);
  } 


  //
  //  save the orientations
  //
  if (strcmp(orifile, "") != 0) {
    
    fp = fopen(orifile, "w");
    if (!fp) {
      printf("Could not open %s ..\n" ,orifile);
    }
    for (i=0; i<nbgenes; i++) {
      
      if (sites_species1_species2[i] == 1) {
	for (j=0; j<sites_species1_nborientations[i]; j++) {
	  fprintf(fp, "%d\n", sites_species1_orientations[i][j]);
	}
      }
      
    }
    fclose(fp);
  } 

  hdestroy();
  
    
  return 0;


}
  

/** REGEXP, returns 1 if a site was found, 0 otherwise **/
int findSites(char* r, char* name, char* seq, int* positions, int* np, int* orientations, int* no, int singlestrand) 
{
  
  int i, j; //, j, n;
  int   score1, score2;
  pcre* re;
  const char* error;
  int   erroffset;
  int   rc;
  int   ovector[30];
  char* cr = 0;
  char*  substring_start;
  int  substring_length;
  int startoffset = 0;
  int  pos;
  char* seq_pos;


  //
  //  store the position where a RE has been found, to avoid counting twice the palindromes
  //
  seq_pos = (char*)calloc(strlen(seq), sizeof(char));

  *np = 0;
  *no = 0;
  

  
  //printf("searching for %s in %s\n", r, seq);


  //
  // ONE STRAND
  //
  
  
  re = pcre_compile(r,
		    0,
		    &error,
		    &erroffset,
		    NULL);
  
  while ( (rc = pcre_exec(re,
		 NULL,
		 seq,
		 strlen(seq),
		 startoffset,
		 0,
		 ovector,
			  30) ) > 0) {
  

    substring_start     = seq         + ovector[0]; 
    substring_length    = ovector[1]  - ovector[0]; 
    pos                 = strlen(seq) - ovector[0] - 1;
    // save the position
    positions[ *np ]    = pos;
    // save the orientation
    orientations[ *no ] = 1;
    startoffset         = ovector[0]  + substring_length;
    seq_pos[ pos ]      = 1;

    

    (*np)++;
    (*no)++;

  }


  //return 0;

  if (singlestrand == 0) {

    // OTHER STRAND
    
    cr = complement(r);
    
    //printf("Complement of %s is %s\n", r, cr);
    
    startoffset = 0;
    
    re = pcre_compile(cr,
		      0,
		      &error,
		      &erroffset,
		      NULL);
    
    while ( (rc = pcre_exec(re,
			    NULL,
			    seq,
			    strlen(seq),
			    startoffset,
			    0,
			    ovector,
			    30)) > 0) {
      substring_start     = seq + ovector[0]; 
      substring_length    = ovector[1] - ovector[0]; 
      pos                 = strlen(seq) - ovector[0] - 1;

      //  if not a palindromic version
      if (seq_pos[ pos ] != 1) {
	
	// save the position
	positions[ *np ] = pos;
	
	
	
	(*np)++;
      }  
      
      // save the orientation in all cases
      orientations[ *no ] = -1;
      (*no)++;
      
      startoffset = ovector[0] + substring_length;
      
      
    }

  }

    
  if (cr != 0) 
    free(cr);
  free(seq_pos);
  pcre_free(re);
  return 0;
}


//
// take a char* and replace every nt by its complement
//
char* complement(char* s) 
{

  int l = strlen(s);
  char* c;
  char  d;
  int i;
  
  c = (char*)calloc(l+1, sizeof(char));

  for (i=l-1; i>=0; i--) {
    if (s[i] == 'A')
      d = 'T';
    else if (s[i] == 'T')
      d = 'A';
    else if (s[i] == 'G')
      d = 'C';
    else if (s[i] == 'C') 
      d = 'G';
    else if (s[i] == '[')
      d = ']';
    else if (s[i] == ']')
      d = '[';
    else
      d = s[i];
    
    strncat(c, &d, 1);
  }
	 
  return c;
  
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

char* nextSequence(FILE* fp, char** name, int* size) 
{
  
  int len = 20000;
  char* seq;
  

  // not yet started, get the first line
  if (nextSequence_started == 0) {
    fgets(nextSequence_currentLine, len, fp);  
    chomp(nextSequence_currentLine);      
    nextSequence_started = 1;
  }
  
  // if file is finished, return 0
  if (nextSequence_ended == 1) {
    return 0;
  }

  // if started, line should be filled with ">..."
  if (nextSequence_currentLine[0] == '>') {
    

    *name = strdup(nextSequence_currentLine + 1);
    
    // create a new line
    seq = (char*)calloc(100000, sizeof(char) ); 

    while (1) {
      
      fgets(nextSequence_currentLine, len, fp);    
      
      if (feof(fp)) {
	nextSequence_ended  = 1;
	return seq;
      }

      chomp(nextSequence_currentLine);      
      if (strlen(nextSequence_currentLine) == 0) { 
	continue;
      }
      
      if (nextSequence_currentLine[0] == '>') {	
	return seq;	
      } else {
	strcat(seq, nextSequence_currentLine);
      }
      
    }

  } else return 0;
  
}


double lcumhyper(int i, int s1, int s2, int N) {
  
  return log(cumhyper(i, s1, s2, N));

}


double cumhyper(int i, int s1, int s2, int N) 
{
  
  int min = (s1<s2?s1:s2);
  double prod = 0.0;
    int a;
    double tmp = 0.0;

    // attention, in R this is i+1 !
  for (a=i; a<=min; a++) {
    tmp = (double)hypergeom(a, s1, s2, N);
    //printf("a=%d, p(X=a)=%15.15f, log(p)=%5.4f\n", a, tmp, log(tmp));
    prod += (double)hypergeom(a, s1, s2, N);
  }
  
  return prod;
}



void nrerror(char str[]) {
  printf("%s\n", str);
}

double bico(int n, int k)                      
     //Returns the binomial coefficient nk as a doubleing-point number.
{
     double factln(int n);

     return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
     //The floor function cleans up roundoff error for smaller values of n and k.
}

double hypergeom(int i, int s1, int s2, int N) {
    double factln(int n);


    return exp(factln(s1) + factln(N - s1) + factln(s2) + factln(N - s2) - 
			 factln(i) - factln(N) - factln(s1 - i) - factln(s2 - i) 
      - factln(N - s1 - s2 + i));
    
}

double factln(int n)
     //Returns ln(n!).
{
  double gammln(double xx);
  void nrerror(char error_text[]);
  static double a[101];                                        
  //A static array is automatically initialized to zero.
  
  if (n < 0) nrerror("Negative factorial in routine factln");
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));                                    
  //In range of table.
  else return gammln(n+1.0);                                  
  //Out of range of table.
}


double gammln(double xx)
     //Returns the value ln[ (xx)] for xx > 0.
{
  //Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
  //accuracy is good enough.
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}



double factrl(int n) 
     //Returns the value n! as a oating-point number. 
{ 
  double gammln(double xx); 
  void nrerror(char error_text[]); 
  static int ntop=4; 
  static double a[33]={1.0,1.0,2.0,6.0,24.0}; 
  //Fill in table only as required. 
  int j; 
  
  if (n < 0) 
    nrerror("Negative factorial in routine factrl"); 
  if (n > 32) 
    return exp(gammln(n+1.0)); 
  // Larger value than size of table is required. Actually, this big a value is going to over ow on many computers, but no harm in trying. 
  while (ntop<n) { 
    // Fill in table up to desired value. 
    j=ntop++; 
    a[ntop]=a[j]*ntop; 
  } 
  
  return a[n]; 

}

char* get_parameter(int argc, char** argv, char* param)
{
  int i = 0;
  while ((i < argc) && (strcmp(param, argv[i])))
    i++;
  if (i<argc)
    return (argv[i+1]);
  else
    return "";
}


void get_overlap_array(int* site1, int* site2, int size, int nbcopies, int* i, int* s1, int* s2, char* sites_species1_species2) 
{

  int overlap_i = 0, overlap_s1 = 0, overlap_s2 = 0;
  int ii;
    
  
  for (ii=0; ii<size; ii++){
    if ( (site1[ii] >= nbcopies) && (site2[ii] >= nbcopies) ) {
      sites_species1_species2[ii] = 1;
      overlap_i++;
    } 
    
    if (site1[ii] >= nbcopies) {
      overlap_s1++;
    }
    
    if (site2[ii] >= nbcopies) {
      overlap_s2++;
    }
    
  }
  
  *i = overlap_i;
  *s1  = overlap_s1;
  *s2 = overlap_s2;


   
}


char *substr(char *string, int start, int length) 
{
  char *substring;
  
  substring = (char*)calloc(length+1, sizeof(char));

  memcpy(substring, string+start, length*sizeof(char));
  substring[length] = '\0';

  return substring;
}
