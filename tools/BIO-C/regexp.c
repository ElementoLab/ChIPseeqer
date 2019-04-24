
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>

#ifdef USEMPATROL
#include <mpatrol.h>
#endif

#ifdef BNS
#include <pcre/pcre.h>
#else
#include <pcre.h>
#endif

#include "dataio.h"
#include "regexp.h"





void findSites_micombine(char* r, char* seq, unsigned short* positions, int* np, unsigned short* orientations, int* no, unsigned short maxnbsites, int singlestrand, int store_pos, int store_ori) 
{
  
  pcre* re;
  const char* error;
  int   erroffset;
  int   rc;
  int   ovector[30];
  char* cr = 0;
  char* substring_start;
  int   substring_length;
  int   startoffset = 0;

  char* seq_pos;
  int   pos;

  //
  //  store the position where a RE has been found, to avoid counting twice the palindromes
  //
  seq_pos = (char*)calloc(strlen(seq), sizeof(char));

  *np = 0;

  if (no != 0) 
    *no = 0;
  

  //
  // FORWARD STRAND
  //
  re = pcre_compile(r,
		    0,
		    &error,
		    &erroffset,
		    NULL);

  startoffset = 0;
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
    pos                 = ovector[0];

    if (store_pos == 1) {
      positions   [ *np ] = (short)pos;
    }

    if (store_ori == 1) {
      orientations[ *no ] = 1;
    }

    startoffset    = pos + 1; //substring_length;
    seq_pos[ pos ] = 1;
    
    (*np)++;
    
    if (no != 0)
      (*no)++;
    
    if (*np == maxnbsites) {
      die("too many sites, not enough allocated memory ..");
    }	

  }

  pcre_free(re);


  //
  // REVERSE STRAND
  //
  if (singlestrand == 0) {

    
    cr = complement(r);
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
      
      substring_start  = seq + ovector[0]; 
      substring_length = ovector[1] - ovector[0]; 
      

      pos = ovector[0];
      
      //  proceed if it is not a palindromic version (in which case we don't want to count it again)
      if (seq_pos[ pos ] != 1) {
	
	if (store_pos == 1) {
	  positions[ *np ] = (short)pos;
	}

	(*np)++;
	
	if (*np == maxnbsites) {
	  die("too many sites, not enough allocated memory ..");
	}	
	
      } 
      
      // save the orientation in all cases
      if (store_ori == 1) {
	orientations[ *no ] = -1;
	(*no)++;
      }

      
      
      
      startoffset = ovector[0] + 1; //substring_length;
      
    
    }
    
    pcre_free(re);
    free(cr);
  }
  
 
  free(seq_pos);
  

}



/** REGEXP, returns 1 if a site was found, 0 otherwise **/
int findAllSites(char* r, char* seq, int* positions, int* np, int* orientations, int* no, int singlestrand) 
{
  
  //int i, j; //, j, n;
  //int   score1, score2;
  pcre* re;
  const char* error;
  int   erroffset;
  int   rc;
  int   ovector[30];
  char* cr = 0;
  char*  substring_start;
  int  substring_length;
  int startoffset = 0;

  char* seq_pos;
  int pos;

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
  

    substring_start    = seq         + ovector[0]; 
    substring_length   = ovector[1]  - ovector[0]; 

    // save the position
    pos                 = ovector[0];

    positions   [ *np ] = pos;
    orientations[ *no ] = 1;    

    startoffset         = pos + substring_length;
    
    // store pos
    seq_pos[ pos ]      = 1;

    (*np)++;
    (*no)++;

  }

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

      substring_start  = seq + ovector[0]; 
      substring_length = ovector[1] - ovector[0]; 

      pos = ovector[0];

      if (seq_pos[ pos ] != 1) { 

	positions   [ *np ] = pos;
	orientations[ *no ] = -1;

	(*np)++;
	(*no)++;
      }

      startoffset        = pos + substring_length;
    }
    
  }
    
  if (cr != 0) 
    free(cr);
  free(seq_pos);
  pcre_free(re);
  return 0;
}





void findSites(char* r, char* seq, int* positions, int* np, int* orientations, int* no, int maxnbsites, int singlestrand, int store_pos_ori) 
{
  
  pcre* re;
  const char* error;
  int   erroffset;
  int   rc;
  int   ovector[30];
  char* cr = 0;
  char*  substring_start;
  int  substring_length;
  int startoffset = 0;

  char* seq_pos;
  int pos;

  //
  //  store the position where a RE has been found, to avoid counting twice the palindromes
  //
  seq_pos = (char*)calloc(strlen(seq), sizeof(char));

  *np = 0;
  *no = 0;
  

  //
  // FORWARD STRAND
  //
  re = pcre_compile(r,
		    0,
		    &error,
		    &erroffset,
		    NULL);

  startoffset = 0;
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
    pos                 = ovector[0];

    if (store_pos_ori == 1) {
      positions   [ *np ] = pos;
      orientations[ *no ] = 1;
    }

    startoffset         = pos + 1; //substring_length;

    seq_pos[ pos ] = 1;

    (*np)++;
    (*no)++;
    
    if (*np == maxnbsites) {
      die("too many sites, not enough allocated memory ..");
    }	

  }

  pcre_free(re);


  //
  // REVERSE STRAND
  //
  if (singlestrand == 0) {

    
    cr = complement(r);
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
      
      substring_start  = seq + ovector[0]; 
      substring_length = ovector[1] - ovector[0]; 
      

      pos = ovector[0];
      
      //  proceed if it is not a palindromic version (in which case we don't want to count it again)
      if (seq_pos[ pos ] != 1) {
	
	if (store_pos_ori == 1) {
	  positions[ *np ] = pos;
	}

	(*np)++;
	
	if (*np == maxnbsites) {
	  die("too many sites, not enough allocated memory ..");
	}	
	
      } 
      
      // save the orientation in all cases
      if (store_pos_ori == 1) {
	orientations[ *no ] = -1;
      }
      (*no)++;
      
      startoffset = ovector[0] + 1; //substring_length;
      
    
    }
    
    pcre_free(re);
    free(cr);
  }
  
 
  free(seq_pos);
  

}




int re_matches(char* r, char* seq, int singlestrand) 
{
  
  pcre* re;
  const char* error;
  int   erroffset;
  int   rc;
  int   ovector[30];
  char* cr = 0;
  int   startoffset = 0;
  
  //
  // FORWARD STRAND
  //
  re = pcre_compile(r,
		    0,
		    &error,
		    &erroffset,
		    NULL);

  startoffset = 0;
  if ((rc = pcre_exec(re,
		      NULL,
		      seq,
		      strlen(seq),
		      startoffset,
		      0,
		      ovector,
		      30) ) > 0) {
    return 1;
  }
  pcre_free(re);


  //
  // REVERSE STRAND
  //
  if (singlestrand != 1) {

    
    cr = complement(r);
    startoffset = 0;
    
    re = pcre_compile(cr,
		      0,
		      &error,
		      &erroffset,
		      NULL);
    
    if ((rc = pcre_exec(re,
			NULL,
			seq,
			strlen(seq),
			startoffset,
			0,
			ovector,
			30)) > 0) {
      return 1;
    }
    
    free(cr);
    pcre_free(re);
  }
  

  
  
  return 0;
}



int raw_re_matches(pcre* re, char* seq, char* seq_c, int singlestrand) 
{
  

  //const char* error;
  //int   erroffset;
  int   rc;
  int   ovector[30];
  //char* cr = 0;
  int   startoffset = 0;
  

  startoffset = 0;
  if ((rc = pcre_exec(re,
		      NULL,
		      seq,
		      strlen(seq),
		      startoffset,
		      0,
		      ovector,
		      30) ) > 0) {
    return 1;
  }
  //pcre_free(re);


  //
  // REVERSE STRAND
  //
  if (singlestrand != 1) {

    
    startoffset = 0;
        
    if ((rc = pcre_exec(re,
			NULL,
			seq_c,
			strlen(seq_c),
			startoffset,
			0,
			ovector,
			30)) > 0) {
      return 1;
    }
  }
  
  return 0;
}


