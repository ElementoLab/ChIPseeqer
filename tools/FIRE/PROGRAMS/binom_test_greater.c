#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "beta.h"


int main(int argc, char** argv) {

  //printf("%f\n", pbinom(9.0, 31.0, 0.046, 1, 0));

  int d;
  double p;
  double p2;

  if (argc < 4) {
    printf("Usage: prg k n p 0/1\n");
    exit(0);
  }

  //printf("'%s'\n", argv[4]);

  if ((atoi(argv[1]) == 0) && (atoi(argv[1]) == 0)) {
    if (d == 0)
      printf("%e\n", 1.0);
    else 
      printf("%e\t%e\n", 1.0, 1.0);
  }
  
  if (argv[4] == 0)
    d = 0;
  else 
    d = atoi(argv[4]);

  p = pbinom( atof(argv[1]) - 1, atof(argv[2]), atof(argv[3]), 0, 0);


  if (d == 0) {
    printf("%e\n", p);
  } else if (d == 1) {
    // follows the R code for binom.test
    //    PVAL <- switch(alternative, less = pbinom(x, n, p), greater = pbinom(x - 
    //                                                 1, n, p, lower = FALSE),
    p2 = pbinom( atof(argv[1]) , atof(argv[2]), atof(argv[3]), 1, 0);
    printf("%e\t%e\n", p, p2); 
  } else {
    printf("don't know what you mean, please indicate 0 or 1 for the 4th argument\n");
  }

  return 0;
}
