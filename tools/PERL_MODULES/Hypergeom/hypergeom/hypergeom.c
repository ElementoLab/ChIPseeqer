#include <math.h>
#include <stdio.h>
#include "hypergeom.h"

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

double bino(int x, int N, double p) 
{

  
  
  double lbico(int n, int x);                      
  
  return exp( lbico(N, x) + log(pow(p, x)) + log(pow(1.0-p, N-x)) );
  
}


double cumbino(int k, int N, double p) 
{

  int x;
  double pb;
  double lbico(int n, int k);                      

  pb = 0;
  for (x=k; x<=N; x++) {
    pb += exp( lbico(N, x) + log(pow(p, x)) + log(pow(1.0-p, N-x)) );
    
  }
  
  return pb;
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

double lbico(int n, int k)                      
     //Returns the binomial coefficient nk as a doubleing-point number.
{
     double factln(int n);

     return factln(n)-factln(k)-factln(n-k);
     //The floor function cleans up roundoff error for smaller values of n and k.
}

double hypergeom(int i, int s1, int s2, int N) {
    double factln(int n);

    //printf("N=%d, s1=%d, s2=%d, i=%d\n", N, s1, s2, i);

    /*
    printf("factln(s1)=%f\n", floor(0.5+exp(factln(s1))));
    printf("factln(N - s1)=%f\n", floor(0.5+exp(factln(N - s1))));
    printf("factln(s2)=%f\n", floor(0.5+exp(factln(s2))));
    printf("factln(N - s2)=%f\n", floor(0.5+exp(factln(N - s2))));
    printf("factln(N)=%f\n", floor(0.5+exp(factln(N))));
    printf("factln(i)=%f\n", floor(0.5+exp(factln(i))));
    printf("factln(s1 - i)=%f\n", floor(0.5+exp(factln(s1 - i))));
    printf("factln(s2 - i)=%f\n", floor(0.5+exp(factln(s2 - i))));
    printf("factln(N - s1 - s2 + i)=%f\n", floor(0.5+exp(factln(N - s1 - s2 + i))));
    */
    /**printf("sum=%f\n", factln(s1) + factln(N - s1) + factln(s2) + factln(N - s2) - 
			 factln(i) - factln(N) - factln(s1 - i) - factln(s2 - i) 
			 - factln(N - s1 - s2 + i));
    **/
	   /*
    return exp(factln(s1) + factln(N - s1) + factln(s2) + factln(N - s2) - 
			 factln(i) - factln(N) - factln(s1 - i) - factln(s2 - i) 
			 - factln(N - s1 - s2 + i)
			 );
	   */
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
  
  if (n < 0) printf("Negative factorial in routine factln %d\n", n);
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

