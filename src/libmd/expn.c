/*							md_expn.c
 *
 *		Exponential integral En
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, md_expn();
 *
 * y = md_expn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                 inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                  1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        5000       2.0e-16     4.6e-17
 *    IEEE      0, 30       10000       1.7e-15     3.6e-16
 *
 */

/*							md_expn.c	*/

/* Cephes Math Library Release 2.8:  June, 2000
   Copyright 1985, 2000 by Stephen L. Moshier */

#include "mconf.h"
#ifdef ANSIPROT
extern double md_pow ( double, double );
extern double md_gamma ( double );
extern double md_log ( double );
extern double md_exp ( double );
extern double md_fabs ( double );
#else
double md_pow(), md_gamma(), md_log(), md_exp(), md_fabs();
#endif
#define EUL 0.57721566490153286060
#define BIG  1.44115188075855872E+17
extern double MAXNUM, MACHEP, MAXLOG;

double md_expn( n, x )
int n;
double x;
{
double ans, r, t, yk, xk;
double pk, pkm1, pkm2, qk, qkm1, qkm2;
double psi, z;
int i, k;
static double big = BIG;

if( n < 0 )
	goto domerr;

if( x < 0 )
	{
domerr:	mtherr( "md_expn", DOMAIN );
	return( MAXNUM );
	}

if( x > MAXLOG )
	return( 0.0 );

if( x == 0.0 )
	{
	if( n < 2 )
		{
		mtherr( "md_expn", SING );
		return( MAXNUM );
		}
	else
		return( 1.0/(n-1.0) );
	}

if( n == 0 )
	return( md_exp(-x)/x );

/*							md_expn.c	*/
/*		Expansion for large n		*/

if( n > 5000 )
	{
	xk = x + n;
	yk = 1.0 / (xk * xk);
	t = n;
	ans = yk * t * (6.0 * x * x  -  8.0 * t * x  +  t * t);
	ans = yk * (ans + t * (t  -  2.0 * x));
	ans = yk * (ans + t);
	ans = (ans + 1.0) * md_exp( -x ) / xk;
	goto done;
	}

if( x > 1.0 )
	goto cfrac;

/*							md_expn.c	*/

/*		Power series expansion		*/

psi = -EUL - md_log(x);
for( i=1; i<n; i++ )
	psi = psi + 1.0/i;

z = -x;
xk = 0.0;
yk = 1.0;
pk = 1.0 - n;
if( n == 1 )
	ans = 0.0;
else
	ans = 1.0/pk;
do
	{
	xk += 1.0;
	yk *= z/xk;
	pk += 1.0;
	if( pk != 0.0 )
		{
		ans += yk/pk;
		}
	if( ans != 0.0 )
		t = md_fabs(yk/ans);
	else
		t = 1.0;
	}
while( t > MACHEP );
k = xk;
t = n;
r = n - 1;
ans = (md_pow(z, r) * psi / md_gamma(t)) - ans;
goto done;

/*							md_expn.c	*/
/*		continued fraction		*/
cfrac:
k = 1;
pkm2 = 1.0;
qkm2 = x;
pkm1 = 1.0;
qkm1 = x + n;
ans = pkm1/qkm1;

do
	{
	k += 1;
	if( k & 1 )
		{
		yk = 1.0;
		xk = n + (k-1)/2;
		}
	else
		{
		yk = x;
		xk = k/2;
		}
	pk = pkm1 * yk  +  pkm2 * xk;
	qk = qkm1 * yk  +  qkm2 * xk;
	if( qk != 0 )
		{
		r = pk/qk;
		t = md_fabs( (ans - r)/r );
		ans = r;
		}
	else
		t = 1.0;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
if( md_fabs(pk) > big )
		{
		pkm2 /= big;
		pkm1 /= big;
		qkm2 /= big;
		qkm1 /= big;
		}
	}
while( t > MACHEP );

ans *= md_exp( -x );

done:
return( ans );
}

