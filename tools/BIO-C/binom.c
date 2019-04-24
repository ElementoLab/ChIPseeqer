#define _GNU_SOURCE
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rcode.h"
#include "binom.h"
#include "statistics.h"





#define min(a,b) (((a) < (b))?(a):(b))
#define max(a,b) (((a) > (b))?(a):(b))

static double bfrac(double, double, double, double, double, double);
static void bgrat(double, double, double, double, double *, double, int *);
static void grat1(double, double, double, double *, double *, double);
static double apser(double, double, double, double);
static double bpser(double, double, double, double);
static double basym(double, double, double, double);
static double fpser(double, double, double, double);
static double bup(double, double, double, double, int, double);
static double exparg(int);
static double psi(double);
static double gam1(double);
static double gamln1(double);
static double betaln(double, double);
static double algdiv(double, double);
static double brcmp1(int, double, double, double, double);
static double brcomp(double, double, double, double);
static double rlog1(double);
static double bcorr(double, double);
static double gamln(double);
static double alnrel(double);
static double esum(int, double);
static double erf__(double);
static double rexp(double);
static double erfc1(int, double);
static double gsumln(double, double);

double Rf_d1mach(int i);
int Rf_i1mach(int i);

double lgammafn(double x);


double cumbinom(int k, int n, double p)
{
  if (k == 0)
    return 1.0;
  else 
    return pbinom( (double)k - 1, (double)n, p, 0, 0);
}



double fmin2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? x : y;
}



double fmax2(double x, double y)
{
#ifdef IEEE_754
        if (ISNAN(x) || ISNAN(y))
                return x + y;
#endif
        return (x < y) ? y : x;
}


double stirlerr(double n)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static double sferr_halves[31] = {
	0.0, /* n=0 - wrong, place holder only */
	0.1534264097200273452913848,  /* 0.5 */
	0.0810614667953272582196702,  /* 1.0 */
	0.0548141210519176538961390,  /* 1.5 */
	0.0413406959554092940938221,  /* 2.0 */
	0.03316287351993628748511048, /* 2.5 */
	0.02767792568499833914878929, /* 3.0 */
	0.02374616365629749597132920, /* 3.5 */
	0.02079067210376509311152277, /* 4.0 */
	0.01848845053267318523077934, /* 4.5 */
	0.01664469118982119216319487, /* 5.0 */
	0.01513497322191737887351255, /* 5.5 */
	0.01387612882307074799874573, /* 6.0 */
	0.01281046524292022692424986, /* 6.5 */
	0.01189670994589177009505572, /* 7.0 */
	0.01110455975820691732662991, /* 7.5 */
	0.010411265261972096497478567, /* 8.0 */
	0.009799416126158803298389475, /* 8.5 */
	0.009255462182712732917728637, /* 9.0 */
	0.008768700134139385462952823, /* 9.5 */
	0.008330563433362871256469318, /* 10.0 */
	0.007934114564314020547248100, /* 10.5 */
	0.007573675487951840794972024, /* 11.0 */
	0.007244554301320383179543912, /* 11.5 */
	0.006942840107209529865664152, /* 12.0 */
	0.006665247032707682442354394, /* 12.5 */
	0.006408994188004207068439631, /* 13.0 */
	0.006171712263039457647532867, /* 13.5 */
	0.005951370112758847735624416, /* 14.0 */
	0.005746216513010115682023589, /* 14.5 */
	0.005554733551962801371038690  /* 15.0 */
    };
    double nn;

    if (n <= 15.0) {
	nn = n + n;
	if (nn == (int)nn) return(sferr_halves[(int)nn]);
	return(lgammafn(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
    }

    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

double log1p(double x)
{
    /* series for log1p on the interval -.375 to .375
     *				     with weighted error   6.35e-32
     *				      log weighted error  31.20
     *			    significant figures required  30.93
     *				 decimal places required  32.01
     */
    const static double alnrcs[43] = {
	+.10378693562743769800686267719098e+1,
	-.13364301504908918098766041553133e+0,
	+.19408249135520563357926199374750e-1,
	-.30107551127535777690376537776592e-2,
	+.48694614797154850090456366509137e-3,
	-.81054881893175356066809943008622e-4,
	+.13778847799559524782938251496059e-4,
	-.23802210894358970251369992914935e-5,
	+.41640416213865183476391859901989e-6,
	-.73595828378075994984266837031998e-7,
	+.13117611876241674949152294345011e-7,
	-.23546709317742425136696092330175e-8,
	+.42522773276034997775638052962567e-9,
	-.77190894134840796826108107493300e-10,
	+.14075746481359069909215356472191e-10,
	-.25769072058024680627537078627584e-11,
	+.47342406666294421849154395005938e-12,
	-.87249012674742641745301263292675e-13,
	+.16124614902740551465739833119115e-13,
	-.29875652015665773006710792416815e-14,
	+.55480701209082887983041321697279e-15,
	-.10324619158271569595141333961932e-15,
	+.19250239203049851177878503244868e-16,
	-.35955073465265150011189707844266e-17,
	+.67264542537876857892194574226773e-18,
	-.12602624168735219252082425637546e-18,
	+.23644884408606210044916158955519e-19,
	-.44419377050807936898878389179733e-20,
	+.83546594464034259016241293994666e-21,
	-.15731559416479562574899253521066e-21,
	+.29653128740247422686154369706666e-22,
	-.55949583481815947292156013226666e-23,
	+.10566354268835681048187284138666e-23,
	-.19972483680670204548314999466666e-24,
	+.37782977818839361421049855999999e-25,
	-.71531586889081740345038165333333e-26,
	+.13552488463674213646502024533333e-26,
	-.25694673048487567430079829333333e-27,
	+.48747756066216949076459519999999e-28,
	-.92542112530849715321132373333333e-29,
	+.17578597841760239233269760000000e-29,
	-.33410026677731010351377066666666e-30,
	+.63533936180236187354180266666666e-31,
    };
    static double xmin; /*was sqrt(d1mach(4)); */

#ifdef NOMORE_FOR_THREADS
    static int nlnrel = 0;

    if (nlnrel == 0) {/* initialize chebychev coefficients */
	nlnrel = chebyshev_init(alnrcs, 43, DBL_EPSILON/20);/*was .1*d1mach(3)*/
    }
#else
# define nlnrel 22
/* 22: for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */
#endif
     
    xmin = -1 + sqrt(1/DBL_EPSILON);
    
    if (x == 0.) return 0.;/* speed */
    if (x == -1) return(ML_NEGINF);
    if (x  < -1) ML_ERR_return_NAN;

    if (fabs(x) <= .375) {
        /* Improve on speed (only);
	   again give result accurate to IEEE double precision: */
	if(fabs(x) < .5 * DBL_EPSILON)
	    return x;

	if( (0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
	    return x * (1 - .5 * x);
	/* else */
	return x * (1 - x * chebyshev_eval(x / .375, alnrcs, nlnrel));
    }
    /* else */
    if (x < xmin) {
	/* answer less than half precision because x too near -1 */
	ML_ERROR(ME_PRECISION);
    }
    return log(1 + x);
}





double expm1(double x)
{
    double y, a = fabs(x);

    if (a < DBL_EPSILON) return x;
    if (a > 0.697) return exp(x) - 1;  /* negligible cancellation */

    if (a > 1e-8)
        y = exp(x) - 1;
    else /* Taylor expansion, more accurate in this range */
        y = (x / 2 + 1) * x;

    /* Newton step for solving   log(1 + y) = x   for y : */
    /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
    y -= (1 + y) * (log1p (y) - x);
    return y;
}




double chebyshev_eval(double x, const double *a, const int n)
{
    double b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) ML_ERR_return_NAN;

    if (x < -1.1 || x > 1.1) ML_ERR_return_NAN;

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}




double gammafn(double x)
{
    const static double gamcs[42] = {
	+.8571195590989331421920062399942e-2,
	+.4415381324841006757191315771652e-2,
	+.5685043681599363378632664588789e-1,
	-.4219835396418560501012500186624e-2,
	+.1326808181212460220584006796352e-2,
	-.1893024529798880432523947023886e-3,
	+.3606925327441245256578082217225e-4,
	-.6056761904460864218485548290365e-5,
	+.1055829546302283344731823509093e-5,
	-.1811967365542384048291855891166e-6,
	+.3117724964715322277790254593169e-7,
	-.5354219639019687140874081024347e-8,
	+.9193275519859588946887786825940e-9,
	-.1577941280288339761767423273953e-9,
	+.2707980622934954543266540433089e-10,
	-.4646818653825730144081661058933e-11,
	+.7973350192007419656460767175359e-12,
	-.1368078209830916025799499172309e-12,
	+.2347319486563800657233471771688e-13,
	-.4027432614949066932766570534699e-14,
	+.6910051747372100912138336975257e-15,
	-.1185584500221992907052387126192e-15,
	+.2034148542496373955201026051932e-16,
	-.3490054341717405849274012949108e-17,
	+.5987993856485305567135051066026e-18,
	-.1027378057872228074490069778431e-18,
	+.1762702816060529824942759660748e-19,
	-.3024320653735306260958772112042e-20,
	+.5188914660218397839717833550506e-21,
	-.8902770842456576692449251601066e-22,
	+.1527474068493342602274596891306e-22,
	-.2620731256187362900257328332799e-23,
	+.4496464047830538670331046570666e-24,
	-.7714712731336877911703901525333e-25,
	+.1323635453126044036486572714666e-25,
	-.2270999412942928816702313813333e-26,
	+.3896418998003991449320816639999e-27,
	-.6685198115125953327792127999999e-28,
	+.1146998663140024384347613866666e-28,
	-.1967938586345134677295103999999e-29,
	+.3376448816585338090334890666666e-30,
	-.5793070335782135784625493333333e-31
    };

    int i, n;
    double y;
    double sinpiy, value;

#ifdef NOMORE_FOR_THREADS
    static int ngam = 0;
    static double xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

    /* Initialize machine dependent constants, the first time gamma() is called.
	FIXME for threads ! */
    if (ngam == 0) {
	ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
	gammalims(&xmin, &xmax);/*-> ./gammalims.c */
	xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
	/*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
	dxrel = sqrt(1/DBL_EPSILON);/*was (1/d1mach(4)) */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 * (xmin, xmax) are non-trivial, see ./gammalims.c
 * xsml = exp(.01)*DBL_MIN
 * dxrel = sqrt(1/DBL_EPSILON) = 2 ^ 26
*/
# define ngam 22
# define xmin -170.5674972726612
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 67108864.
#endif

    if(ISNAN(x)) return x;

    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x == 0 || (x < 0 && x == (long)x)) {
	ML_ERROR(ME_RANGE);
	return ML_NAN;
    }

    y = fabs(x);

    if (y <= 10) {

	/* Compute gamma(x) for -10 <= x <= 10
	 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	 * first of all. */

	n = x;
	if(x < 0) --n;
	y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
	--n;
	value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
	if (n == 0)
	    return value;/* x = 1.dddd = 1+y */

	if (n < 0) {
	    /* compute gamma(x) for -10 <= x < 1 */

	    /* exact 0 or "-n" checked already above */

	    /* The answer is less than half precision */
	    /* because x too near a negative integer. */
	    if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
		ML_ERROR(ME_PRECISION);
	    }

	    /* The argument is so close to 0 that the result would overflow. */
	    if (y < xsml) {
		ML_ERROR(ME_RANGE);
		if(x > 0) return ML_POSINF;
		else return ML_NEGINF;
	    }

	    n = -n;

	    for (i = 0; i < n; i++) {
		value /= (x + i);
	    }
	    return value;
	}
	else {
	    /* gamma(x) for 2 <= x <= 10 */

	    for (i = 1; i <= n; i++) {
		value *= (y + i);
	    }
	    return value;
	}
    }
    else {
	/* gamma(x) for	 y = |x| > 10. */

	if (x > xmax) {			/* Overflow */
	    ML_ERROR(ME_RANGE);
	    return ML_POSINF;
	}

	if (x < xmin) {			/* Underflow */
	    ML_ERROR(ME_UNDERFLOW);
	    return ML_UNDERFLOW;
	}

	if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
	    value = 1.;
	    for (i = 2; i < y; i++) value *= i;
	}
	else { /* normal case */
	    value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
			((2*y == (int)2*y)? stirlerr(y) : lgammacor(y)));
	}
	if (x > 0)
	    return value;

	if (fabs((x - (int)(x - 0.5))/x) < dxrel){

	    /* The answer is less than half precision because */
	    /* the argument is too near a negative integer. */

	    ML_ERROR(ME_PRECISION);
	}

	sinpiy = sin(M_PI * y);
	if (sinpiy == 0) {		/* Negative integer arg - overflow */
	    ML_ERROR(ME_RANGE);
	    return ML_POSINF;
	}

	return -M_PI / (y * sinpiy * value);
    }
}


int signgam;

double lgammafn(double x)
{
    double ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
    static double xmax = 0.;
    static double dxrel = 0.;

    if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
	xmax = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
	dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
   xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
   dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
 */
#define xmax  2.5327372760800758e+305
#define dxrel 1.490116119384765696e-8
#endif

    signgam = 1;

#ifdef IEEE_754
    if(ISNAN(x)) return x;
#endif

    if (x < 0 && fmod(floor(-x), 2.) == 0)
	signgam = -1;

#ifdef HAVE_TRUNC
    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
#else
    if (x <= 0 && x == rint(x)) { /* Negative integer argument */      
#endif
	ML_ERROR(ME_RANGE);
	return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = fabs(x);

    if (y <= 10)
	return log(fabs(gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax) {
	ML_ERROR(ME_RANGE);
	return ML_POSINF;
    }

    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
	if(x > 1e17)
	    return(x*(log(x) - 1.));
	else if(x > 4934720.)
	    return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
	else
#endif
	    return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sin(M_PI * y));

    if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
	MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
	ML_ERR_return_NAN;
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

#ifdef HAVE_TRUNC
    if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

	/* The answer is less than half precision because
	 * the argument is too near a negative integer. */

	ML_ERROR(ME_PRECISION);
    }
#endif

    return ans;
}



double lgammacor(double x)
{
    const static double algmcs[15] = {
	+.1666389480451863247205729650822e+0,
	-.1384948176067563840732986059135e-4,
	+.9810825646924729426157171547487e-8,
	-.1809129475572494194263306266719e-10,
	+.6221098041892605227126015543416e-13,
	-.3399615005417721944303330599666e-15,
	+.2683181998482698748957538846666e-17,
	-.2868042435334643284144622399999e-19,
	+.3962837061046434803679306666666e-21,
	-.6831888753985766870111999999999e-23,
	+.1429227355942498147573333333333e-24,
	-.3547598158101070547199999999999e-26,
	+.1025680058010470912000000000000e-27,
	-.3401102254316748799999999999999e-29,
	+.1276642195630062933333333333333e-30
    };

    double tmp;

#ifdef NOMORE_FOR_THREADS
    static int nalgm = 0;
    static double xbig = 0, xmax = 0;

    /* Initialize machine dependent constants, the first time gamma() is called.
	FIXME for threads ! */
    if (nalgm == 0) {
	/* For IEEE double precision : nalgm = 5 */
	nalgm = chebyshev_init(algmcs, 15, DBL_EPSILON/2);/*was d1mach(3)*/
	xbig = 1 / sqrt(DBL_EPSILON/2); /* ~ 94906265.6 for IEEE double */
	xmax = exp(fmin2(log(DBL_MAX / 12), -log(12 * DBL_MIN)));
	/*   = DBL_MAX / 48 ~= 3.745e306 for IEEE double */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xbig = 2 ^ 26.5
 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
# define nalgm 5
# define xbig  94906265.62425156
# define xmax  3.745194030963158e306
#endif

    if (x < 10)
	ML_ERR_return_NAN
    else if (x >= xmax) {
	ML_ERROR(ME_UNDERFLOW);
	return ML_UNDERFLOW;
    }
    else if (x < xbig) {
	tmp = 10 / x;
	return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    else return 1 / (x * 12);
}


double lbeta(double a, double b)
{
    double corr, p, q;

    p = q = a;
    if(b < p) p = b;/* := min(a,b) */
    if(b > q) q = b;/* := max(a,b) */

#ifdef IEEE_754
    if(ISNAN(a) || ISNAN(b))
	return a + b;
#endif

    /* both arguments must be >= 0 */

    if (p < 0)
	ML_ERR_return_NAN
    else if (p == 0) {
	return ML_POSINF;
    }
    else if (!R_FINITE(q)) {
	return ML_NEGINF;
    }

    if (p >= 10) {
	/* p and q are big. */
	corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
	return log(q) * -0.5 + M_LN_SQRT_2PI + corr
		+ (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
    }
    else if (q >= 10) {
	/* p is small, but q is big. */
	corr = lgammacor(q) - lgammacor(p + q);
	return lgammafn(p) + corr + p - p * log(p + q)
		+ (q - 0.5) * log1p(-p / (p + q));
    }
    else
	/* p and q are small: p <= q < 10. */
	return log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
}


double pbeta_raw(double x, double pin, double qin, int lower_tail)
{
    double ans, c, finsum, p, ps, p1, q, term, xb, xi, y;
    int n, i, ib, swap_tail;

    const static double eps = .5*DBL_EPSILON;
    const static double sml = DBL_MIN;
    const double lneps = log(eps);
    const double lnsml = log(sml);

    /* Switch to TOMS 708 if p or q is large */
    if (pin > 15 || qin > 15) {
	double x1 = 1 - x, w, wc;
	int ierr;
	bratio(pin, qin, x, x1, &w, &wc, &ierr);
	if(ierr)
	  printf("pbeta_raw() -> bratio() gave error code %d", ierr);
	return lower_tail ? w : wc;
    }
    
    /* swap tails if x is greater than the mean */
    if (pin / (pin + qin) < x) {
	swap_tail = 1;
	y = 1 - x;
	p = qin;
	q = pin;
    }
    else {
	swap_tail = 0;
	y = x;
	p = pin;
	q = qin;
    }

    if ((p + q) * y / (p + 1) < eps) {

	/* tail approximation */

	xb = p * log(fmax2(y, sml)) - log(p) - lbeta(p, q);
	if (xb > lnsml && y != 0) {
	    ans = (swap_tail == lower_tail) ? -expm1(xb) : exp(xb);
	} else {
	    ans = (swap_tail == lower_tail) ? 1. : 0;
	}
    } else {
	/* evaluate the infinite sum first.  term will equal */
	/* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */

	/* Ly := log(y) */
	double Ly = swap_tail ? log1p(-x) : log(y);

	ps = q - floor(q);
	xb = p * Ly;
	if (ps == 0)
	    ps = 1; /*==> lbeta(ps,p)= log Beta(1,p) = log(1/p) = -log(p) */
	else
	    xb -= (lbeta(ps, p) + log(p));
	ans = 0;
	if (xb >= lnsml) {
	    ans = exp(xb);
	    term = ans * p;
	    if (ps != 1) {
		n = fmax2(lneps/Ly, 4.0);
		for(i=1 ; i <= n ; i++) {
		    xi = i;
		    term *= (xi - ps) * y / xi;
		    ans += term / (p + xi);
		}
	    }
	}

	/* now evaluate the finite sum, maybe. */

	if (q > 1) {

	    double liy;/* == log(1-y) */
	    if(swap_tail) {
		c = 1./x;/* == 1/(1 - y) */
		liy = log(x);
	    }
	    else {
		c = 1./(1. - y);
		liy = log1p(-y);
	    }
	    xb = p * Ly + q * liy - lbeta(p, q) - log(q);
	    ib = fmax2(xb / lnsml, 0.0);
	    term = exp(xb - ib * lnsml);
	    p1 = q * c / (p + q - 1);

	    finsum = 0;
	    n = q;
	    if (q == n)
		n--;
	    for(i= 1; i <= n; i++) {
		if (p1 <= 1 && term / eps <= finsum)
		    break;
		xi = i;
		term = (q - xi + 1) * c * term / (p + q - xi);
		if (term > 1) {
		    ib--;
		    term *= sml;
		}
		if (ib == 0)
		    finsum += term;
	    }
	    ans += finsum;
	}
	if (swap_tail == lower_tail)
	    ans = 1 - ans;
	ans = fmax2(fmin2(ans, 1.), 0.);
    }
    return ans;
} /* pbeta_raw() */

double pbeta(double x, double pin, double qin, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(pin) || ISNAN(qin))
	return x + pin + qin;
#endif

    if (pin <= 0 || qin <= 0) ML_ERR_return_NAN;

    if (x <= 0)
	return R_DT_0;
    if (x >= 1)
	return R_DT_1;
    /* FIXME: this doesn't allow an extended value range when log_p = TRUE */
    return R_D_val(pbeta_raw(x, pin, qin, lower_tail));
}



double pbinom(double x, double n, double p, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n) || ISNAN(p))
	return x + n + p;
    if (!R_FINITE(n) || !R_FINITE(p)) ML_ERR_return_NAN;

#endif
    if(R_D_nonint(n)) ML_ERR_return_NAN;
    n = R_D_forceint(n);
    if(n <= 0 || p < 0 || p > 1) ML_ERR_return_NAN;

    x = floor(x + 1e-7);
    if (x < 0.0) return R_DT_0;
    if (n <= x) return R_DT_1;
    return pbeta(p, x + 1, n - x, !lower_tail, log_p);
}
/* Based on C translation of ACM TOMS 708
   Please do not change this, e.g. to use R's versions of the 
   ancillary routines, without investigating the error analysis as we
   do need very high relative accuracy.  This version has about
   14 digits accuracy.
*/


/*      ALGORITHM 708, COLLECTED ALGORITHMS FROM ACM. */
/*      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 360-373z. */
void bratio(double a, double b, double x, double y, double *w,
	    double *w1, int *ierr)
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    int n;
    double t, z__, a0, b0, x0, y0;
    int ind;
    double eps;
    int ierr1;
    double lambda;

/* ----------------------------------------------------------------------- */

/*            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B) */

/*                     -------------------- */

/*     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1 */
/*     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES */

/*                      W  = IX(A,B) */
/*                      W1 = 1 - IX(A,B) */

/*     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. */
/*     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND */
/*     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED, */
/*     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO */
/*     ONE OF THE FOLLOWING VALUES ... */

/*        IERR = 1  IF A OR B IS NEGATIVE */
/*        IERR = 2  IF A = B = 0 */
/*        IERR = 3  IF X .LT. 0 OR X .GT. 1 */
/*        IERR = 4  IF Y .LT. 0 OR Y .GT. 1 */
/*        IERR = 5  IF X + Y .NE. 1 */
/*        IERR = 6  IF X = A = 0 */
/*        IERR = 7  IF Y = B = 0 */

/* -------------------- */
/*     WRITTEN BY ALFRED H. MORRIS, JR. */
/*        NAVAL SURFACE WARFARE CENTER */
/*        DAHLGREN, VIRGINIA */
/*     REVISED ... NOV 1991 */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST */
/*            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0 */

    eps = 2.0 * Rf_d1mach(3);

/* ----------------------------------------------------------------------- */
    *w = 0.0;
    *w1 = 0.0;
    if (a < 0.0 || b < 0.0) {
	goto L300;
    }
    if (a == 0.0 && b == 0.0) {
	goto L310;
    }
    if (x < 0.0 || x > 1.0) {
	goto L320;
    }
    if (y < 0.0 || y > 1.0) {
	goto L330;
    }
    z__ = x + y - 0.5 - 0.5;
    if (fabs(z__) > eps * 3.0) {
	goto L340;
    }

    *ierr = 0;
    if (x == 0.0) {
	goto L200;
    }
    if (y == 0.0) {
	goto L210;
    }
    if (a == 0.0) {
	goto L211;
    }
    if (b == 0.0) {
	goto L201;
    }

    eps = max(eps,1e-15);
    if (max(a,b) < eps * .001) {
	goto L230;
    }

    ind = 0;
    a0 = a;
    b0 = b;
    x0 = x;
    y0 = y;
    if (min(a0,b0) > 1.0) {
	goto L30;
    }

/*             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1 */

    if (x <= 0.5) {
	goto L10;
    }
    ind = 1;
    a0 = b;
    b0 = a;
    x0 = y;
    y0 = x;

L10:
/* Computing MIN */
    if (b0 < min(eps, eps * a0)) {
	goto L80;
    }
/* Computing MIN */
    if (a0 < min(eps, eps * b0) && b0 * x0 <= 1.0) {
	goto L90;
    }
    if (max(a0,b0) > 1.0) {
	goto L20;
    }
    if (a0 >= min(.2,b0)) {
	goto L100;
    }
    if (pow(x0, a0) <= .9f) {
	goto L100;
    }
    if (x0 >= .3f) {
	goto L110;
    }
    n = 20;
    goto L130;

L20:
    if (b0 <= 1.0) {
	goto L100;
    }
    if (x0 >= .3f) {
	goto L110;
    }
    if (x0 >= .1f) {
	goto L21;
    }
    if (pow(x0*b0, a0) <= .7f) {
	goto L100;
    }
L21:
    if (b0 > 15.0) {
	goto L131;
    }
    n = 20;
    goto L130;

/*             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1 */

L30:
    if (a > b) {
	goto L31;
    }
    lambda = a - (a + b) * x;
    goto L32;
L31:
    lambda = (a + b) * y - b;
L32:
    if (lambda >= 0.0) {
	goto L40;
    }
    ind = 1;
    a0 = b;
    b0 = a;
    x0 = y;
    y0 = x;
    lambda = fabs(lambda);

L40:
    if (b0 < 40.0 && b0 * x0 <= .7f) {
	goto L100;
    }
    if (b0 < 40.0) {
	goto L140;
    }
    if (a0 > b0) {
	goto L50;
    }
    if (a0 <= 100.0) {
	goto L120;
    }
    if (lambda > a0 * .03f) {
	goto L120;
    }
    goto L180;
L50:
    if (b0 <= 100.0) {
	goto L120;
    }
    if (lambda > b0 * .03f) {
	goto L120;
    }
    goto L180;

/*            EVALUATION OF THE APPROPRIATE ALGORITHM */

L80:
    *w = fpser(a0, b0, x0, eps);
    *w1 = 0.5 - *w + 0.5;
    goto L220;

L90:
    *w1 = apser(a0, b0, x0, eps);
    *w = 0.5 - *w1 + 0.5;
    goto L220;

L100:
    *w = bpser(a0, b0, x0, eps);
    *w1 = 0.5 - *w + 0.5;
    goto L220;

L110:
    *w1 = bpser(b0, a0, y0, eps);
    *w = 0.5 - *w1 + 0.5;
    goto L220;

L120:
    d__1 = eps * 15.0;
    *w = bfrac(a0, b0, x0, y0, lambda, d__1);
    *w1 = 0.5 - *w + 0.5;
    goto L220;

L130:
    *w1 = bup(b0, a0, y0, x0, n, eps);
    b0 += n;
L131:
    bgrat(b0, a0, y0, x0, w1, 15*eps, &ierr1);
    *w = 0.5 - *w1 + 0.5;
    goto L220;

L140:
    n = (int) b0;
    b0 -= n;
    if (b0 != 0.0) {
	goto L141;
    }
    --n;
    b0 = 1.0;
L141:
    *w = bup(b0, a0, y0, x0, n, eps);
    if (x0 > .7f) {
	goto L150;
    }
    *w += bpser(a0, b0, x0, eps);
    *w1 = 0.5 - *w + 0.5;
    goto L220;

L150:
    if (a0 > 15.0) {
	goto L151;
    }
    n = 20;
    *w += bup(a0, b0, x0, y0, n, eps);
    a0 += n;
L151:
    bgrat(a0, b0, x0, y0, w, 15*eps, &ierr1);
    *w1 = 0.5 - *w + 0.5;
    goto L220;

L180:
    d__1 = eps * 100.0;
    *w = basym(a0, b0, lambda, d__1);
    *w1 = 0.5 - *w + 0.5;
    goto L220;

/*               TERMINATION OF THE PROCEDURE */

L200:
    if (a == 0.0) {
	goto L350;
    }
L201:
    *w = 0.0;
    *w1 = 1.0;
    return;

L210:
    if (b == 0.0) {
	goto L360;
    }
L211:
    *w = 1.0;
    *w1 = 0.0;
    return;

L220:
    if (ind == 0) {
	return;
    }
    t = *w;
    *w = *w1;
    *w1 = t;
    return;

/*           PROCEDURE FOR A AND B .LT. 1.D-3*EPS */

L230:
    *w = b / (a + b);
    *w1 = a / (a + b);
    return;

/*                       ERROR RETURN */

L300:
    *ierr = 1;
    return;
L310:
    *ierr = 2;
    return;
L320:
    *ierr = 3;
    return;
L330:
    *ierr = 4;
    return;
L340:
    *ierr = 5;
    return;
L350:
    *ierr = 6;
    return;
L360:
    *ierr = 7;
    return;
} /* bratio */

double fpser(double a, double b, double x, double eps)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    double c__, s, t, an, tol;

/* ----------------------------------------------------------------------- */

/*                 EVALUATION OF I (A,B) */
/*                                X */

/*          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5. */

/* ----------------------------------------------------------------------- */

/*                  SET  FPSER = X**A */

    ret_val = 1.0;
    if (a <= eps * .001f) {
	goto L10;
    }
    ret_val = 0.0;
    t = a * log(x);
    if (t < exparg(1)) {
	return ret_val;
    }
    ret_val = exp(t);

/*                NOTE THAT 1/B(A,B) = B */

L10:
    ret_val = b / a * ret_val;
    tol = eps / a;
    an = a + 1.0;
    t = x;
    s = t / an;
L20:
    an += 1.0;
    t = x * t;
    c__ = t / an;
    s += c__;
    if (fabs(c__) > tol) {
	goto L20;
    }

    ret_val *= a * s + 1.0;
    return ret_val;
} /* fpser */

static double apser(double a, double b, double x, double eps)
{
    /* Initialized data */

    static double g = .577215664901533;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double c__, j, s, t, aj, bx;
    double tol;

/* ----------------------------------------------------------------------- */
/*     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR */
/*     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN */
/*     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED. */
/* ----------------------------------------------------------------------- */
/* -------------------- */
/* -------------------- */
    bx = b * x;
    t = x - bx;
    if (b * eps > .02f) {
	goto L10;
    }
    c__ = log(x) + psi(b) + g + t;
    goto L20;
L10:
    c__ = log(bx) + g + t;

L20:
    tol = eps * 5.0 * fabs(c__);
    j = 1.0;
    s = 0.0;
L30:
    j += 1.0;
    t *= x - bx / j;
    aj = t / j;
    s += aj;
    if (fabs(aj) > tol) {
	goto L30;
    }

    ret_val = -(a) * (c__ + s);
    return ret_val;
} /* apser */

static double bpser(double a, double b, double x, double eps)
{
    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    double c__;
    int i__, m;
    double n, t, u, w, z__, a0, b0, apb, tol, sum;

/* ----------------------------------------------------------------------- */
/*     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1 */
/*     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED. */
/* ----------------------------------------------------------------------- */

    ret_val = 0.0;
    if (x == 0.0) {
	return ret_val;
    }
/* ----------------------------------------------------------------------- */
/*            COMPUTE THE FACTOR X**A/(A*BETA(A,B)) */
/* ----------------------------------------------------------------------- */
    a0 = min(a,b);
    if (a0 < 1.0) {
	goto L10;
    }
    z__ = a * log(x) - betaln(a, b);
    ret_val = exp(z__) / a;
    goto L70;
L10:
    b0 = max(a,b);
    if (b0 >= 8.0) {
	goto L60;
    }
    if (b0 > 1.0) {
	goto L40;
    }

/*            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1 */

    ret_val = pow(x, a);
    if (ret_val == 0.0) {
	return ret_val;
    }

    apb = a + b;
    if (apb > 1.0) {
	goto L20;
    }
    z__ = gam1(apb) + 1.0;
    goto L30;
L20:
    u = (double) (a) + (double) (b) - 1.;
    z__ = (gam1(u) + 1.0) / apb;

L30:
    c__ = (gam1(a) + 1.0) * (gam1(b) + 1.0) / z__;
    ret_val = ret_val * c__ * (b / apb);
    goto L70;

/*         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8 */

L40:
    u = gamln1(a0);
    m = b0 - 1.0;
    if (m < 1) {
	goto L50;
    }
    c__ = 1.0;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b0 += -1.0;
/* L41: */
	c__ *= b0 / (a0 + b0);
    }
    u = log(c__) + u;

L50:
    z__ = a * log(x) - u;
    b0 += -1.0;
    apb = a0 + b0;
    if (apb > 1.0) {
	goto L51;
    }
    t = gam1(apb) + 1.0;
    goto L52;
L51:
    u = a0 + b0 - 1.;
    t = (gam1(u) + 1.0) / apb;
L52:
    ret_val = exp(z__) * (a0 / a) * (gam1(b0) + 1.0) / t;
    goto L70;

/*            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8 */

L60:
    u = gamln1(a0) + algdiv(a0, b0);
    z__ = a * log(x) - u;
    ret_val = a0 / a * exp(z__);
L70:
    if (ret_val == 0.0 || a <= eps * .1f) {
	return ret_val;
    }
/* ----------------------------------------------------------------------- */
/*                     COMPUTE THE SERIES */
/* ----------------------------------------------------------------------- */
    sum = 0.0;
    n = 0.0;
    c__ = 1.0;
    tol = eps / a;
L100:
    n += 1.0;
    c__ = c__ * (0.5 - b / n + 0.5) * x;
    w = c__ / (a + n);
    sum += w;
    if (fabs(w) > tol) {
	goto L100;
    }
    ret_val *= a * sum + 1.0;
    return ret_val;
} /* bpser */

static double bup(double a, double b, double x, double y, int n, double eps)
{
    /* System generated locals */
    int i__1;
    double ret_val, d__1;

    /* Local variables */
    double d__;
    int i__, k;
    double l, r__, t, w;
    int mu;
    double ap1;
    int nm1, kp1;
    double apb;

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INT. */
/*     EPS IS THE TOLERANCE USED. */
/* ----------------------------------------------------------------------- */

/*          OBTAIN THE SCALING FACTOR EXP(-MU) AND */
/*             EXP(MU)*(X**A*Y**B/BETA(A,B))/A */

    apb = a + b;
    ap1 = a + 1.0;
    mu = 0;
    d__ = 1.0;
    if (n == 1 || a < 1.0) {
	goto L10;
    }
    if (apb < ap1 * 1.1f) {
	goto L10;
    }
    mu = (d__1 = exparg(1), (int) fabs(d__1));
    k = (int) exparg(0);
    if (k < mu) {
	mu = k;
    }
    t = (double) mu;
    d__ = exp(-t);

L10:
    ret_val = brcmp1(mu, a, b, x, y) / a;
    if (n == 1 || ret_val == 0.0) {
	return ret_val;
    }
    nm1 = n - 1;
    w = d__;

/*          LET K BE THE INDEX OF THE MAXIMUM TERM */

    k = 0;
    if (b <= 1.0) {
	goto L40;
    }
    if (y > 1e-4) {
	goto L20;
    }
    k = nm1;
    goto L30;
L20:
    r__ = (b - 1.0) * x / y - a;
    if (r__ < 1.0) {
	goto L40;
    }
    k = nm1;
    t = (double) nm1;
    if (r__ < t) {
	k = (int) r__;
    }

/*          ADD THE INCREASING TERMS OF THE SERIES */

L30:
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = (double) (i__ - 1);
	d__ = (apb + l) / (ap1 + l) * x * d__;
	w += d__;
/* L31: */
    }
    if (k == nm1) {
	goto L50;
    }

/*          ADD THE REMAINING TERMS OF THE SERIES */

L40:
    kp1 = k + 1;
    i__1 = nm1;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	l = (double) (i__ - 1);
	d__ = (apb + l) / (ap1 + l) * x * d__;
	w += d__;
	if (d__ <= eps * w) {
	    goto L50;
	}
/* L41: */
    }

/*               TERMINATE THE PROCEDURE */

L50:
    ret_val *= w;
    return ret_val;
} /* bup */

static double bfrac(double a, double b, double x, double y, double lambda, double eps)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double c__, e, n, p, r__, s, t, w, c0, c1, r0, an, bn, yp1, anp1, bnp1,
	    beta, alpha;

/* ----------------------------------------------------------------------- */
/*     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1. */
/*     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B. */
/* ----------------------------------------------------------------------- */
/* -------------------- */
    ret_val = brcomp(a, b, x, y);
    if (ret_val == 0.0) {
	return ret_val;
    }

    c__ = lambda + 1.0;
    c0 = b / a;
    c1 = 1.0 / a + 1.0;
    yp1 = y + 1.0;

    n = 0.0;
    p = 1.0;
    s = a + 1.0;
    an = 0.0;
    bn = 1.0;
    anp1 = 1.0;
    bnp1 = c__ / c1;
    r__ = c1 / c__;

/*        CONTINUED FRACTION CALCULATION */

L10:
    n += 1.0;
    t = n / a;
    w = n * (b - n) * x;
    e = a / s;
    alpha = p * (p + c0) * e * e * (w * x);
    e = (t + 1.0) / (c1 + t + t);
    beta = n + w / s + e * (c__ + n * yp1);
    p = t + 1.0;
    s += 2.0;

/*        UPDATE AN, BN, ANP1, AND BNP1 */

    t = alpha * an + beta * anp1;
    an = anp1;
    anp1 = t;
    t = alpha * bn + beta * bnp1;
    bn = bnp1;
    bnp1 = t;

    r0 = r__;
    r__ = anp1 / bnp1;
    if ((d__1 = r__ - r0, fabs(d__1)) <= eps * r__) {
	goto L20;
    }

/*        RESCALE AN, BN, ANP1, AND BNP1 */

    an /= bnp1;
    bn /= bnp1;
    anp1 = r__;
    bnp1 = 1.0;
    goto L10;

/*                 TERMINATION */

L20:
    ret_val *= r__;
    return ret_val;
} /* bfrac */

static double brcomp(double a, double b, double x, double y)
{
    /* Initialized data */

    static double const__ = .398942280401433;

    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    double c__, e, h__;
    int i__, n;
    double t, u, v, z__, a0, b0, x0, y0, apb, lnx, lny;
    double lambda;

/* ----------------------------------------------------------------------- */
/*               EVALUATION OF X**A*Y**B/BETA(A,B) */
/* ----------------------------------------------------------------------- */
/* ----------------- */
/*     CONST = 1/SQRT(2*PI) */
/* ----------------- */

    ret_val = 0.0;
    if (x == 0.0 || y == 0.0) {
	return ret_val;
    }
    a0 = min(a, b);
    if (a0 >= 8.0) {
	goto L100;
    }

    if (x > .375) {
	goto L10;
    }
    lnx = log(x);
    lny = alnrel(-x);
    goto L20;
L10:
    if (y > .375) {
	goto L11;
    }
    lnx = alnrel(-y);
    lny = log(y);
    goto L20;
L11:
    lnx = log(x);
    lny = log(y);

L20:
    z__ = a * lnx + b * lny;
    if (a0 < 1.0) {
	goto L30;
    }
    z__ -= betaln(a, b);
    ret_val = exp(z__);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .LT. 1 OR B .LT. 1 */
/* ----------------------------------------------------------------------- */
L30:
    b0 = max(a, b);
    if (b0 >= 8.0) {
	goto L80;
    }
    if (b0 > 1.0) {
	goto L60;
    }

/*                   ALGORITHM FOR B0 .LE. 1 */

    ret_val = exp(z__);
    if (ret_val == 0.0) {
	return ret_val;
    }

    apb = a + b;
    if (apb > 1.0) {
	goto L40;
    }
    z__ = gam1(apb) + 1.0;
    goto L50;
L40:
    u = a + b - 1.;
    z__ = (gam1(u) + 1.0) / apb;

L50:
    c__ = (gam1(a) + 1.0) * (gam1(b) + 1.0) / z__;
    ret_val = ret_val * (a0 * c__) / (a0 / b0 + 1.0);
    return ret_val;

/*                ALGORITHM FOR 1 .LT. B0 .LT. 8 */

L60:
    u = gamln1(a0);
    n = b0 - 1.0;
    if (n < 1) {
	goto L70;
    }
    c__ = 1.0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b0 += -1.0;
	c__ *= b0 / (a0 + b0);
/* L61: */
    }
    u = log(c__) + u;

L70:
    z__ -= u;
    b0 += -1.0;
    apb = a0 + b0;
    if (apb > 1.0) {
	goto L71;
    }
    t = gam1(apb) + 1.0;
    goto L72;
L71:
    u = a0 + b0 - 1.;
    t = (gam1(u) + 1.0) / apb;
L72:
    ret_val = a0 * exp(z__) * (gam1(b0) + 1.0) / t;
    return ret_val;

/*                   ALGORITHM FOR B0 .GE. 8 */

L80:
    u = gamln1(a0) + algdiv(a0, b0);
    ret_val = a0 * exp(z__ - u);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .GE. 8 AND B .GE. 8 */
/* ----------------------------------------------------------------------- */
L100:
    if (a > b) {
	goto L101;
    }
    h__ = a / b;
    x0 = h__ / (h__ + 1.0);
    y0 = 1.0 / (h__ + 1.0);
    lambda = a - (a + b) * x;
    goto L110;
L101:
    h__ = b / a;
    x0 = 1.0 / (h__ + 1.0);
    y0 = h__ / (h__ + 1.0);
    lambda = (a + b) * y - b;

L110:
    e = -lambda / a;
    if (fabs(e) > .6) {
	goto L111;
    }
    u = rlog1(e);
    goto L120;
L111:
    u = e - log(x / x0);

L120:
    e = lambda / b;
    if (fabs(e) > .6) {
	goto L121;
    }
    v = rlog1(e);
    goto L130;
L121:
    v = e - log(y / y0);

L130:
    z__ = exp(-(a * u + b * v));
    ret_val = const__ * sqrt(b * x0) * z__ * exp(-bcorr(a, b));
    return ret_val;
} /* brcomp */

static double brcmp1(int mu, double a, double b, double x, double y)
{
    /* Initialized data */

    static double const__ = .398942280401433;

    /* System generated locals */
    int i__1;
    double ret_val, r__1;

    /* Local variables */
    double c__, e, h__;
    int i__, n;
    double t, u, v, z__, a0, b0, x0, y0, apb, lnx, lny;
    double lambda;

/* ----------------------------------------------------------------------- */
/*          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B)) */
/* ----------------------------------------------------------------------- */
/* ----------------- */
/*     CONST = 1/SQRT(2*PI) */
/* ----------------- */

    a0 = min(a,b);
    if (a0 >= 8.0) {
	goto L100;
    }

    if (x > .375) {
	goto L10;
    }
    lnx = log(x);
    lny = alnrel(-x);
    goto L20;
L10:
    if (y > .375) {
	goto L11;
    }
    lnx = alnrel(-y);
    lny = log(y);
    goto L20;
L11:
    lnx = log(x);
    lny = log(y);

L20:
    z__ = a * lnx + b * lny;
    if (a0 < 1.0) {
	goto L30;
    }
    z__ -= betaln(a, b);
    ret_val = esum(mu, z__);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .LT. 1 OR B .LT. 1 */
/* ----------------------------------------------------------------------- */
L30:
    b0 = max(a,b);
    if (b0 >= 8.0) {
	goto L80;
    }
    if (b0 > 1.0) {
	goto L60;
    }

/*                   ALGORITHM FOR B0 .LE. 1 */

    ret_val = esum(mu, z__);
    if (ret_val == 0.0) {
	return ret_val;
    }

    apb = a + b;
    if (apb > 1.0) {
	goto L40;
    }
    z__ = gam1(apb) + 1.0;
    goto L50;
L40:
    u = a + b - 1.;
    z__ = (gam1(u) + 1.0) / apb;

L50:
    c__ = (gam1(a) + 1.0) * (gam1(b) + 1.0) / z__;
    ret_val = ret_val * (a0 * c__) / (a0 / b0 + 1.0);
    return ret_val;

/*                ALGORITHM FOR 1 .LT. B0 .LT. 8 */

L60:
    u = gamln1(a0);
    n = b0 - 1.0;
    if (n < 1) {
	goto L70;
    }
    c__ = 1.0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b0 += -1.0;
	c__ *= b0 / (a0 + b0);
/* L61: */
    }
    u = log(c__) + u;

L70:
    z__ -= u;
    b0 += -1.0;
    apb = a0 + b0;
    if (apb > 1.0) {
	goto L71;
    }
    t = gam1(apb) + 1.0;
    goto L72;
L71:
    u = a0 + b0 - 1.;
    t = (gam1(u) + 1.0) / apb;
L72:
    ret_val = a0 * esum(mu, z__) * (gam1(b0) + 1.0) / t;
    return ret_val;

/*                   ALGORITHM FOR B0 .GE. 8 */

L80:
    u = gamln1(a0) + algdiv(a0, b0);
    ret_val = a0 * esum(mu, z__ - u);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .GE. 8 AND B .GE. 8 */
/* ----------------------------------------------------------------------- */
L100:
    if (a > b) {
	goto L101;
    }
    h__ = a / b;
    x0 = h__ / (h__ + 1.0);
    y0 = 1.0 / (h__ + 1.0);
    lambda = a - (a + b) * x;
    goto L110;
L101:
    h__ = b / a;
    x0 = 1.0 / (h__ + 1.0);
    y0 = h__ / (h__ + 1.0);
    lambda = (a + b) * y - b;

L110:
    e = -lambda / a;
    if (fabs(e) > .6f) {
	goto L111;
    }
    u = rlog1(e);
    goto L120;
L111:
    u = e - log(x / x0);

L120:
    e = lambda / b;
    if (fabs(e) > .6f) {
	goto L121;
    }
    v = rlog1(e);
    goto L130;
L121:
    v = e - log(y / y0);

L130:
    r__1 = -(a * u + b * v);
    z__ = esum(mu, r__1);
    ret_val = const__ * sqrt(b * x0) * z__ * exp(-bcorr(a, b));
    return ret_val;
} /* brcmp1 */

static void bgrat(double a, double b, double x, double y, double *w,
		  double eps, int *ierr)
{
    /* System generated locals */
    int i__1;
    double r__1;

    /* Local variables */
    double c__[30], d__[30];
    int i__;
    double j, l;
    int n;
    double p, q, r__, s, t, u, v, z__, n2, t2, dj, cn, nu, bm1;
    int nm1;
    double lnx, sum;
    double bp2n, coef;

/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B. */
/*     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED */
/*     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED. */
/*     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. */
/* ----------------------------------------------------------------------- */

    bm1 = b - 0.5 - 0.5;
    nu = a + bm1 * 0.5;
    if (y > .375f) {
	goto L10;
    }
    lnx = alnrel(-y);
    goto L11;
L10:
    lnx = log(x);
L11:
    z__ = -nu * lnx;
    if (b * z__ == 0.0) {
	goto L100;
    }

/*                 COMPUTATION OF THE EXPANSION */
/*                 SET R = EXP(-Z)*Z**B/GAMMA(B) */

    r__ = b * (gam1(b) + 1.0) * exp(b * log(z__));
    r__ = r__ * exp(a * lnx) * exp(bm1 * 0.5 * lnx);
    u = algdiv(b, a) + b * log(nu);
    u = r__ * exp(-u);
    if (u == 0.0) {
	goto L100;
    }
    grat1(b, z__, r__, &p, &q, eps);

/* Computing 2nd power */
    r__1 = 1.0 / nu;
    v = r__1 * r__1 * .25;
    t2 = lnx * .25f * lnx;
    l = *w / u;
    j = q / r__;
    sum = j;
    t = 1.0;
    cn = 1.0;
    n2 = 0.0;
    for (n = 1; n <= 30; ++n) {
	bp2n = b + n2;
	j = (bp2n * (bp2n + 1.0) * j + (z__ + bp2n + 1.0) * t) * v;
	n2 += 2.0;
	t *= t2;
	cn /= n2 * (n2 + 1.0);
	c__[n - 1] = cn;
	s = 0.0;
	if (n == 1) {
	    goto L21;
	}
	nm1 = n - 1;
	coef = b - n;
	i__1 = nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s += coef * c__[i__ - 1] * d__[n - i__ - 1];
/* L20: */
	    coef += b;
	}
L21:
	d__[n - 1] = bm1 * cn + s / n;
	dj = d__[n - 1] * j;
	sum += dj;
	if (sum <= 0.0) {
	    goto L100;
	}
	if (fabs(dj) <= eps * (sum + l)) {
	    goto L30;
	}
/* L22: */
    }

/*                    ADD THE RESULTS TO W */

L30:
    *ierr = 0;
    *w += u * sum;
    return;

/*               THE EXPANSION CANNOT BE COMPUTED */

L100:
    *ierr = 1;
    return;
} /* bgrat */

static void grat1(double a, double x, double r, double *p, double *q,
		  double eps)
{
    /* Local variables */
    double c__, g, h__, j, l, t, w, z__, an, am0, an0, a2n, b2n, cma;
    double tol, sum;
    double a2nm1, b2nm1;

/* ----------------------------------------------------------------------- */
/*        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS */
/*                      P(A,X) AND Q(A,X) */

/*     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED. */
/*     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A). */
/* ----------------------------------------------------------------------- */
    if (a * x == 0.0) {
	goto L130;
    }
    if (a == 0.5) {
	goto L120;
    }
    if (x < 1.1f) {
	goto L10;
    }
    goto L50;

/*             TAYLOR SERIES FOR P(A,X)/X**A */

L10:
    an = 3.0;
    c__ = x;
    sum = x / (a + 3.0);
    tol = eps * .1f / (a + 1.0);
L11:
    an += 1.0;
    c__ = -c__ * (x / an);
    t = c__ / (a + an);
    sum += t;
    if (fabs(t) > tol) {
	goto L11;
    }
    j = a * x * ((sum / 6.0 - 0.5 / (a + 2.0)) * x + 1.0 / (a + 1.0));

    z__ = a * log(x);
    h__ = gam1(a);
    g = h__ + 1.0;
    if (x < .25f) {
	goto L20;
    }
    if (a < x / 2.59f) {
	goto L40;
    }
    goto L30;
L20:
    if (z__ > -.13394f) {
	goto L40;
    }

L30:
    w = exp(z__);
    *p = w * g * (0.5 - j + 0.5);
    *q = 0.5 - *p + 0.5;
    return;

L40:
    l = rexp(z__);
    w = l + 0.5 + 0.5;
    *q = (w * j - l) * g - h__;
    if (*q < 0.0) {
	goto L110;
    }
    *p = 0.5 - *q + 0.5;
    return;

/*              CONTINUED FRACTION EXPANSION */

L50:
    a2nm1 = 1.0;
    a2n = 1.0;
    b2nm1 = x;
    b2n = x + (1.0 - a);
    c__ = 1.0;
L51:
    a2nm1 = x * a2n + c__ * a2nm1;
    b2nm1 = x * b2n + c__ * b2nm1;
    am0 = a2nm1 / b2nm1;
    c__ += 1.0;
    cma = c__ - a;
    a2n = a2nm1 + cma * a2n;
    b2n = b2nm1 + cma * b2n;
    an0 = a2n / b2n;
    if (fabs(an0 - am0) >= eps * an0) {
	goto L51;
    }
    *q = r * an0;
    *p = 0.5 - *q + 0.5;
    return;

/*                SPECIAL CASES */

L100:
    *p = 0.0;
    *q = 1.0;
    return;

L110:
    *p = 1.0;
    *q = 0.0;
    return;

L120:
    if (x >= .25f) {
	goto L121;
    }
    *p = erf__(sqrt(x));
    *q = 0.5 - *p + 0.5;
    return;
L121:
    *q = erfc1(0, sqrt(x));
    *p = 0.5 - *q + 0.5;
    return;

L130:
    if (x <= a) {
	goto L100;
    }
    goto L110;
} /* grat1 */

static double basym(double a, double b, double lambda, double eps)
{
    /* Initialized data */

    static int num = 20;
    static double e0 = 1.12837916709551;
    static double e1 = .353553390593274;

    /* System generated locals */
    int i__1, i__2, i__3, i__4;
    double ret_val, r__1, r__2;

    /* Local variables */
    double c__[21], d__[21], f, h__;
    int i__, j, m, n;
    double r__, s, t, u, w, z__, a0[21], b0[21], j0, j1, h2, r0, r1, t0, t1, w0,
	     z0, z2, hn, zn;
    int im1, mm1, np1, imj, mmj;
    double sum, znm1, bsum, dsum;

/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B. */
/*     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED. */
/*     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT */
/*     A AND B ARE GREATER THAN OR EQUAL TO 15. */
/* ----------------------------------------------------------------------- */
/* ------------------------ */
/*     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP */
/*            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. */
/*            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1. */

/* ------------------------ */
/*     E0 = 2/SQRT(PI) */
/*     E1 = 2**(-3/2) */
/* ------------------------ */
/* ------------------------ */
    ret_val = 0.0;
    if (a >= b) {
	goto L10;
    }
    h__ = a / b;
    r0 = 1.0 / (h__ + 1.0);
    r1 = (b - a) / b;
    w0 = 1.0 / sqrt(a * (h__ + 1.0));
    goto L20;
L10:
    h__ = b / a;
    r0 = 1.0 / (h__ + 1.0);
    r1 = (b - a) / a;
    w0 = 1.0 / sqrt(b * (h__ + 1.0));

L20:
    r__1 = -lambda / a;
    r__2 = lambda / b;
    f = a * rlog1(r__1) + b * rlog1(r__2);
    t = exp(-f);
    if (t == 0.0) {
	return ret_val;
    }
    z0 = sqrt(f);
    z__ = z0 / e1 * 0.5;
    z2 = f + f;

    a0[0] = r1 * .66666666666666663;
    c__[0] = a0[0] * -0.5;
    d__[0] = -c__[0];
    j0 = 0.5 / e0 * erfc1(1, z0);
    j1 = e1;
    sum = j0 + d__[0] * w0 * j1;

    s = 1.0;
    h2 = h__ * h__;
    hn = 1.0;
    w = w0;
    znm1 = z__;
    zn = z2;
    i__1 = num;
    for (n = 2; n <= i__1; n += 2) {
	hn = h2 * hn;
	a0[n - 1] = r0 * 2.0 * (h__ * hn + 1.0) / (n + 2.0);
	np1 = n + 1;
	s += hn;
	a0[np1 - 1] = r1 * 2.0 * s / (n + 3.0);

	i__2 = np1;
	for (i__ = n; i__ <= i__2; ++i__) {
	    r__ = (i__ + 1.0) * -0.5;
	    b0[0] = r__ * a0[0];
	    i__3 = i__;
	    for (m = 2; m <= i__3; ++m) {
		bsum = 0.0;
		mm1 = m - 1;
		i__4 = mm1;
		for (j = 1; j <= i__4; ++j) {
		    mmj = m - j;
/* L30: */
		    bsum += (j * r__ - mmj) * a0[j - 1] * b0[mmj - 1];
		}
/* L31: */
		b0[m - 1] = r__ * a0[m - 1] + bsum / m;
	    }
	    c__[i__ - 1] = b0[i__ - 1] / (i__ + 1.0);

	    dsum = 0.0;
	    im1 = i__ - 1;
	    i__3 = im1;
	    for (j = 1; j <= i__3; ++j) {
		imj = i__ - j;
/* L40: */
		dsum += d__[imj - 1] * c__[j - 1];
	    }
/* L41: */
	    d__[i__ - 1] = -(dsum + c__[i__ - 1]);
	}

	j0 = e1 * znm1 + (n - 1.0) * j0;
	j1 = e1 * zn + n * j1;
	znm1 = z2 * znm1;
	zn = z2 * zn;
	w = w0 * w;
	t0 = d__[n - 1] * w * j0;
	w = w0 * w;
	t1 = d__[np1 - 1] * w * j1;
	sum += t0 + t1;
	if (fabs(t0) + fabs(t1) <= eps * sum) {
	    goto L60;
	}
/* L50: */
    }

L60:
    u = exp(-bcorr(a, b));
    ret_val = e0 * t * u * sum;
    return ret_val;
} /* basym_ */


static double exparg(int l)
{
    /* Local variables */
    int m;
    double lnb;

/* -------------------------------------------------------------------- */
/*     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH */
/*     EXP(W) CAN BE COMPUTED. */

/*     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR */
/*     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO. */

/*     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED. */
/* -------------------------------------------------------------------- */

    lnb = .69314718055995;
    if (l == 0) {
	m = Rf_i1mach(16);
	return m * lnb * .99999;
    }
    m = Rf_i1mach(15) - 1;
    return m * lnb * .99999;
} /* exparg */

static double esum(int mu, double x)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    double w;

/* ----------------------------------------------------------------------- */
/*                    EVALUATION OF EXP(MU + X) */
/* ----------------------------------------------------------------------- */
    if (x > 0.0) {
	goto L10;
    }

    if (mu < 0) {
	goto L20;
    }
    w = mu + x;
    if (w > 0.0) {
	goto L20;
    }
    ret_val = exp(w);
    return ret_val;

L10:
    if (mu > 0) {
	goto L20;
    }
    w = mu + x;
    if (w < 0.0) {
	goto L20;
    }
    ret_val = exp(w);
    return ret_val;

L20:
    w = (double) (mu);
    ret_val = exp(w) * exp(x);
    return ret_val;
} /* esum */

double rexp(double x)
{
    /* Initialized data */

    static double p1 = 9.14041914819518e-10;
    static double p2 = .0238082361044469;
    static double q1 = -.499999999085958;
    static double q2 = .107141568980644;
    static double q3 = -.0119041179760821;
    static double q4 = 5.95130811860248e-4;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double w;

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION EXP(X) - 1 */
/* ----------------------------------------------------------------------- */
/* ----------------------- */
    if (fabs(x) > .15f) {
	goto L10;
    }
    ret_val = x * (((p2 * x + p1) * x + 1.0) / ((((q4 * x + q3) * x + q2)
	     * x + q1) * x + 1.0));
    return ret_val;

L10:
    w = exp(x);
    if (x > 0.0) {
	goto L20;
    }
    ret_val = w - 0.5 - 0.5;
    return ret_val;
L20:
    ret_val = w * (0.5 - 1.0 / w + 0.5);
    return ret_val;
} /* rexp */

static double alnrel(double a)
{
    /* Initialized data */

    static double p1 = -1.29418923021993;
    static double p2 = .405303492862024;
    static double p3 = -.0178874546012214;
    static double q1 = -1.62752256355323;
    static double q2 = .747811014037616;
    static double q3 = -.0845104217945565;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double t, w, x, t2;

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION LN(1 + A) */
/* ----------------------------------------------------------------------- */
/* -------------------------- */
    if (fabs(a) > .375f) {
	goto L10;
    }
    t = a / (a + 2.0);
    t2 = t * t;
    w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.0) / (((q3 * t2 + q2) * t2 + q1)
	    * t2 + 1.0);
    ret_val = t * 2.0 * w;
    return ret_val;

L10:
    x = a + 1.;
    ret_val = log(x);
    return ret_val;
} /* alnrel */

static double rlog1(double x)
{
    /* Initialized data */

    static double a = .0566749439387324;
    static double b = .0456512608815524;
    static double p0 = .333333333333333;
    static double p1 = -.224696413112536;
    static double p2 = .00620886815375787;
    static double q1 = -1.27408923933623;
    static double q2 = .354508718369557;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double h__, r__, t, w, w1;

/* ----------------------------------------------------------------------- */
/*             EVALUATION OF THE FUNCTION X - LN(1 + X) */
/* ----------------------------------------------------------------------- */
/* ------------------------ */
/* ------------------------ */
    if (x < -.39f || x > .57f) {
	goto L100;
    }
    if (x < -.18f) {
	goto L10;
    }
    if (x > .18f) {
	goto L20;
    }

/*              ARGUMENT REDUCTION */

    h__ = x;
    w1 = 0.0;
    goto L30;

L10:
    h__ = x + .3;
    h__ /= .7;
    w1 = a - h__ * .3;
    goto L30;

L20:
    h__ = x * .75 - .25;
    w1 = b + h__ / 3.0;

/*               SERIES EXPANSION */

L30:
    r__ = h__ / (h__ + 2.0);
    t = r__ * r__;
    w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0);
    ret_val = t * 2.0 * (1.0 / (1.0 - r__) - r__ * w) + w1;
    return ret_val;


L100:
    w = x + 0.5 + 0.5;
    ret_val = x - log(w);
    return ret_val;
} /* rlog1 */

static double erf__(double x)
{
    /* Initialized data */

    static double c__ = .564189583547756;
    static double a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static double b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static double p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static double q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static double r__[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static double s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double t, x2, ax, bot, top;

/* ----------------------------------------------------------------------- */
/*             EVALUATION OF THE REAL ERROR FUNCTION */
/* ----------------------------------------------------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
    ax = fabs(x);
    if (ax > 0.5) {
	goto L10;
    }
    t = x * x;
    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
    ret_val = x * (top / bot);
    return ret_val;

L10:
    if (ax > 4.0) {
	goto L20;
    }
    top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
	    + p[5]) * ax + p[6]) * ax + p[7];
    bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
	    + q[5]) * ax + q[6]) * ax + q[7];
    ret_val = 0.5 - exp(-x * x) * top / bot + 0.5;
    if (x < 0.0) {
	ret_val = -ret_val;
    }
    return ret_val;

L20:
    if (ax >= 5.8) {
	goto L30;
    }
    x2 = x * x;
    t = 1.0 / x2;
    top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
    ret_val = (c__ - top / (x2 * bot)) / ax;
    ret_val = 0.5 - exp(-x2) * ret_val + 0.5;
    if (x < 0.0) {
	ret_val = -ret_val;
    }
    return ret_val;

L30:
    ret_val = x > 0 ? 1 : -1;
    return ret_val;
} /* erf__ */

static double erfc1(int ind, double x)
{
    /* Initialized data */

    static double c__ = .564189583547756;
    static double a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static double b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static double p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static double q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static double r__[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static double s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double e, t, w, ax, bot, top;

/* ----------------------------------------------------------------------- */
/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
/* ----------------------------------------------------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */

/*                     ABS(X) .LE. 0.5 */

    ax = fabs(x);
    if (ax > 0.5) {
	goto L10;
    }
    t = x * x;
    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
    ret_val = 0.5 - x * (top / bot) + 0.5;
    if (ind != 0) {
	ret_val = exp(t) * ret_val;
    }
    return ret_val;

/*                  0.5 .LT. ABS(X) .LE. 4 */

L10:
    if (ax > 4.0) {
	goto L20;
    }
    top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
	    + p[5]) * ax + p[6]) * ax + p[7];
    bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
	    + q[5]) * ax + q[6]) * ax + q[7];
    ret_val = top / bot;
    goto L40;

/*                      ABS(X) .GT. 4 */

L20:
    if (x <= -5.6f) {
	goto L50;
    }
    if (ind != 0) {
	goto L30;
    }
    if (x > 100.0) {
	goto L60;
    }
    if (x * x > -exparg(1)) {
	goto L60;
    }

L30:
/* Computing 2nd power */
    d__1 = 1.0 / x;
    t = d__1 * d__1;
    top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
    ret_val = (c__ - t * top / bot) / ax;

/*                      FINAL ASSEMBLY */

L40:
    if (ind == 0) {
	goto L41;
    }
    if (x < 0.0) {
	ret_val = exp(x * x) * 2.0 - ret_val;
    }
    return ret_val;
L41:
    w = (double) (x) * (double) (x);
    t = w;
    e = w - t;
    ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
    if (x < 0.0) {
	ret_val = 2.0 - ret_val;
    }
    return ret_val;

/*             LIMIT VALUE FOR LARGE NEGATIVE X */

L50:
    ret_val = 2.0;
    if (ind != 0) {
	ret_val = exp(x * x) * 2.0;
    }
    return ret_val;

/*             LIMIT VALUE FOR LARGE POSITIVE X */
/*                       WHEN IND = 0 */

L60:
    ret_val = 0.0;
    return ret_val;
} /* erfc1 */

static double gam1(double a)
{
    /* Initialized data */

    static double p[7] = { .577215664901533,-.409078193005776,
	    -.230975380857675,.0597275330452234,.0076696818164949,
	    -.00514889771323592,5.89597428611429e-4 };
    static double q[5] = { 1.,.427569613095214,.158451672430138,
	    .0261132021441447,.00423244297896961 };
    static double r__[9] = { -.422784335098468,-.771330383816272,
	    -.244757765222226,.118378989872749,9.30357293360349e-4,
	    -.0118290993445146,.00223047661158249,2.66505979058923e-4,
	    -1.32674909766242e-4 };
    static double s1 = .273076135303957;
    static double s2 = .0559398236957378;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double d__, t, w, bot, top;

/*     ------------------------------------------------------------------ */
/*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5 */
/*     ------------------------------------------------------------------ */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
    t = a;
    d__ = a - 0.5;
    if (d__ > 0.0) {
	t = d__ - 0.5;
    }
    if (t < 0.0) {
	goto L30;
    } else if (t == 0) {
	goto L10;
    } else {
	goto L20;
    }

L10:
    ret_val = 0.0;
    return ret_val;

L20:
    top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]) * t + p[1]
	    ) * t + p[0];
    bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.0;
    w = top / bot;
    if (d__ > 0.0) {
	goto L21;
    }
    ret_val = a * w;
    return ret_val;
L21:
    ret_val = t / a * (w - 0.5 - 0.5);
    return ret_val;

L30:
    top = (((((((r__[8] * t + r__[7]) * t + r__[6]) * t + r__[5]) * t + r__[4]
	    ) * t + r__[3]) * t + r__[2]) * t + r__[1]) * t + r__[0];
    bot = (s2 * t + s1) * t + 1.0;
    w = top / bot;
    if (d__ > 0.0) {
	goto L31;
    }
    ret_val = a * (w + 0.5 + 0.5);
    return ret_val;
L31:
    ret_val = t * w / a;
    return ret_val;
} /* gam1 */

static double gamln1(double a)
{
    /* Initialized data */

    static double p0 = .577215664901533;
    static double p1 = .844203922187225;
    static double p2 = -.168860593646662;
    static double p3 = -.780427615533591;
    static double p4 = -.402055799310489;
    static double p5 = -.0673562214325671;
    static double p6 = -.00271935708322958;
    static double q1 = 2.88743195473681;
    static double q2 = 3.12755088914843;
    static double q3 = 1.56875193295039;
    static double q4 = .361951990101499;
    static double q5 = .0325038868253937;
    static double q6 = 6.67465618796164e-4;
    static double r0 = .422784335098467;
    static double r1 = .848044614534529;
    static double r2 = .565221050691933;
    static double r3 = .156513060486551;
    static double r4 = .017050248402265;
    static double r5 = 4.97958207639485e-4;
    static double s1 = 1.24313399877507;
    static double s2 = .548042109832463;
    static double s3 = .10155218743983;
    static double s4 = .00713309612391;
    static double s5 = 1.16165475989616e-4;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double w, x;

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25 */
/* ----------------------------------------------------------------------- */
/* ---------------------- */
/* ---------------------- */
    if (a >= .6f) {
	goto L10;
    }
    w = ((((((p6 * a + p5) * a + p4) * a + p3) * a + p2) * a + p1) * a
	    + p0) / ((((((q6 * a + q5) * a + q4) * a + q3) * a + q2) * a
	    + q1) * a + 1.0);
    ret_val = -(a) * w;
    return ret_val;

L10:
    x = a - 0.5 - 0.5;
    w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) / (((((s5 *
	    x + s4) * x + s3) * x + s2) * x + s1) * x + 1.0);
    ret_val = x * w;
    return ret_val;
} /* gamln1 */

static double psi(double x)
{
    /* Initialized data */

    static double piov4 = .785398163397448;
    static double dx0 = 1.461632144968362341262659542325721325;
    static double p1[7] = { .0089538502298197,4.77762828042627,
	    142.441585084029,1186.45200713425,3633.51846806499,
	    4138.10161269013,1305.60269827897 };
    static double q1[6] = { 44.8452573429826,520.752771467162,
	    2210.0079924783,3641.27349079381,1908.310765963,
	    6.91091682714533e-6 };
    static double p2[4] = { -2.12940445131011,-7.01677227766759,
	    -4.48616543918019,-.648157123766197 };
    static double q2[4] = { 32.2703493791143,89.2920700481861,
	    54.6117738103215,7.77788548522962 };

    /* System generated locals */
    double ret_val, d__1, d__2;

    /* Local variables */
    int i__, m, n;
    double w, z__;
    int nq;
    double den, aug, sgn, xmx0, xmax1, upper;
    double xsmall;

/* --------------------------------------------------------------------- */

/*                 EVALUATION OF THE DIGAMMA FUNCTION */

/*                           ----------- */

/*     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT */
/*     BE COMPUTED. */

/*     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV */
/*     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY */
/*     CODY, STRECOK AND THACHER. */

/* --------------------------------------------------------------------- */
/*     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK */
/*     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY */
/*     A.H. MORRIS (NSWC). */
/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     PIOV4 = PI/4 */
/*     DX0 = ZERO OF PSI TO EXTENDED PRECISION */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0 */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0 */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     MACHINE DEPENDENT CONSTANTS ... */

/*        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT */
/*                 WITH ENTIRELY INT REPRESENTATION.  ALSO USED */
/*                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE */
/*                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH */
/*                 PSI MAY BE REPRESENTED AS LOG(X). */

/*        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X) */
/*                 MAY BE REPRESENTED BY 1/X. */

/* --------------------------------------------------------------------- */
    xmax1 = Rf_d1mach(4) - 1.0;
/* Computing MIN */
    d__1 = xmax1, d__2 = 0.5 / Rf_d1mach(3);;
    xmax1 = min(d__1,d__2);
    xsmall = 1e-9;
/* --------------------------------------------------------------------- */
    aug = 0.0;
    if (x >= 0.5) {
	goto L200;
    }
/* --------------------------------------------------------------------- */
/*     X .LT. 0.5,  USE REFLECTION FORMULA */
/*     PSI(1-X) = PSI(X) + PI * COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    if (fabs(x) > xsmall) {
	goto L100;
    }
    if (x == 0.f) {
	goto L400;
    }
/* --------------------------------------------------------------------- */
/*     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE */
/*     FOR  PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    aug = -1.0 / x;
    goto L150;
/* --------------------------------------------------------------------- */
/*     REDUCTION OF ARGUMENT FOR COTAN */
/* --------------------------------------------------------------------- */
L100:
    w = -x;
    sgn = piov4;
    if (w > 0.0) {
	goto L120;
    }
    w = -w;
    sgn = -sgn;
/* --------------------------------------------------------------------- */
/*     MAKE AN ERROR EXIT IF X .LE. -XMAX1 */
/* --------------------------------------------------------------------- */
L120:
    if (w >= xmax1) {
	goto L400;
    }
    nq = (int) w;
    w -= (double) nq;
    nq = (int) (w * 4.0);
    w = (w - (double) nq * .25f) * 4.0;
/* --------------------------------------------------------------------- */
/*     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X. */
/*     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST */
/*     QUADRANT AND DETERMINE SIGN */
/* --------------------------------------------------------------------- */
    n = nq / 2;
    if (n + n != nq) {
	w = 1.0 - w;
    }
    z__ = piov4 * w;
    m = n / 2;
    if (m + m != n) {
	sgn = -sgn;
    }
/* --------------------------------------------------------------------- */
/*     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    n = (nq + 1) / 2;
    m = n / 2;
    m += m;
    if (m != n) {
	goto L140;
    }
/* --------------------------------------------------------------------- */
/*     CHECK FOR SINGULARITY */
/* --------------------------------------------------------------------- */
    if (z__ == 0.0) {
	goto L400;
    }
/* --------------------------------------------------------------------- */
/*     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND */
/*     SIN/COS AS A SUBSTITUTE FOR TAN */
/* --------------------------------------------------------------------- */
    aug = sgn * (cos(z__) / sin(z__) * 4.0);
    goto L150;
L140:
    aug = sgn * (sin(z__) / cos(z__) * 4.0);
L150:
    x = 1.0 - x;
L200:
    if (x > 3.0) {
	goto L300;
    }
/* --------------------------------------------------------------------- */
/*     0.5 .LE. X .LE. 3.0 */
/* --------------------------------------------------------------------- */
    den = x;
    upper = p1[0] * x;

    for (i__ = 1; i__ <= 5; ++i__) {
	den = (den + q1[i__ - 1]) * x;
	upper = (upper + p1[i__]) * x;
/* L210: */
    }

    den = (upper + p1[6]) / (den + q1[5]);
    xmx0 = x - dx0;
    ret_val = den * xmx0 + aug;
    return ret_val;
/* --------------------------------------------------------------------- */
/*     IF X .GE. XMAX1, PSI = LN(X) */
/* --------------------------------------------------------------------- */
L300:
    if (x >= xmax1) {
	goto L350;
    }
/* --------------------------------------------------------------------- */
/*     3.0 .LT. X .LT. XMAX1 */
/* --------------------------------------------------------------------- */
    w = 1.0 / (x * x);
    den = w;
    upper = p2[0] * w;

    for (i__ = 1; i__ <= 3; ++i__) {
	den = (den + q2[i__ - 1]) * w;
	upper = (upper + p2[i__]) * w;
/* L310: */
    }

    aug = upper / (den + q2[3]) - 0.5 / x + aug;
L350:
    ret_val = aug + log(x);
    return ret_val;
/* --------------------------------------------------------------------- */
/*     ERROR RETURN */
/* --------------------------------------------------------------------- */
L400:
    ret_val = 0.0;
    return ret_val;
} /* psi */

static double betaln(double a0, double b0)
{
    /* Initialized data */

    static double e = .918938533204673;

    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    double a, b, c__, h__;
    int i__, n;
    double u, v, w, z__;

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION */
/* ----------------------------------------------------------------------- */
/*     E = 0.5*LN(2*PI) */
/* -------------------------- */
/* -------------------------- */
    a = min(a0 ,b0);
    b = max(a0, b0);
    if (a >= 8.0) {
	goto L60;
    }
    if (a >= 1.0) {
	goto L20;
    }
/* ----------------------------------------------------------------------- */
/*                   PROCEDURE WHEN A .LT. 1 */
/* ----------------------------------------------------------------------- */
    if (b >= 8.0) {
	goto L10;
    }
    ret_val = gamln(a) + (gamln(b) - gamln(a+b));
    return ret_val;
L10:
    ret_val = gamln(a) + algdiv(a, b);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*                PROCEDURE WHEN 1 .LE. A .LT. 8 */
/* ----------------------------------------------------------------------- */
L20:
    if (a > 2.0) {
	goto L30;
    }
    if (b > 2.0) {
	goto L21;
    }
    ret_val = gamln(a) + gamln(b) - gsumln(a, b);
    return ret_val;
L21:
    w = 0.0;
    if (b < 8.0) {
	goto L40;
    }
    ret_val = gamln(a) + algdiv(a, b);
    return ret_val;

/*                REDUCTION OF A WHEN B .LE. 1000 */

L30:
    if (b > 1e3f) {
	goto L50;
    }
    n = a - 1.0;
    w = 1.0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a += -1.0;
	h__ = a / b;
	w *= h__ / (h__ + 1.0);
/* L31: */
    }
    w = log(w);
    if (b < 8.0) {
	goto L40;
    }
    ret_val = w + gamln(a) + algdiv(a, b);
    return ret_val;

/*                 REDUCTION OF B WHEN B .LT. 8 */

L40:
    n = b - 1.0;
    z__ = 1.0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b += -1.0;
	z__ *= b / (a + b);
/* L41: */
    }
    ret_val = w + log(z__) + (gamln(a) + (gamln(b) - gsumln(a, b)));
    return ret_val;

/*                REDUCTION OF A WHEN B .GT. 1000 */

L50:
    n = a - 1.0;
    w = 1.0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a += -1.0;
	w *= a / (a / b + 1.0);
/* L51: */
    }
    ret_val = log(w) - n * log(b) + (gamln(a) + algdiv(a, b));
    return ret_val;
/* ----------------------------------------------------------------------- */
/*                   PROCEDURE WHEN A .GE. 8 */
/* ----------------------------------------------------------------------- */
L60:
    w = bcorr(a, b);
    h__ = a / b;
    c__ = h__ / (h__ + 1.0);
    u = -(a - 0.5) * log(c__);
    v = b * alnrel(h__);
    if (u <= v) {
	goto L61;
    }
    ret_val = log(b) * -0.5 + e + w - v - u;
    return ret_val;
L61:
    ret_val = log(b) * -0.5 + e + w - u - v;
    return ret_val;
} /* betaln */

static double gsumln(double a, double b)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    double x;

/* ----------------------------------------------------------------------- */
/*          EVALUATION OF THE FUNCTION LN(GAMMA(A + B)) */
/*          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2 */
/* ----------------------------------------------------------------------- */
    x = a + b - 2.;
    if (x > 0.25) {
	goto L10;
    }
    ret_val = gamln1(x + 1.0);
    return ret_val;
L10:
    if (x > 1.25) {
	goto L20;
    }
    ret_val = gamln1(x) + alnrel(x);
    return ret_val;
L20:
    ret_val = gamln1(x - 1.0) + log(x * (x + 1.0));
    return ret_val;
} /* gsumln */

static double bcorr(double a0, double b0)
{
    /* Initialized data */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    /* System generated locals */
    double ret_val, r__1;

    /* Local variables */
    double a, b, c__, h__, t, w, x, s3, s5, x2, s7, s9, s11;

/* ----------------------------------------------------------------------- */

/*     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE */
/*     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A). */
/*     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8. */

/* ----------------------------------------------------------------------- */
/* ------------------------ */
    a = min(a0, b0);
    b = max(a0, b0);

    h__ = a / b;
    c__ = h__ / (h__ + 1.0);
    x = 1.0 / (h__ + 1.0);
    x2 = x * x;

/*                SET SN = (1 - X**N)/(1 - X) */

    s3 = x + x2 + 1.0;
    s5 = x + x2 * s3 + 1.0;
    s7 = x + x2 * s5 + 1.0;
    s9 = x + x2 * s7 + 1.0;
    s11 = x + x2 * s9 + 1.0;

/*                SET W = DEL(B) - DEL(A + B) */

/* Computing 2nd power */
    r__1 = 1.0 / b;
    t = r__1 * r__1;
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
	    s3) * t + c0;
    w *= c__ / b;

/*                   COMPUTE  DEL(A) + W */

/* Computing 2nd power */
    r__1 = 1.0 / a;
    t = r__1 * r__1;
    ret_val = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a +
	    w;
    return ret_val;
} /* bcorr */

static double algdiv(double a, double b)
{
    /* Initialized data */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    /* System generated locals */
    double ret_val, r__1;

    /* Local variables */
    double c__, d__, h__, t, u, v, w, x, s3, s5, x2, s7, s9, s11;

/* ----------------------------------------------------------------------- */

/*     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B .GE. 8 */

/*                         -------- */

/*     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY */
/*     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X). */

/* ----------------------------------------------------------------------- */
/* ------------------------ */
    if (a <= b) {
	goto L10;
    }
    h__ = b / a;
    c__ = 1.0 / (h__ + 1.0);
    x = h__ / (h__ + 1.0);
    d__ = a + (b - 0.5);
    goto L20;
L10:
    h__ = a / b;
    c__ = h__ / (h__ + 1.0);
    x = 1.0 / (h__ + 1.0);
    d__ = b + (a - 0.5);

/*                SET SN = (1 - X**N)/(1 - X) */

L20:
    x2 = x * x;
    s3 = x + x2 + 1.0;
    s5 = x + x2 * s3 + 1.0;
    s7 = x + x2 * s5 + 1.0;
    s9 = x + x2 * s7 + 1.0;
    s11 = x + x2 * s9 + 1.0;

/*                SET W = DEL(B) - DEL(A + B) */

/* Computing 2nd power */
    r__1 = 1.0 / b;
    t = r__1 * r__1;
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
	    s3) * t + c0;
    w *= c__ / b;

/*                    COMBINE THE RESULTS */

    r__1 = a / b;
    u = d__ * alnrel(r__1);
    v = a * (log(b) - 1.0);
    if (u <= v) {
	goto L30;
    }
    ret_val = w - v - u;
    return ret_val;
L30:
    ret_val = w - u - v;
    return ret_val;
} /* algdiv */

static double gamln(double a)
{
    /* Initialized data */

    static double d__ = .418938533204673;
    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    /* System generated locals */
    int i__1;
    double ret_val, r__1;

    /* Local variables */
    int i__, n;
    double t, w;

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A */
/* ----------------------------------------------------------------------- */
/*     WRITTEN BY ALFRED H. MORRIS */
/*          NAVAL SURFACE WARFARE CENTER */
/*          DAHLGREN, VIRGINIA */
/* -------------------------- */
/*     D = 0.5*(LN(2*PI) - 1) */
/* -------------------------- */
/* -------------------------- */
/* ----------------------------------------------------------------------- */
    if (a > .8f) {
	goto L10;
    }
    ret_val = gamln1(a) - log(a);
    return ret_val;
L10:
    if (a > 2.25f) {
	goto L20;
    }
    t = a - 0.5 - 0.5;
    ret_val = gamln1(t);
    return ret_val;

L20:
    if (a >= 10.0) {
	goto L30;
    }
    n = a - 1.25;
    t = a;
    w = 1.0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t += -1.0;
/* L21: */
	w = t * w;
    }
    r__1 = t - 1.0;
    ret_val = gamln1(r__1) + log(w);
    return ret_val;

L30:
/* Computing 2nd power */
    r__1 = 1.0 / a;
    t = r__1 * r__1;
    w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;
    ret_val = d__ + w + (a - 0.5) * (log(a) - 1.0);
    return ret_val;
} /* gamln */
/*
 *  Mathlib - A Mathematical Function Library
 *  Copyright (C) 1998  Ross Ihaka
 *  Copyright (C) 2000 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* NaNs propagated correctly */


/*-- FIXME:  Eliminate calls to these
 *   =====   o   from C code when
 *	     o   it is only used to initialize "static" variables (threading)
 *  and use the DBL_... constants instead
 */


double Rf_d1mach(int i)
{
    switch(i) {
    case 1: return DBL_MIN;
    case 2: return DBL_MAX;

    case 3: /* = FLT_RADIX  ^ - DBL_MANT_DIG
	      for IEEE:  = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON */
	return pow((double)i1mach(10), -(double)i1mach(14));

    case 4: /* = FLT_RADIX  ^ (1- DBL_MANT_DIG) =
	      for IEEE:  = 2^52 = 4503599627370496 = 1/DBL_EPSILON */
	return pow((double)i1mach(10), 1-(double)i1mach(14));

    case 5: return log10(2.0);/* = M_LOG10_2 in Rmath.h */


    default: return 0.0;
    }
}

/*
 *  Mathlib - A Mathematical Function Library
 *  Copyright (C) 1998  Ross Ihaka
 *  Copyright (C) 2000 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 */


int Rf_i1mach(int i)
{
    switch(i) {

    case  1: return 5;
    case  2: return 6;
    case  3: return 0;
    case  4: return 0;

    case  5: return CHAR_BIT * sizeof(int);
    case  6: return sizeof(int)/sizeof(char);

    case  7: return 2;
    case  8: return CHAR_BIT * sizeof(int) - 1;
    case  9: return INT_MAX;

    case 10: return FLT_RADIX;

    case 11: return FLT_MANT_DIG;
    case 12: return FLT_MIN_EXP;
    case 13: return FLT_MAX_EXP;

    case 14: return DBL_MANT_DIG;
    case 15: return DBL_MIN_EXP;
    case 16: return DBL_MAX_EXP;

    default: return 0;
    }
}



/*
  Mathlib : A C Library of Special Functions
  Copyright (C) 1999-2007  The R Development Core Team

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or (at
  your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, a copy is available at
  http://www.r-project.org/Licenses/

  SYNOPSIS

    #include <Rmath.h>
    double dwilcox(double x, double m, double n, int give_log)
    double pwilcox(double x, double m, double n, int lower_tail, int log_p)
    double qwilcox(double x, double m, double n, int lower_tail, int log_p);
    double rwilcox(double m, double n)

  DESCRIPTION

    dwilcox	The density of the Wilcoxon distribution.
    pwilcox	The distribution function of the Wilcoxon distribution.
    qwilcox	The quantile function of the Wilcoxon distribution.
    rwilcox	Random variates from the Wilcoxon distribution.

 */

/* 
   Note: the checks here for R_CheckInterrupt also do stack checking.

   calloc/free are remapped for use in R, so allocation checks are done there.
   freeing is completed by an on.exit action in the R wrappers.
*/

//#include "nmath.h"
//#include "dpq.h"

//#ifndef MATHLIB_STANDALONE
//#include <R_ext/Utils.h>
//#endif

static double ***w; /* to store  cwilcox(i,j,k) -> w[i][j][k] */
static int allocated_m, allocated_n;

static void
w_free(int m, int n)
{
    int i, j;

    for (i = m; i >= 0; i--) {
	for (j = n; j >= 0; j--) {
	    if (w[i][j] != 0)
		free((void *) w[i][j]);
	}
	free((void *) w[i]);
    }
    free((void *) w);
    w = 0; allocated_m = allocated_n = 0;
}

static void
w_init_maybe(int m, int n)
{
    int i;

    if (m > n) {
	i = n; n = m; m = i;
    }
    if (w && (m > allocated_m || n > allocated_n))
	w_free(allocated_m, allocated_n); /* zeroes w */

    if (!w) { /* initialize w[][] */
	m = imax2(m, WILCOX_MAX);
	n = imax2(n, WILCOX_MAX);
	w = (double ***) calloc(m + 1, sizeof(double **));
#ifdef MATHLIB_STANDALONE
	if (!w) MATHLIB_ERROR(_("wilcox allocation error %d"), 1);
#endif
	for (i = 0; i <= m; i++) {
	    w[i] = (double **) calloc(n + 1, sizeof(double *));
#ifdef MATHLIB_STANDALONE
	    /* the apparent leak here in the in-R case should be
	       swept up by the on.exit action */
	    if (!w[i]) {
		/* first free all earlier allocations */
		w_free(i-1, n);
		MATHLIB_ERROR(_("wilcox allocation error %d"), 2);
	    }
#endif
	}
	allocated_m = m; allocated_n = n;
    }
}

static void
w_free_maybe(int m, int n)
{
    if (m > WILCOX_MAX || n > WILCOX_MAX)
	w_free(m, n);
}


/* This counts the number of choices with statistic = k */
static double
cwilcox(int k, int m, int n)
{
    int c, u, i, j, l;

    
    u = m * n;
    if (k < 0 || k > u)
	return(0);
    c = (int)(u / 2);
    if (k > c)
	k = u - k; /* hence  k <= floor(u / 2) */
    if (m < n) {
	i = m; j = n;
    } else {
	i = n; j = m;
    } /* hence  i <= j */

    if (j == 0) /* and hence i == 0 */
	return (k == 0);


    /* We can simplify things if k is small.  Consider the Mann-Whitney 
       definition, and sort y.  Then if the statistic is k, no more 
       than k of the y's can be <= any x[i], and since they are sorted 
       these can only be in the first k.  So the count is the same as
       if there were just k y's. 
    */
    if (j > 0 && k < j) return cwilcox(k, i, k);    
    
    if (w[i][j] == 0) {
	w[i][j] = (double *) calloc(c + 1, sizeof(double));
#ifdef MATHLIB_STANDALONE
	if (!w[i][j]) MATHLIB_ERROR(_("wilcox allocation error %d"), 3);
#endif
	for (l = 0; l <= c; l++)
	    w[i][j][l] = -1;
    }
    if (w[i][j][k] < 0) {
	if (j == 0) /* and hence i == 0 */
	    w[i][j][k] = (k == 0);
	else
	    w[i][j][k] = cwilcox(k - j, i - 1, j) + cwilcox(k, i, j - 1);

    }
    return(w[i][j][k]);
}

double dwilcox(double x, double m, double n, int give_log)
{
    double d;

#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(m) || ISNAN(n))
	return(x + m + n);
#endif
    m = floor(m + 0.5);
    n = floor(n + 0.5);
    if (m <= 0 || n <= 0)
	ML_ERR_return_NAN;

    if (fabs(x - floor(x + 0.5)) > 1e-7)
	return(R_D__0);
    x = floor(x + 0.5);
    if ((x < 0) || (x > m * n))
	return(R_D__0);

    w_init_maybe(m, n);
    d = give_log ?
	log(cwilcox(x, m, n)) - lchoose(m + n, n) :
	    cwilcox(x, m, n)  /	 choose(m + n, n);

    return(d);
}

/* args have the same meaning as R function pwilcox */
double pwilcox(double q, double m, double n, int lower_tail, int log_p)
{
    int i;
    double c, p;

#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(m) || ISNAN(n))
	return(q + m + n);
#endif
    if (!R_FINITE(m) || !R_FINITE(n))
	ML_ERR_return_NAN;
    m = floor(m + 0.5);
    n = floor(n + 0.5);
    if (m <= 0 || n <= 0)
	ML_ERR_return_NAN;

    q = floor(q + 1e-7);

    if (q < 0.0)
	return(R_DT_0);
    if (q >= m * n)
	return(R_DT_1);

    w_init_maybe(m, n);
    c = choose(m + n, n);
    p = 0;
    /* Use summation of probs over the shorter range */
    if (q <= (m * n / 2)) {
	for (i = 0; i <= q; i++)
	    p += cwilcox(i, m, n) / c;
    }
    else {
	q = m * n - q;
	for (i = 0; i < q; i++)
	    p += cwilcox(i, m, n) / c;
	lower_tail = !lower_tail; /* p = 1 - p; */
    }

    return(R_DT_val(p));
} /* pwilcox */

/* x is 'p' in R function qwilcox */

double qwilcox(double x, double m, double n, int lower_tail, int log_p)
{
    double c, p, q;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(m) || ISNAN(n))
	return(x + m + n);
#endif
    if(!R_FINITE(x) || !R_FINITE(m) || !R_FINITE(n))
	ML_ERR_return_NAN;
    //R_Q_P01_check(x);

    m = floor(m + 0.5);
    n = floor(n + 0.5);
    if (m <= 0 || n <= 0)
	ML_ERR_return_NAN;

    if (x == R_DT_0)
	return(0);
    if (x == R_DT_1)
	return(m * n);

    if(log_p || !lower_tail)
	x = R_DT_qIv(x); /* lower_tail,non-log "p" */

    w_init_maybe(m, n);
    c = choose(m + n, n);
    p = 0;
    q = 0;
    if (x <= 0.5) {
	x = x - 10 * DBL_EPSILON;
	for (;;) {
	    p += cwilcox(q, m, n) / c;
	    if (p >= x)
		break;
	    q++;
	}
    }
    else {
	x = 1 - x + 10 * DBL_EPSILON;
	for (;;) {
	    p += cwilcox(q, m, n) / c;
	    if (p > x) {
		q = m * n - q;
		break;
	    }
	    q++;
	}
    }

    return(q);
}



void wilcox_free(void)
{
    w_free_maybe(allocated_m, allocated_n);
}

int imax2(int x, int y)
{
    return (x < y) ? y : x;
}

 
double wilcoxtest(double* x, int nx, double* y, int ny, int alter)
{
  
  double* r;
  int ties  = 0;
  int exact = 0;
  double stat = 0;
  int i;
  int n = nx+ny;
  double p = -1;
  int correct = 1;

  if ((nx == 0) || (ny == 0)) 
    return 1.0;

  if ((nx < 50) && (ny < 50)) {
    exact = 1;
  }

  combinedRanks(x, nx, y, ny, &r, &ties);
  
  //for (i=0; i<n; i++) {
  //  printf(" %3.2f", r[i]);
  //}
  
  stat = 0.0;
  for (i=0; i<nx; i++) {
    stat += r[i];
  }
  stat = stat - nx * (nx + 1) / (double)2.0;
  /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/
  if ((exact == 1) && (ties == 0)) {
    
    if (alter == 2) {
      
      if (stat > (nx * ny / (double)2.0)) {
	p = pwilcox(stat - 1.0, nx, ny, 0, 0);
      } else {
	p = min(2*pwilcox(stat, nx, ny, 1, 0), 1);
      }      
    } // alter == 2 (two-tailed)
    else if (alter == 1) {      
      p = pwilcox(stat - 1.0, nx, ny, 0, 0);
    } else if (alter == 0) {
      p = pwilcox(stat, nx, ny, 1, 0);
    }
  } else {
    // 
    int* t; int nt;
    table(r, n, &t, &nt);
    //printf("nt=%d\n", nt);
    double wsum = 0.0;
    for (i=0; i<nt; i++) {
      wsum += ( t[i]* t[i]* t[i] - t[i] );
    }
    double z = stat - nx * ny / (double)2.0;
    double sigma = sqrt((nx*ny/(double)12.0) * 
			((nx + ny + 1)
			 - wsum 
			 / ((nx+ny) * ((double)nx+ny-1) )));
    
    double cor = 0.0;
    if (correct == 1) {
      if (alter == 2)
	cor = sign(z) * 0.5;
      else if (alter == 1)
	cor = 0.5;
      else 
	cor = -0.5;
    }

    z = ( z - cor )  / sigma;
    
    if (alter == 2) {
      double p1 = pnorm(z, 0, 1, 1, 0);
      double p2 = pnorm(z, 0, 1, 0, 0);
      p = 2.0 * min( p1, p2 );
    } else if (alter == 1)
      p = pnorm(z, 0, 1, 0, 0);
    else
      p = pnorm(z, 0, 1, 1, 0);
    
    free(t);
  } // correction
    
  

  free(r);

  return p;
}


double pnorm(double x, double mu, double sigma, int lower_tail, int log_p)
{
  return pnorm5(x,mu, sigma, lower_tail, log_p);
}

double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p)
{
    double p, cp;

    /* Note: The structure of these checks has been carefully thought through.
     * For example, if x == mu and sigma == 0, we get the correct answer 1.
     */
#ifdef IEEE_754
    if(ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
    if(!R_FINITE(x) && mu == x) return ML_NAN;/* x-mu is NaN */
    if (sigma <= 0) {
	if(sigma < 0) ML_ERR_return_NAN;
	/* sigma = 0 : */
	return (x < mu) ? R_DT_0 : R_DT_1;
    }
    p = (x - mu) / sigma;
    if(!R_FINITE(p))
	return (x < mu) ? R_DT_0 : R_DT_1;
    x = p;

    pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

    return(lower_tail ? p : cp);
}

#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)
{
/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/
    const static double a[5] = {
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
    };
    const static double b[4] = {
	47.20258190468824187,
	976.09855173777669322,
	10260.932208618978205,
	45507.789335026729956
    };
    const static double c[9] = {
	0.39894151208813466764,
	8.8831497943883759412,
	93.506656132177855979,
	597.27027639480026226,
	2494.5375852903726711,
	6848.1904505362823326,
	11602.651437647350124,
	9842.7148383839780218,
	1.0765576773720192317e-8
    };
    const static double d[8] = {
	22.266688044328115691,
	235.38790178262499861,
	1519.377599407554805,
	6485.558298266760755,
	18615.571640885098091,
	34900.952721145977266,
	38912.003286093271411,
	19685.429676859990727
    };
    const static double p[6] = {
	0.21589853405795699,
	0.1274011611602473639,
	0.022235277870649807,
	0.001421619193227893466,
	2.9112874951168792e-5,
	0.02307344176494017303
    };
    const static double q[5] = {
	1.28426009614491121,
	0.468238212480865118,
	0.0659881378689285515,
	0.00378239633202758244,
	7.29751555083966205e-5
    };

    double xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
    double min = DBL_MIN;
#endif
    int i, lower, upper;

#ifdef IEEE_754
    if(ISNAN(x)) { *cum = *ccum = x; return; }
#endif

    /* Consider changing these : */
    eps = DBL_EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = fabs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	if (y > eps) {
	    xsq = x * x;
	    xnum = a[4] * xsq;
	    xden = xsq;
	    for (i = 0; i < 3; ++i) {
		xnum = (xnum + a[i]) * xsq;
		xden = (xden + b[i]) * xsq;
	    }
	} else xnum = xden = 0.0;

	temp = x * (xnum + a[3]) / (xden + b[3]);
	if(lower)  *cum = 0.5 + temp;
	if(upper) *ccum = 0.5 - temp;
	if(log_p) {
	    if(lower)  *cum = log(*cum);
	    if(upper) *ccum = log(*ccum);
	}
    }
    else if (y <= M_SQRT_32) {

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

	xnum = c[8] * y;
	xden = y;
	for (i = 0; i < 7; ++i) {
	    xnum = (xnum + c[i]) * y;
	    xden = (xden + d[i]) * y;
	}
	temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = trunc(X * SIXTEN) / SIXTEN;				\
	del = (X - xsq) * (X + xsq);					\
	if(log_p) {							\
	    *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);	\
	    if((lower && x > 0.) || (upper && x <= 0.))			\
		  *ccum = log1p(-exp(-xsq * xsq * 0.5) *		\
				exp(-del * 0.5) * temp);		\
	}								\
	else {								\
	    *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	    *ccum = 1.0 - *cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
	}

	do_del(y);
	swap_tail;
    }

/* else	  |x| > sqrt(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !

 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
    else if((log_p && y < 1e170) /* avoid underflow below */
	/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
	 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

	 xsq = x*x;

	 if(xsq * DBL_EPSILON < 1.)
	    del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
	 else
	    del = 0.;
	 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
	 *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

 	 swap_tail;

	 [Yes, but xsq might be infinite.]

	*/
	    || (lower && -37.5193 < x  &&  x < 8.2924)
	    || (upper && -8.2924  < x  &&  x < 37.5193)
	) {

	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
	xnum = p[5] * xsq;
	xden = xsq;
	for (i = 0; i < 4; ++i) {
	    xnum = (xnum + p[i]) * xsq;
	    xden = (xden + q[i]) * xsq;
	}
	temp = xsq * (xnum + p[4]) / (xden + q[4]);
	temp = (M_1_SQRT_2PI - temp) / y;

	do_del(x);
	swap_tail;
    } else { /* large x such that probs are 0 or 1 */
	if(x > 0) {	*cum = R_D__1; *ccum = R_D__0;	}
	else {	        *cum = R_D__0; *ccum = R_D__1;	}
    }


#ifdef NO_DENORMS
    /* do not return "denormalized" -- we do in R */
    if(log_p) {
	if(*cum > -min)	 *cum = -0.;
	if(*ccum > -min)*ccum = -0.;
    }
    else {
	if(*cum < min)	 *cum = 0.;
	if(*ccum < min)	*ccum = 0.;
    }
#endif
    return;
}
