#include <math.h>
#include <stdio.h>

#define R_DT_0 (lower_tail ? R_D__0 : R_D__1)          /* 0 */
#define R_DT_1 (lower_tail ? R_D__1 : R_D__0)          /* 1 */
#define R_Log1_Exp(x) ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
#define DBL_EPSILON 2.2204460492503131e-16
#define DBL_MIN 2.2250738585072014e-308
#define DBL_MAX 1.7976931348623157e+308
#define DBL_MANT_DIG 53
#define DBL_MIN_EXP (-1021)
#define DBL_MAX_EXP 1024
#define ML_UNDERFLOW (DBL_MIN * DBL_MIN)
#define M_LN2 0.693147180559945309417232121458
#define ML_NEGINF       ((-1.0) / 0.0)
#define R_D__0    (log_p ? ML_NEGINF : 0.)                /* 0 */
#define R_D__1    (log_p ? 0. : 1.)                       /* 1 */
#ifndef NAN
static const double NAN = (0.0 / 0.0);
#endif


int signgam;

double ppois(double x, double lambda, int lower_tail, int log_p);
double pgamma(double x, double alph, double scale, int lower_tail, int log_p);
static double pgamma_raw (double x, double alph, int lower_tail, int log_p);
static double pd_upper_series (double x, double y, int log_p);
static double logcf (double x, double i, double d);
double log1pmx (double x);
double lgamma1p (double a);
double logspace_add (double logx, double logy);
double logspace_sub (double logx, double logy);
static double dpois_wrap (double x_plus_1, double lambda, int give_log);
static double pgamma_smallx (double x, double alph, int lower_tail, int log_p);
static double pd_lower_cf (double i, double d);
static double pd_lower_series (double lambda, double y);
static double ppois_asymp (double x, double lambda, int lower_tail, int log_p);

/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define SQR(x) ((x)*(x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

double ppois(double x, double lambda, int lower_tail, int log_p)
{

    x = floor(x + 1e-7);

    if (x < 0)          return R_DT_0;
    if (lambda == 0.)   return R_DT_1;

    return pgamma(lambda, x + 1, 1., !lower_tail, log_p);
}

double pgamma(double x, double alph, double scale, int lower_tail, int log_p)
{

    if(alph <= 0. || scale <= 0.)
	return NAN;

    x /= scale;

    if (x <= 0.) /* also for scale=Inf and finite x */
	return R_DT_0;

    return pgamma_raw (x, alph, lower_tail, log_p);
}



static double
pgamma_raw (double x, double alph, int lower_tail, int log_p)
{
    double res;

    if (x < 1) {
	res = pgamma_smallx (x, alph, lower_tail, log_p);
    } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {/* incl. large alph */
	double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
	double d = dpois_wrap (alph, x, log_p);

	if (!lower_tail)
	    res = log_p
		? R_Log1_Exp (d + sum)
		: 1 - d * sum;
	else
	    res = log_p ? sum + d : sum * d;
    } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {/* incl. large x */
	double sum;
	double d = dpois_wrap (alph, x, log_p);

	if (alph < 1) {
	    if (x * DBL_EPSILON > 1 - alph)
		sum = R_D__1;
	    else {
		double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
		/* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
		sum = log_p ? log (f) : f;
	    }
	} else {
	    sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
	    sum = log_p ? log1p (sum) : 1 + sum;
	}

	if (!lower_tail)
	    res = log_p ? sum + d : sum * d;
	else
	    res = log_p
		? R_Log1_Exp (d + sum)
		: 1 - d * sum;
    } else {

	res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
    }

    /*
     * We lose a fair amount of accuracy to underflow in the cases
     * where the final result is very close to DBL_MIN.	 In those
     * cases, simply redo via log space.
     */
    if (!log_p && res < DBL_MIN / DBL_EPSILON) {
	/* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */

	return exp (pgamma_raw (x, alph, lower_tail, 1));
    } else
	return res;
}


static double pd_upper_series (double x, double y, int log_p)
{
    double term = x / y;
    double sum = term;

    do {
	y++;
	term *= x / y;
	sum += term;
    } while (term > sum * DBL_EPSILON);

    /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
     *     =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
     *     =  x/y * (1 + \sum_{n=1}^oo  x^n / ((y+1)*...*(y+n)))
     *     ~  x/y +  o(x/y)   {which happens when alph -> Inf}
     */
    return log_p ? log (sum) : sum;
}



/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]  =~=  -x */
static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ...
 *
 * auxilary in log1pmx() and lgamma1p()
 */
static double logcf (double x, double i, double d)
{
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;
    const double cfVSmall = 1.0e-14;/* ~= relative tolerance */

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    while (fabs (a2 * b1 - a1 * b2) > fabs (cfVSmall * b1 * b2)) {
	double c3 = c2*c2*x;
	c2 += d;
	c4 += d;
	a1 = c4 * a2 - c3 * a1;
	b1 = c4 * b2 - c3 * b1;

	c3 = c1 * c1 * x;
	c1 += d;
	c4 += d;
	a2 = c4 * a1 - c3 * a2;
	b2 = c4 * b1 - c3 * b2;

	if (fabs (b2) > scalefactor) {
	    a1 /= scalefactor;
	    b1 /= scalefactor;
	    a2 /= scalefactor;
	    b2 /= scalefactor;
	} else if (fabs (b2) < 1 / scalefactor) {
	    a1 *= scalefactor;
	    b1 *= scalefactor;
	    a2 *= scalefactor;
	    b2 *= scalefactor;
	}
    }

    return a2 / b2;
}

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double log1pmx (double x)
{
    static const double minLog1Value = -0.79149064;
    static const double two = 2;

    if (x > 1 || x < minLog1Value)
	return log1p(x) - x;
    else { /* expand in	 [x/(2+x)]^2 */
	double term = x / (2 + x);
	double y = term * term;
	if (fabs(x) < 1e-2)
	    return term * ((((two / 9 * y + two / 7) * y + two / 5) * y +
			    two / 3) * y - x);
	else
	    return term * (2 * y * logcf (y, 3, 2) - x);
    }
}


/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double lgamma1p (double a)
{
    const double eulers_const =	 0.5772156649015328606065120900824024;

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 1:N, N = 40 : */
    const int N = 40;
    static const double coeffs[40] = {
	0.3224670334241132182362075833230126e-0,
	0.6735230105319809513324605383715000e-1,
	0.2058080842778454787900092413529198e-1,
	0.7385551028673985266273097291406834e-2,
	0.2890510330741523285752988298486755e-2,
	0.1192753911703260977113935692828109e-2,
	0.5096695247430424223356548135815582e-3,
	0.2231547584535793797614188036013401e-3,
	0.9945751278180853371459589003190170e-4,
	0.4492623673813314170020750240635786e-4,
	0.2050721277567069155316650397830591e-4,
	0.9439488275268395903987425104415055e-5,
	0.4374866789907487804181793223952411e-5,
	0.2039215753801366236781900709670839e-5,
	0.9551412130407419832857179772951265e-6,
	0.4492469198764566043294290331193655e-6,
	0.2120718480555466586923135901077628e-6,
	0.1004322482396809960872083050053344e-6,
	0.4769810169363980565760193417246730e-7,
	0.2271109460894316491031998116062124e-7,
	0.1083865921489695409107491757968159e-7,
	0.5183475041970046655121248647057669e-8,
	0.2483674543802478317185008663991718e-8,
	0.1192140140586091207442548202774640e-8,
	0.5731367241678862013330194857961011e-9,
	0.2759522885124233145178149692816341e-9,
	0.1330476437424448948149715720858008e-9,
	0.6422964563838100022082448087644648e-10,
	0.3104424774732227276239215783404066e-10,
	0.1502138408075414217093301048780668e-10,
	0.7275974480239079662504549924814047e-11,
	0.3527742476575915083615072228655483e-11,
	0.1711991790559617908601084114443031e-11,
	0.8315385841420284819798357793954418e-12,
	0.4042200525289440065536008957032895e-12,
	0.1966475631096616490411045679010286e-12,
	0.9573630387838555763782200936508615e-13,
	0.4664076026428374224576492565974577e-13,
	0.2273736960065972320633279596737272e-13,
	0.1109139947083452201658320007192334e-13
    };

    const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    double lgam;
    int i;

    if (fabs (a) >= 0.5)
	return lgammafn (a + 1);

    /* Abramowitz & Stegun 6.1.33,
     * also  http://functions.wolfram.com/06.11.06.0008.01 */
    lgam = c * logcf (-a / 2, N + 2, 1);
    for (i = N - 1; i >= 0; i--)
	lgam = coeffs[i] - a * lgam;

    return (a * lgam - eulers_const) * a - log1pmx (a);
} /* lgamma1p */



/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_add (double logx, double logy)
{
    return fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
}


/*
 * Compute the log of a difference from logs of terms, i.e.,
 *
 *     log (exp (logx) - exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_sub (double logx, double logy)
{
    return logx + log1p (-exp (logy - logx));
}




/* dpois_wrap (x_P_1,  lambda, g_log) ==
 *   dpois (x_P_1 - 1, lambda, g_log)
*/
static double dpois_wrap (double x_plus_1, double lambda, int give_log)
{
  
  //if (!R_FINITE(lambda))
  //return R_D__0;
    if (x_plus_1 > 1)
	return dpois_raw (x_plus_1 - 1, lambda, give_log);
    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
	return R_D_exp(-lambda - lgammafn(x_plus_1));
    else {
	double d = dpois_raw (x_plus_1, lambda, give_log);

	return give_log
	    ? d + log (x_plus_1 / lambda)
	    : d * (x_plus_1 / lambda);
    }
}

/*
 * Abramowitz and Stegun 6.5.29 [right]
 */
static double pgamma_smallx (double x, double alph, int lower_tail, int log_p)
{
    double sum = 0, c = alph, n = 0, term;

    /*
     * Relative to 6.5.29 all terms have been multiplied by alph
     * and the first, thus being 1, is omitted.
     */

    do {
	n++;
	c *= -x / n;
	term = c / (alph + n);
	sum += term;
    } while (fabs (term) > DBL_EPSILON * fabs (sum));


    if (lower_tail) {
	double f1 = log_p ? log1p (sum) : 1 + sum;
	double f2;
	if (alph > 1) {
	    f2 = dpois_raw (alph, x, log_p);
	    f2 = log_p ? f2 + x : f2 * exp (x);
	} else if (log_p)
	    f2 = alph * log (x) - lgamma1p (alph);
	else
	    f2 = pow (x, alph) / exp (lgamma1p (alph));

	return log_p ? f1 + f2 : f1 * f2;
    } else {
	double lf2 = alph * log (x) - lgamma1p (alph);

	if (log_p)
	    return R_Log1_Exp (log1p (sum) + lf2);
	else {
	    double f1m1 = sum;
	    double f2m1 = expm1 (lf2);
	    return -(f1m1 + f2m1 + f1m1 * f2m1);
	}
    }
} /* pgamma_smallx() */


/* Continued fraction for calculation of
 *    ???
 *  =  (i / d)  +  o(i/d)
 */
static double pd_lower_cf (double i, double d)
{
    double f = 0, of;

    double c1 = 0, c2, c3, c4;
    double a1 = 0, b1 = 1;
    double a2 = i, b2 = d;

#define	NEEDED_SCALE				\
	  (b2 > scalefactor) {			\
	    a1 /= scalefactor;			\
	    b1 /= scalefactor;			\
	    a2 /= scalefactor;			\
	    b2 /= scalefactor;			\
	}

#define max_it 200000

#ifdef DEBUG_p
    printf("pd_lower_cf(i=%.14g, d=%.14g)\n", i, d);
#endif

    while NEEDED_SCALE

    if(a2 == 0)
	return 0;/* when   d >>>> i  originally */

    c2 = a2;
    c4 = b2;

    while (c1 < max_it) {
	c1++;
	c2--;
	c3 = c1 * c2;
	c4 += 2;
	a1 = c4 * a2 + c3 * a1;
	b1 = c4 * b2 + c3 * b1;

	c1++;
	c2--;
	c3 = c1 * c2;
	c4 += 2;
	a2 = c4 * a1 + c3 * a2;
	b2 = c4 * b1 + c3 * b2;

	if NEEDED_SCALE

	if (b2 != 0) {
	    of = f;
	    f = a2 / b2;
	    /* convergence check: relative; absolute for small f : */
	    if (fabs (f - of) <= DBL_EPSILON * fmax2(1., fabs(f)))
		return f;
	}
    }

    printf(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n", f);
    return f;/* should not happen ... */
} /* pd_lower_cf() */
#undef NEEDED_SCALE


static double pd_lower_series (double lambda, double y)
{
    double term = 1, sum = 0;


    while (y >= 1 && term > sum * DBL_EPSILON) {
	term *= y / lambda;
	sum += term;
	y--;
    }
    /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
     *     =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n
     *     ~  y/lambda + o(y/lambda)
     */

    if (y != floor (y)) {
	/*
	 * The series does not converge as the terms start getting
	 * bigger (besides flipping sign) for y < -lambda.
	 */
	double f;

	f = pd_lower_cf (y, lambda + 1 - y);

	sum += term * f;
    }

    return sum;
} /* pd_lower_series() */

/*
 * Asymptotic expansion to calculate the probability that poisson variate
 * has value <= x.
 */
static double ppois_asymp (double x, double lambda, int lower_tail, int log_p)
{
    static const double coef15 = 1/12.;
    static const double coef25 = 1/288.;
    static const double coef35 = -139/51840.;
    static const double coef45 = -571/2488320.;
    static const double coef55 = 163879/209018880.;
    static const double coef65 =  5246819/75246796800.;
    static const double coef75 = -534703531/902961561600.;
    static const double coef1 = 2/3.;
    static const double coef2 = -4/135.;
    static const double coef3 = 8/2835.;
    static const double coef4 = 16/8505.;
    static const double coef5 = -8992/12629925.;
    static const double coef6 = -334144/492567075.;
    static const double coef7 = 698752/1477701225.;
    static const double two = 2;

    double dfm, pt_,s2pt,res1,res2,elfb,term;
    double ig2,ig3,ig4,ig5,ig6,ig7,ig25,ig35,ig45,ig55,ig65,ig75;
    double f, np, nd;

    dfm = lambda - x;
    pt_ = -x * log1pmx (dfm / x);
    s2pt = sqrt (2 * pt_);
    if (dfm < 0) s2pt = -s2pt;

    ig2 = 1.0 + pt_;
    term = pt_ * pt_ * 0.5;
    ig3 = ig2 + term;
    term *= pt_ / 3;
    ig4 = ig3 + term;
    term *= pt_ / 4;
    ig5 = ig4 + term;
    term *= pt_ / 5;
    ig6 = ig5 + term;
    term *= pt_ / 6;
    ig7 = ig6 + term;

    term = pt_ * (two / 3);
    ig25 = 1.0 + term;
    term *= pt_ * (two / 5);
    ig35 = ig25 + term;
    term *= pt_ * (two / 7);
    ig45 = ig35 + term;
    term *= pt_ * (two / 9);
    ig55 = ig45 + term;
    term *= pt_ * (two / 11);
    ig65 = ig55 + term;
    term *= pt_ * (two / 13);
    ig75 = ig65 + term;

    elfb = ((((((coef75/x + coef65)/x + coef55)/x + coef45)/x + coef35)/x +
	     coef25)/x + coef15) + x;
    res1 = ((((((ig7*coef7/x + ig6*coef6)/x + ig5*coef5)/x + ig4*coef4)/x +
	      ig3*coef3)/x + ig2*coef2)/x + coef1)*sqrt(x);
    res2 = ((((((ig75*coef75/x + ig65*coef65)/x + ig55*coef55)/x + ig45*coef45)/
	      x + ig35*coef35)/x + ig25*coef25)/x + coef15)*s2pt;

    if (!lower_tail) elfb = -elfb;
    f = (res1 + res2) / elfb;

    np = pnorm (s2pt, 0.0, 1.0, !lower_tail, log_p);
    nd = dnorm (s2pt, 0.0, 1.0, log_p);

#ifdef DEBUG_p
    printf ("pp*_asymp(): f=%.14g np=%.14g nd=%.14g  f*nd=%.14g\n",
	      f, np, nd, f * nd);
#endif

    if (log_p)
	return (f >= 0)
	    ? logspace_add (np, log (fabs (f)) + nd)
	    : logspace_sub (np, log (fabs (f)) + nd);
    else
	return np + f * nd;
} /* ppois_asymp() */




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

    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
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

    if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

	/* The answer is less than half precision because
	 * the argument is too near a negative integer. */

	ML_ERROR(ME_PRECISION);
    }

    return ans;
}
