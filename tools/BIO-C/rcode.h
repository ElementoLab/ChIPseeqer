
#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi)) */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2)) */
#endif



#define DBL_EPSILON 2.2204460492503131e-16
#define DBL_MIN 2.2250738585072014e-308
#define DBL_MAX 1.7976931348623157e+308
#define FLT_RADIX 2
#define FLT_MANT_DIG 24
#define FLT_MIN_EXP (-125)
#define FLT_MAX_EXP 128
#define DBL_MANT_DIG 53
#define DBL_MIN_EXP (-1021)
#define DBL_MAX_EXP 1024

#define MATHLIB_ERROR(fmt,x)	{ printf(fmt,x); exit(1); }
#define MATHLIB_WARNING(fmt,x)		printf(fmt,x)
#define MATHLIB_WARNING2(fmt,x,x2)	printf(fmt,x,x2)
#define MATHLIB_WARNING3(fmt,x,x2,x3)	printf(fmt,x,x2,x3)
#define MATHLIB_WARNING4(fmt,x,x2,x3,x4) printf(fmt,x,x2,x3,x4)

#define ISNAN(x) (isnan(x)!=0)
#define R_FINITE(x)    finite(x)

#define ML_POSINF	(1.0 / 0.0)
#define ML_NEGINF	((-1.0) / 0.0)
#define ML_NAN		(0.0 / 0.0)


#define ML_ERROR(x)	/* nothing */
#define ML_UNDERFLOW	(DBL_MIN * DBL_MIN)
#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN); return ML_NAN; }

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

/* Wilcoxon Signed Rank Distribution */

#define SIGNRANK_MAX 50

/* Formerly private part of Mathlib.h */

/* always remap internal functions */
#define bd0       	Rf_bd0
#define chebyshev_eval	Rf_chebyshev_eval
#define chebyshev_init	Rf_chebyshev_init
#define i1mach		Rf_i1mach
#define gammalims	Rf_gammalims
#define lfastchoose	Rf_lfastchoose
#define lgammacor	Rf_lgammacor
#define stirlerr       	Rf_stirlerr

	/* Chebyshev Series */

int	chebyshev_init(double*, int, double);
double	chebyshev_eval(double, const double *, const int);

	/* Gamma and Related Functions */

void	gammalims(double*, double*);
double	lgammacor(double); /* log(gamma) correction */
double  stirlerr(double);  /* Stirling expansion "error" */

double	lfastchoose(double, double);

double  bd0(double, double);

/* Consider adding these two to the API (Rmath.h): */
double	dbinom_raw(double, double, double, double, int);
double	dpois_raw (double, double, int);
double  pnchisq_raw(double, double, double, double, double, int);
double  pgamma_raw(double, double, int, int);

int	i1mach(int);

/* From toms708.c */
void bratio(double a, double b, double x, double y,
	    double *w, double *w1, int *ierr);





#define give_log log_p
                                                        /* "DEFAULT" */
                                                        /* --------- */
#define R_D__0  (log_p ? ML_NEGINF : 0.)                /* 0 */
#define R_D__1  (log_p ? 0. : 1.)                       /* 1 */
#define R_DT_0  (lower_tail ? R_D__0 : R_D__1)          /* 0 */
#define R_DT_1  (lower_tail ? R_D__1 : R_D__0)          /* 1 */

#define R_D_Lval(p)     (lower_tail ? (p) : (1 - (p)))  /*  p  */
#define R_D_Cval(p)     (lower_tail ? (1 - (p)) : (p))  /*  1 - p */

#define R_D_val(x)      (log_p  ? log(x) : (x))         /*  x  in pF(x,..) */

#define R_D_qIv(p)      (log_p  ? exp(p) : (p))         /*  p  in qF(p,..) */
#define R_D_exp(x)      (log_p  ?  (x)   : exp(x))      /* exp(x) */
#define R_D_log(p)      (log_p  ?  (p)   : log(p))      /* log(p) */
#define R_D_Clog(p)     (log_p  ? log1p(-(p)) : (1 - (p)))/* [log](1-p) */

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

#define R_DT_val(x)     (lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)    (lower_tail ? R_D_Clog(x) : R_D_val(x))


#define R_DT_qIv(p)     (log_p ? (lower_tail ? exp(p) : - expm1(p)) \
                               : R_D_Lval(p))

#define R_DT_CIv(p)     (log_p ? (lower_tail ? -expm1(p) : exp(p)) \
                               : R_D_Cval(p))
#define R_D_nonint(x)      (fabs((x) - floor((x)+0.5)) > 1e-7)
#define R_D_forceint(x)   floor((x) + 0.5)
