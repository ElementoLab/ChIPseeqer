#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#include "hypergeom/hypergeom.h"

#include "const-c.inc"

MODULE = Hypergeom		PACKAGE = Hypergeom		

INCLUDE: const-xs.inc

double
cumhyper(i,s1,s2,N)
	int             i
	int             s1
	int             s2
	int             N
     OUTPUT:
	RETVAL

double
lcumhyper(i,s1,s2,N)
	int             i
	int             s1
	int             s2
	int             N
     OUTPUT:
	RETVAL

double
cumbino(k,N,p)
	int		k
	int		N
	double		p
     OUTPUT:
	RETVAL

double
bino(x,N,p)
	int		x
	int		N
	double		p
     OUTPUT:
	RETVAL

