bin_PROGRAMS = genregexp mi_find mi_optimize mi_dist mi_signif mi_optimize_slow binom_test_greater recompare mi_combine myhypergeom MyCompareAce mi_motif_compare mi_optimize_motif_WM

# includes and libraries for all
AM_CPPFLAGS = -I`pwd`/../modules/include -g

if HAVE_SUNMATH
LDADD_SUNMATH = -lsunmath
else
LDADD_SUNMATH = 
endif

LDADD = libfire.a -L`pwd`/../modules/lib -lpcre -lm $(LDADD_SUNMATH)

# programs
mi_find_SOURCES = mi_find.c

mi_optimize_SOURCES = mi_optimize.c

mi_optimize_slow_SOURCES = mi_optimize_slow.c

mi_signif_SOURCES = mi_signif.c

mi_dist_SOURCES = mi_dist.c

mi_combine_SOURCES = mi_combine.c

mi_motif_compare_SOURCES = mi_motif_compare.c

mi_optimize_motif_WM_SOURCES = mi_optimize_motif_WM.c

binom_test_greater_SOURCES = binom_test_greater.c beta.c toms708.c d1mach.c i1mach.c rcode.h beta.h 

genregexp_SOURCES = genregexp.c
recompare_SOURCES = recompare.c

myhypergeom_SOURCES = myhypergeom.c
myhypergeom_LDADD = -lm

MyCompareAce_SOURCES = MyCompareAce.c
MyCompareAce_LDADD = -lm

noinst_LIBRARIES = libfire.a
libfire_a_SOURCES = statistics.h statistics.c dataio.h dataio.c prefix.c prefix.h information.c information.h mi_library.c mi_library.h sequences.c sequences.h regexp.c regexp.h

# to include in distribution
EXTRA_DIST = pcre-7.4
