#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.10)
AC_INIT(FIRE, 1.1, elemento@princeton.edu, FIRE)

AM_CONFIG_HEADER([PROGRAMS/config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])


AC_CANONICAL_HOST

# m4_pattern_allow
# Checks for programs.
# AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# searching for essential programs
AC_PATH_PROGS(FIND, find)
AC_PATH_PROGS(PERL, perl)
AC_PATH_PROGS(CP, cp)

AC_HAVE_LIBRARY([pcre], [PCRELIB=-lpcre], [PCRELIB=])
AC_SUBST(PCRELIB)

AC_SUBST(host_os)
case "$host_os" in
	*cygwin* ) AC_CHECK_FILE([$(pwd)/modules/lib/libpcre.dll.a], [PCRELOCALLIB=yes], [PCRELOCALLIB=])
		   ;;	
	*darwin* ) AC_CHECK_FILE([$(pwd)/modules/lib/libpcre.dylib], [PCRELOCALLIB=yes], [PCRELOCALLIB=])
		   ;;
	*) AC_CHECK_FILE([$(pwd)/modules/lib/libpcre.so], [PCRELOCALLIB=yes], [PCRELOCALLIB=])
esac
AC_SUBST(PCRELOCALLIB)



 
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stddef.h stdlib.h string.h sys/timeb.h])

# added for Solaris
AC_CHECK_HEADERS(sunmath.h, HAVE_SUNMATH=yes)
AM_CONDITIONAL([HAVE_SUNMATH], [test x"$HAVE_SUNMATH" != x])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor ftime pow sqrt strdup trunc nearbyintf])

# 
AC_ARG_ENABLE([myhypergeom], ["uses myhypergeom instead of Perl module"], [USE_MYHYPERGEOM=yes])
# echo myhp=$USE_MYHYPERGEOM
AM_CONDITIONAL([USE_MYHYPERGEOM], [test "$USE_MYHYPERGEOM" = yes])

# AC_PREFIX_DEFAULT([$(PWD)/modules])
AC_CONFIG_FILES([Makefile
                 PROGRAMS/Makefile
                 SCRIPTS/Makefile])
AC_OUTPUT

if test x"$PCRELIB" = x -a x"$PCRELOCALLIB" = x; then
  echo ""
  echo "Attention: you don't appear to have libpcre installed on your machine."
  echo "You can either install it yourself from http://www.pcre.org/"
  echo "Or execute 'sh install-pcre.sh' to install the pcre-7.4 library we provide with FIRE."
  exit -1	
fi
