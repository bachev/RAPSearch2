#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([RAPSearch2], [2.23])
AC_CONFIG_SRCDIR([Src/main.cpp])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([-Wall foreign])

# init compile flags from environment
CXXFLAGS="${CXXFLAGS} "
# get the system name (for later)
systemuname=`uname`


AC_ARG_ENABLE(clang,
	AC_HELP_STRING([--enable-clang],
	[use clang instead of gcc]))

if test "${enable_clang}" = "yes" ; then
  CC=clang
fi

AC_ARG_ENABLE(static,
	AC_HELP_STRING([--enable-static],
	[enables building of static binaries, default is no]))



AC_LANG([C++])

# ---------------------------------------------------------------------------
# Check compiler family
# ---------------------------------------------------------------------------
bccheck1="[gcc]"
bccheck2="[g++]"

case "${CC}" in
*gcc*)
  bccheck1="[gcc]"
  bccheck2="[g++]"
;;
*clang*)
  bccheck1="[clang]"
  bccheck2="[clang++]"
;;
esac

AC_PROG_CC($bccheck1)
AC_PROG_CXX($bccheck2)

bachCCOMPILER=""
case "${CXX}" in
*g++*)
  bachCCOMPILER=gcc
;;
*clang++*)
  bachCCOMPILER=clang
;;
*)
  bachCCOMPILER=unknown
;;
esac

if test "$bachCCOMPILER" != "gcc" ; then
  # No, autoconf, no! Bad autoconf!
  # We do not have GCC! I'm telling you!
  # (though, to be fair, it's not entirely your fault)
  GCC="no"
fi

# End compiler family
# ---------------------------------------------------------------------------


# Checks for programs.

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create], [], [AC_MSG_ERROR(pthread not found.)])
AC_CHECK_LIB([boost_serialization], [main], [], [AC_MSG_ERROR(boost_serialization not found.)])
AC_CHECK_LIB([boost_system], [main], [], [AC_MSG_ERROR(boost_system not found.)])
AC_CHECK_LIB([boost_thread], [main], [], [AC_MSG_ERROR(boost_thread not found.)])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h])
AC_CHECK_HEADERS([boost/archive/binary_iarchive.hpp boost/archive/binary_oarchive.hpp boost/bind.hpp boost/lambda/lambda.hpp boost/lexical_cast.hpp boost/serialization/vector.hpp boost/thread/mutex.hpp boost/thread/thread.hpp])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset pow sqrt])


# ---------------------------------------------------------------------------
# Set optimisation flags
# ---------------------------------------------------------------------------
if test "${enable_optimisations+set}" != "set"; then
  enable_optimisations="yes"
fi
AC_ARG_ENABLE(optimisations,
	AC_HELP_STRING([--enable-optimisations],
	[enables compiler optimisations, default is yes]))

AC_ARG_ENABLE(native,
	AC_HELP_STRING([--enable-native],
	[compile only for your specific processor, default is no]))

if test "${enable_native+set}" != "set"; then
  enable_native="no"
fi


if test "${enable_optimisations}" = "yes" ; then
  AC_MSG_CHECKING(for compiler optimisation settings)

  if test $bachCCOMPILER = gcc; then
    # For GCC, be conservative for all plattforms except the ones where
    #  -O3 is known to work

    OPTFLAGS="-O2"
    if test "${systemuname}" = "Linux"; then
      OPTFLAGS="-O3 -funroll-loops"
    fi

    if test "${systemuname}" = "Darwin"; then
      OPTFLAGS="-O3 -funroll-loops"
    fi

    # native?
    if test "${enable_native}" = "yes" ; then
      OPTFLAGS+=" -march=native -mtune=native"
    fi

  else
    # default
    OPTFLAGS="-O2"
  fi

  CXXFLAGS="${CXXFLAGS} ${OPTFLAGS}"

  AC_MSG_RESULT(${OPTFLAGS})
fi
# End of optimisation tests
# ---------------------------------------------------------------------------



if test "${enable_static}" = "yes" ; then
  CXXFLAGS="--static $CXXFLAGS"
fi


AC_CONFIG_FILES([Makefile
                 Src/Makefile])
AC_OUTPUT
