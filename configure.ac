#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.55])
AC_INIT([Semi-Markov], [0.1], [ajd27@cornell.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_SRCDIR([src/stochnet.h])

AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX([clang++ g++])
AX_BOOST_BASE([1.54])
# Find the flavor of the Boost library.
BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/-L\(.\+\)/\1/'`

AC_CONFIG_FILES([Makefile])


AC_ARG_ENABLE([variant],
  [AS_HELP_STRING([--enable-variant],
    [Build with the boost library that ends in -mt])],
  [boost_variant="$enableval"],
  [boost_variant="yes"])


AC_ARG_ENABLE([prand],
  [AS_HELP_STRING([--enable-prand],
    [Build the prand wrapper program @<:@no@:>@])],
  [build_prand="$enableval"],
  [build_prand="yes"])
AM_CONDITIONAL([BUILD_PRAND], [test x"$build_PRAND" != xno])

AC_CHECK_LIB([pthread], [exit])

# Boost libraries can end in -mt, bare, or other endings.
FLAVOR=$boost_variant
LDFLAGS_SAVE=$LDFLAGS
LDFLAGS=$BOOST_LDFLAGS $LDFLAGS
AC_CHECK_LIB([boost_system${FLAVOR}], [exit])
AC_CHECK_LIB([boost_random${FLAVOR}], [exit])
AC_CHECK_LIB([boost_program_options${FLAVOR}], [exit])
AC_CHECK_LIB([boost_filesystem${FLAVOR}], [exit])
AC_CHECK_LIB([boost_date_time${FLAVOR}], [exit])
AC_CHECK_LIB([boost_log${FLAVOR}], [exit])
AC_CHECK_LIB([boost_log_setup${FLAVOR}], [exit])
LDFLAGS=$LDFLAGS_SAVE

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.

AC_OUTPUT
