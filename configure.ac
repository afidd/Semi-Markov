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
AC_CONFIG_FILES([Makefile])
#AC_PROG_LIBTOOL

# Boost libraries can end in -mt, bare, or other endings.
FLAVOR=-mt
AC_CHECK_LIB([pthread], [exit])
AC_CHECK_LIB([boost_system${FLAVOR}], [exit])
AC_CHECK_LIB([boost_random${FLAVOR}], [exit])
AC_CHECK_LIB([boost_program_options${FLAVOR}], [exit])
AC_CHECK_LIB([boost_filesystem${FLAVOR}], [exit])
AC_CHECK_LIB([boost_date_time${FLAVOR}], [exit])
AC_CHECK_LIB([boost_log${FLAVOR}], [exit])
AC_CHECK_LIB([boost_log_setup${FLAVOR}], [exit])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.

AC_OUTPUT
