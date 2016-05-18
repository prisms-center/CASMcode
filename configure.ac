#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([autocasm], [0.0.0], [casm@casm.edu])
AC_CONFIG_SRCDIR([include/casm/CASM_global_definitions.hh])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign subdir-objects])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgzstream':
AC_CHECK_LIB([gzstream], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])

#I added this
LT_INIT
AX_CXX_COMPILE_STDCXX_11
INCLUDES= -src
AX_BOOST_FILESYSTEM
AX_BOOST_SYSTEM
AX_BOOST_REGEX
AX_BOOST_PROGRAM_OPTIONS

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_HEADER_MAJOR
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor iswprint memmove memset mkdir pow rmdir select sqrt strchr strpbrk strrchr strstr strtol])

AC_CONFIG_FILES([Makefile
                 apps/auto/Makefile
                 src/casm/Makefile])
AC_OUTPUT
