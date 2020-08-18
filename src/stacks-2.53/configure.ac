# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Stacks], [2.53])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign parallel-tests subdir-objects])
AC_CONFIG_SRCDIR([src/ustacks.cc])
AC_CONFIG_SRCDIR([htslib/hts.c])
AC_CONFIG_HEADERS([config.h])
m4_pattern_allow([AC_OPENMP])

#
# Get BAM, aka HTSLib library and include locations if requested
#
#AC_ARG_ENABLE([bam],
#    AS_HELP_STRING([--disable-bam], [Disable use of BAM files through HTSLib.]))
#AS_IF([test "x$enable_bam" != "xno"], [
#AC_DEFINE([HAVE_BAM], [1], [Enable compilation with Samtools BAM library])
#BAM_CFLAGS='-I./htslib/htslib'
#AC_SUBST([BAM_CFLAGS])
#BAM_LIBS='./htslib/libhts.a'
#AC_SUBST([BAM_LIBS])
#])

# Checks for programs.
AC_PROG_CXX
AM_PROG_CC_C_O
AX_CXX_COMPILE_STDCXX(11,, [mandatory])

# Checks for ar-based static libs.
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Checks for libraries.
AC_CHECK_LIB([gomp], [omp_set_num_threads],, [AC_MSG_WARN([Unable to locate OpenMP library, you should probably specify '--disable-openmp'.])])
AC_CHECK_LIB([z], [gzread],, [AC_MSG_ERROR([Zlib not found, reading gzipped files will not be possible.])]) 

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h])
AC_CHECK_HEADERS([unistd.h])

# Check for OpenMP parallel execution support
AC_OPENMP

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_CLOSEDIR_VOID
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt])
AC_CHECK_FUNCS([clock_gettime])

# `clock_gettime()`, declared in unistd.h, may be in librt (e.g. for glibc).
if test $ac_cv_func_clock_gettime = no; then
	AC_CHECK_LIB(rt, clock_gettime,
		[LIBS="-lrt $LIBS"
		AC_DEFINE(HAVE_CLOCK_GETTIME, 1, [Define to 1 if you have the `clock_gettime' function.])
		])
fi

#
# Test if regular expressions are functional. Early regex implementations have nonfunctional stubs,
# so the check for C++11 succeeds, but the compiled code will throw an exception.
#
AC_LANG(C++)
AC_DEFUN([AC_REGEX_FUNC], [
AC_CACHE_CHECK([for functional regular expressions], [stacks_cv_header_regex_func],
[
stacks_cv_header_regex_func=no
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <regex>
int main() {
  const std::regex regex(".*");
  const std::string string = "This should match!";
  const auto result = std::regex_search(string, regex);
  return result ? 0 : 1;
}
]])],
    [stacks_cv_header_regex_func=yes],
    [AC_MSG_ERROR([Regular expressions are not functional, you need g++ 4.9.0 or greater.])],
    [AC_MSG_WARN([Requires g++ 4.9.0 or greater.])]
)])
if test "$stacks_cv_header_regex_func" = yes; then
    AC_DEFINE([REGEX_FUNCTIONAL], [1], [Regular expressions are defined and implemented])
fi
])

AC_REGEX_FUNC

# For test harness
AC_PROG_AWK

AC_CONFIG_FILES([Makefile htslib/Makefile])
AC_OUTPUT
