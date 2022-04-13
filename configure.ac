#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([qmckl_bench], [0.0.1], [])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([include/config.h])
AM_INIT_AUTOMAKE([subdir-objects color-tests parallel-tests silent-rules 1.11])

AM_PROG_AR
AM_MAINTAINER_MODE()

# Checks for programs.
AC_LANG(C)
AC_PROG_F77
AC_PROG_CC
AC_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_LN_S

AC_ARG_WITH(qmckl, [AS_HELP_STRING([--with-qmckl],[Specify installation directory of QMCKL])],
            with_qmckl=$withval, with_qmckl=yes)

# Checks for libraries.
## BLAS
AX_BLAS([], [AC_MSG_ERROR([BLAS was not found.])])

## LAPACK
AX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found.])])

# Checks for library functions.
AC_CHECK_LIB([m], [sqrt])
AS_IF([test "$with_qmckl" = "yes"],
  [
    AC_CHECK_LIB([qmckl], [qmckl_context_create], [],
      [ AC_MSG_ERROR([QMCKL not found]) ]
    )
  ], [
    QMCKL_LIBS="${with_qmckl}/lib/libqmckl.a"
    QMCKL_CFLAGS="-I${with_qmckl}/include"
    QMCKL_LDFLAGS=""
  ])

LIBS="${QMCKL_LIBS} ${LIBS}"
LDFLAGS="${QMCKL_LDFLAGS} ${LDFLAGS}"
CFLAGS="${QMCKL_CFLAGS} ${CFLAGS} "
AC_DEFINE_UNQUOTED(DATA_DIR, ["`cd ${srcdir}/data && echo $(pwd -P)`"], [Path to the data files])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h qmckl.h])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
