#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Copter], [0.8.7], [jwgcarlson AT berkeley.edu])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE([disable])
LT_PREREQ([2.2])
LT_INIT
AC_CONFIG_SRCDIR([src/Cosmology.h])
AC_CONFIG_MACRO_DIR([m4])

# Decide if this is a debug or optimized build.
AC_ARG_ENABLE(
    debug,
    [AS_HELP_STRING([--enable-debug], [compile with debugging symbols (default no)])]
)

# Set default compiler flags (if none were set by the user).
if test x"$enable_debug" = xyes; then
    # Debugging build
    AC_DEFINE([DEBUG], 1, [Enable debugging output])
    CFLAGS=${CFLAGS:--g -O2}
    CXXFLAGS=${CXXFLAGS:-$CFLAGS}
else
    enable_debug=no
    # Optimized build
    AX_CC_MAXOPT

    # This is broken for autoconf 2.64 and up, and I don't understand why.  It
    # gives errors about AC_PROG_CXX being expanded before required, even
    # though AX_CXX_MAXOPT is identical to AX_CC_MAXOPT...
    #AX_CXX_MAXOPT
    CXXFLAGS=${CXXFLAGS-$CFLAGS}
fi

# Decide whether or not to keep assert statements.
AC_ARG_ENABLE(
    assert,
    [AS_HELP_STRING([--disable-assert], [ignore assert statements])]
)
if test x"$enable_assert" = xno; then
    AC_DEFINE([NDEBUG], 1, [Ignore assert statements])
fi

# Decide whether or not to display verbose runtime messages.
AC_ARG_ENABLE(
    verbose,
    [AS_HELP_STRING([--enable-verbose], [display verbose runtime messages (default no)])]
)
if test x"$enable_verbose" = xyes; then
    AC_DEFINE([VERBOSE], [1], [Enable verbose output]),
else
    enable_verbose=no
fi

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CXXCPP

# Checks for common header files.
#AC_CHECK_HEADERS([limits.h pthread.h sys/stat.h unistd.h])
AC_CHECK_HEADERS([limits.h sys/stat.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_C_CONST
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_C_VOLATILE
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STAT
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([memset strerror strtol pow sqrt])

# Check for GSL.
#AX_PATH_GSL(1.10)
#CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"
#LIBS="$LIBS $GSL_LIBS"

# Check for FFTW 3.
AC_CHECK_HEADER([fftw3.h])
AC_CHECK_LIB(fftw3, fftw_execute)

# Check for OpenMP support.
AC_OPENMP

# Check for boost.
#acx_have_boost=no
#AX_BOOST_BASE([], [acx_have_boost=yes], [])
#if test x"$acx_have_boost" != xyes; then
#    AC_MSG_ERROR([could not find boost])
#fi
#AM_CONDITIONAL([HAVE_AM_BOOST], [test $acx_have_boost = yes])

# Decide whether to build Python bindings.
AC_ARG_ENABLE(
    [python],
    [AS_HELP_STRING([--enable-python], [build Python bindings (default no)])],
    ,
    enable_python=no
)

AC_ARG_WITH(
    [swig],
    [AS_HELP_STRING([--with-swig=PATH], [path to SWIG executable])]
)

if test "$enable_python" = yes; then
    AM_PATH_PYTHON([2.3])

    # Allow user to override default search path for SWIG
    if test -n "$with_swig"; then
        if test -x "$with_swig"; then
            SWIG="$with_swig"
        else
            AC_MSG_ERROR([swig path "$with_swig" does not exist or is not executable])
        fi
    fi

    # Find SWIG executable and Python development headers
    AC_PROG_SWIG([1.3.35])
    SWIG_ENABLE_CXX
    SWIG_PYTHON

    # Enable compilation of Python bindings
    AC_SUBST([PYCOPTER_PY], [pycopter.py])
    AC_SUBST([_PYCOPTER_LA], [_pycopter.la])
fi

# Get the default number of threads to use for parallelized code.
#AC_ARG_WITH(
#    [numthreads],
#    [AS_HELP_STRING([--with-numthreads=N],
#                    [number of threads to use for parallelized code])],
#    [OMP_NUM_THREADS="$with_numthreads"]
#)
#AC_SUBST(OMP_NUM_THREADS)

# Bypass unnecessary invocations of autoheader.  (Beware of hack!)
#AUTOHEADER="echo autoheader ignored"

# Write Makefiles and configuration header.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 data/Makefile
                 doc/Makefile
                 examples/Makefile
                 python/Makefile
                 src/Makefile
                 tests/Makefile])
AC_OUTPUT

# Summarize configuration.
cat << EOF
-------------------------------------------------------------------------------
Copter configuration summary:

Debug?  $enable_debug
Verbose?  $enable_verbose
Python bindings?  $enable_python

Installation prefix = $prefix
CC = $CC
CXX = $CXX
CPPFLAGS = $CPPFLAGS
CFLAGS = $CFLAGS
CXXFLAGS = $CXXFLAGS
LDFLAGS = $LDFLAGS
LIBS = $LIBS
OPENMP_CFLAGS = $OPENMP_CFLAGS
SWIG = $SWIG
PYTHON_CPPFLAGS = $PYTHON_CPPFLAGS
PYTHON_LDFLAGS = $PYTHON_LDFLAGS
-------------------------------------------------------------------------------
EOF
