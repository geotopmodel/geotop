## Process this file to produce a configure script:
##   aclocal -I build-aux/m4 && autoheader && autoconf && automake
##

## Preamble - used to set up meta paths, meta-information, etc.
#
# require a minimum version of AutoConf
AC_PREREQ([2.65])
# software name, version, contact address
AC_INIT([geotop],[1.2.10],[stefano.end@gmail.com])
# if this file does not exist, `configure` was invoked in the wrong directory
AC_CONFIG_SRCDIR([src/geotop/geotop.c])
# directories (relative to top-level) to look into for AutoConf/AutoMake files
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux])
# enable AutoMake
AM_INIT_AUTOMAKE([1.10])
# all defined C macros (HAVE_*) will be saved to this file
AC_CONFIG_HEADERS([config.h])

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_INSTALL

## Initialize GNU LibTool
#
# GNU LibTool provides a portable way to build libraries.  AutoMake
# knows how to use it; you just need to activate it.
LT_INIT([static])


AC_ARG_WITH([static-geotop],
    [AS_HELP_STRING([--with-static-geotop=yes|no],
            [Compile a static executable.])],
    [with_static_geotop=yes],
    [])
# The `AM_CONDITIONAL` macro makes the results of the test available into `Makefile.am`
AM_CONDITIONAL([WITH_STATIC_GEOTOP], [test "_$with_static_geotop" = '_yes'])


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

# Substitute all conditionals in these files; this is normally used to
# create `Makefile`s but could also be used for scripts, include
# files, etc.
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/geotop/Makefile
                 src/libraries/Makefile
                 src/libraries/ascii/Makefile
                 src/libraries/math/Makefile
                 src/libraries/fluidturtle//Makefile
                 src/libraries/geomorphology/Makefile])
AC_OUTPUT
