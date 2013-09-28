## Process this file to produce a configure script:
##   aclocal -I build-aux/m4 && autoheader && autoconf && automake
##

## Preamble - used to set up meta paths, meta-information, etc.
#
# require a minimum version of AutoConf
AC_PREREQ([2.63])
# software name, version, contact address
AC_INIT([geotop],[1.2.10],[stefano.end@gmail.com])
# if this file does not exist, `configure` was invoked in the wrong directory
AC_CONFIG_SRCDIR([geotop/geotop.cc])
# directories (relative to top-level) to look into for AutoConf/AutoMake files
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux])
# enable AutoMake
AM_INIT_AUTOMAKE([1.10])
# all defined C macros (HAVE_*) will be saved to this file
AC_CONFIG_HEADERS([config.h])

## Initialize GNU LibTool
#
# GNU LibTool provides a portable way to build libraries.  AutoMake
# knows how to use it; you just need to activate it.
LT_INIT([static])

## Option to compile geotop statically ##
AC_ARG_WITH([static-geotop],
    [AS_HELP_STRING([--with-static-geotop],
            [Compile a static executable.] Default: no)],
    [with_static_geotop=yes],
    [])

## Enable debug compilation ##
AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],
               [enable debugging, setting -g3 and -O0. Default: no]),
[case "${enableval}" in
             yes) debug=true;;
             no)  debug=false ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],
[debug=false])

# This test has to run *before* AC_PROC_CC
if test x"$debug" = x"true"; then
   AC_DEFINE([DEBUG], [1], [Enable debug compilation])
   CFLAGS=${CFLAGS:="-g3 -O0"}
   CXXFLAGS=${CXXFLAGS:="-g3 -O0"}
fi

AM_CONDITIONAL(DEBUG, test x"$debug" = x"true")
## end debug compilation option ##

AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_PROG_INSTALL


## Use internal interpolation or MeteoIO one ##
AC_ARG_ENABLE([internal-meteodistr],
    [AS_HELP_STRING([--enable-internal-meteodistr],
        [Use internal meteo interpolation, routine instead of using MeteoIO])],
    [AC_DEFINE([USE_INTERNAL_METEODISTR], [1], [Decide which interpolation to use])],
    [])

if test x"$enable_internal_meteodistr" = x"yes"
then
    METEODISTR_RESULTS_DIR=results.meteodistr
else
    METEODISTR_RESULTS_DIR=results
fi
AC_SUBST([METEODISTR_RESULTS_DIR])

## Enable NetCDF ##
AC_ARG_ENABLE([netcdf],
    [AS_HELP_STRING([--enable-netcdf],
        [Use NetCDF for ...?])],
    [AC_DEFINE([USE_NETCDF], [1], [Use NetCDF])],
    [])

AM_CONDITIONAL([INCLUDE_NETCDF_SUBDIR], test "x$enable_netcdf" = "xyes")

## If netcdf is enabled, check if it's actually available
AS_IF([test "x$enable_netcdf" = "xyes"],
    [AC_CHECK_LIB([netcdf], [nc_inq_libvers],
        [], [AC_MSG_ERROR([No NetCDF library found])])])
## end netcdf check

## Configure meteoio ##
AC_ARG_WITH([meteoio],
    [AS_HELP_STRING([--with-meteoio],
        [location of the MeteoIO library])],
    [METEOIO_CPPFLAGS="-I$withval/include" METEOIO_LDFLAGS="-L$withval/lib/"],
    [])

AC_CHECK_FILE($with_meteoio/lib/libmeteoio.a,
        [HAVE_STATIC_LIBMETEOIO=1],
        [HAVE_STATIC_LIBMETEOIO=0],
)

AC_CHECK_FILE($with_meteoio/lib/libmeteoio.a,
        [HAVE_STATIC_LIBMETEOIO=1],
        [HAVE_STATIC_LIBMETEOIO=0],
)

# Configure meteoio linker flags
METEOIO_LIBS="-lmeteoio"

# There *has* to be a more-autoconf way to check if
# --with-static-geotop _and_ if HAVE_STATIC_LIBMETEOIO=0
if test x$with_static_geotop = xyes; then
   if test $HAVE_STATIC_LIBMETEOIO = 0; then
     AC_MSG_ERROR([Missing MeteoIO static library. Specify the correct path with `--with-meteoio` and ensure that MeteoIO has been compiled with STATIC library support])
   else
      METEOIO_LIBS="$with_meteoio/lib/libmeteoio.a"
      # Very ugly hack to find libproj!
      for path in /usr/lib /usr/local/lib /opt/local/lib
      do
          if test -f $path/libproj.a; then
            METEOIO_LIBS+=" $path/libproj.a"
            break
          fi
      done
      METEOIO_LIBS+=" -lpthread"
   fi
fi

AC_SUBST([METEOIO_CPPFLAGS])
AC_SUBST([METEOIO_LDFLAGS])
AC_SUBST([METEOIO_LIBS])

CPPFLAGS+=" $METEOIO_CPPFLAGS"
LDFLAGS+=" $METEOIO_LDFLAGS"

#### Check MeteoIO library ####
# save old values for compiling
OLD_LIBS="$LIBS"
LIBS+=" $METEOIO_LIBS"

# Set language
AC_LANG([C++])
# Test if we can compile a program
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include <meteoio/MeteoIO.h>
    ]],
        [[mio::getLibVersion();]])],
    [AC_MSG_RESULT([yes])
        FOUND_METEOIO=1;],
    [AC_MSG_RESULT([no])
        FOUND_METEOIO=0;])

# reset original *FLAGS
LIBS="$OLD_LIBS"

if test $FOUND_METEOIO != 1; then
    AC_MSG_NOTICE([])
    AC_MSG_NOTICE([The MeteoIO library was not found!])
    AC_MSG_NOTICE([])
    AC_MSG_ERROR([ Please provide the path to the MeteoIO installation using `--with-meteoio` and re-run the configure script.])
fi

#### END checking of MeteoIO ####


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
                 geotop/Makefile
                 libraries/Makefile
                 libraries/ascii/Makefile
                 libraries/fluidturtle/Makefile
                 wrf/Makefile
                 netCDF/Makefile
                 meteoio_plugin/Makefile
                 tests/test_sample_run/test_runner.py
                 ])
AC_OUTPUT
