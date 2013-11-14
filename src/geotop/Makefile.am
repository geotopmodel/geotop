bin_PROGRAMS = geotop

geotop_SOURCES = blowingsnow.c channels.c clouds.c deallocate.c energy.balance.c geotop.c indices.c input.c meteo.c meteodata.c meteodistr.c output.c parameters.c PBSM.c pedo.funct.c radiation.c recovering.c snow.c tables.c times.c turbulence.c vegetation.c water.balance.c transport.model.c

AM_CPPFLAGS = -I$(srcdir) \
              -I$(top_srcdir)/src/libraries/ascii \
              -I$(top_srcdir)/src/libraries/math \
              -I$(top_srcdir)/src/libraries/fluidturtle \
              -I$(top_srcdir)/src/libraries/geomorphology
AM_CFLAGS =

AM_LDFLAGS = -static

if WITH_STATIC_GEOTOP
    AM_LDFLAGS += -all-static
endif

geotop_LDADD = $(top_srcdir)/src/libraries/ascii/libascii.la \
               $(top_srcdir)/src/libraries/math/libmath.la \
               $(top_srcdir)/src/libraries/geomorphology/libgeomorphology.la \
               $(top_srcdir)/src/libraries/fluidturtle/libfluidturtle.la \
	           -lm

#geotop_LIBTOOLFLAGS = -static
