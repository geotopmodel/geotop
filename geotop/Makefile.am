bin_PROGRAMS = geotop

geotop_SOURCES = blowingsnow.cc geotop_common.cc keywords.cc channels.cc clouds.cc deallocate.cc energy.balance.cc geotop.cc indices.cc input.cc meteo.cc meteodata.cc meteodistr.cc output.cc output_nc.cc parameters.cc PBSM.cc pedo.funct.cc pure_gt.cc radiation.cc recovering.cc snow.cc tables.cc times.cc turbulence.cc vegetation.cc water.balance.cc 

AM_CPPFLAGS = -I$(srcdir) \
              -I$(top_srcdir)/libraries/ascii \
              -I$(top_srcdir)/libraries/fluidturtle
AM_CFLAGS =

AM_LDFLAGS = -static

if WITH_STATIC_GEOTOP
    AM_LDFLAGS += -all-static
endif

geotop_LDADD = $(top_srcdir)/libraries/ascii/libascii.la \
               $(top_srcdir)/libraries/fluidturtle/libfluidturtle.la \
               # $(top_srcdir)/wrf/libwrfio.la \
               $(top_srcdir)/meteoio_plugin/libmeteoioplugin.la \
	           -lm $(METEOIO_LIBS)

if INCLUDE_NETCDF_SUBDIR
  geotop_LDADD += $(top_srcdir)/netCDF/libnetcdfio.la
endif

#geotop_LIBTOOLFLAGS = -static
