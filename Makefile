################################################################################
# Makefile for GEOtop program version 1.45
#################################################################################
# Author: Matteo Dall'Amico insipired by Thomas Egger
# Date:   3 feb 2012
# Comment: The binary will be created in the same folder where the makefile is.

HM          	= .
LHM		= $(HM)
BINPATH 	= $(LHM)
NAME		= GEOtop_4
BINS		= $(BINPATH)/$(NAME)

SCRSPATH1	= $(HM)/geotop
LIBPATH1	= $(HM)/libraries/fluidturtle
LIBPATH2	= $(HM)/libraries/ascii
LIBPATH3	= $(HM)/meteoio_plugin
LIBPATH4	= $(HM)/netCDF

OBJ		= 	$(SCRSPATH1)/blowingsnow.o $(SCRSPATH1)/channels.o\
		  	$(SCRSPATH1)/clouds.o $(SCRSPATH1)/deallocate.o\
		  	$(SCRSPATH1)/energy.balance.o $(SCRSPATH1)/geotop.o\
			$(SCRSPATH1)/water.balance.o $(SCRSPATH1)/output_nc.o\
		  	$(SCRSPATH1)/indices.o $(SCRSPATH1)/input.o\
			$(SCRSPATH1)/meteo.o $(SCRSPATH1)/meteodata.o\
			$(SCRSPATH1)/output.o $(SCRSPATH1)/parameters.o\
			$(SCRSPATH1)/PBSM.o $(SCRSPATH1)/pedo.funct.o\
			$(SCRSPATH1)/radiation.o $(SCRSPATH1)/recovering.o\
			$(SCRSPATH1)/snow.o $(SCRSPATH1)/tables.o\
			$(SCRSPATH1)/times.o $(SCRSPATH1)/turbulence.o\
			$(SCRSPATH1)/vegetation.o\
			$(LIBPATH1)/t_alloc.o $(LIBPATH1)/t_io.o\
			$(LIBPATH1)/tensors3D.o \
			$(LIBPATH2)/rw_maps.o $(LIBPATH2)/tabs.o\
			$(LIBPATH3)/meteoioplugin.o

NETCDFOBJ	= 	$(LIBPATH4)/netcdfIO.o $(LIBPATH4)/read_command_line.o\
			$(LIBPATH4)/read_command_line_netcdf.o

HPATH1 	= $(LIBPATH1)
HPATH2	= $(LIBPATH2)
HPATH3  = $(LIBPATH3)
HPATH0  = $(SCRSPATH1)
INCMETEOIO = /usr/local/include/meteoio

CPPFLAGS = -g 
INCLUDE = -I$(HPATH1) -I$(HPATH2) -I$(HPATH3) -I$(HPATH4) -I$(HPATH5) -I$(HPATH0) -I$(INCMETEOIO)

DEBUG   = -g -Wall
CPP	= g++

.cc.o: $*.cc $*.h
	$(CPP) $(CPPFLAGS) -c $< $(INCLUDE) -o $@

all: geotop

geotop: $(OBJ)
	$(CPP) -o $(BINS) $(OBJ) -rdynamic -lm -lstdc++ -lmeteoio -ldl 

netcdf: CPPFLAGS = -g -DUSE_NETCDF
netcdf: $(OBJ) $(NETCDFOBJ)
	$(CPP) -o $(BINS) $(OBJ) $(NETCDFOBJ) -rdynamic -lm -lstdc++ -lmeteoio -ldl -lnetcdf

.PHONY: tests

tests: geotop
	nosetests -v

clean:
	rm -rf *.o *~ $(OBJ) tests/test_sample_run/*/output_*/*
