################################################################################
# Makefile for GEOtop program version 1.45
#################################################################################
# Author: Matteo Dall'Amico insipired by Thomas Egger
# Date:   29 aug 2011
# Comment: The binary will be created in the same folder where the makefile is.

HM          	= .
LHM		= $(HM)
BINPATH 	= $(LHM)/
NAME		= GEOtop1.45
BINS		= $(BINPATH)/$(NAME)

SCRSPATH1	= $(HM)/src/geotop
LIBPATH1	= $(HM)/src/libraries/fluidturtle
LIBPATH2	= $(HM)/src/libraries/ascii
LIBPATH3	= $(HM)/src/libraries/geomorphology
LIBPATH4	= $(HM)/src/libraries/math
LIBPATH5	= $(HM)/src/gt_utilities
LIBPATH6	= $(HM)/src/MeteoIO_plug

OBJ		= 	$(SCRSPATH1)/blowingsnow.o $(SCRSPATH1)/channels.o\
		  	$(SCRSPATH1)/clouds.o $(SCRSPATH1)/deallocate.o\
		  	$(SCRSPATH1)/energy.balance.o $(SCRSPATH1)/geotop.o\
		  	$(SCRSPATH1)/indices.o $(SCRSPATH1)/input.o\
			$(SCRSPATH1)/meteo.o\
			$(SCRSPATH1)/meteodata.o $(SCRSPATH1)/meteodistr.o\
			$(SCRSPATH1)/output.o $(SCRSPATH1)/parameters.o\
			$(SCRSPATH1)/PBSM.o $(SCRSPATH1)/pedo.funct.o\
			$(SCRSPATH1)/radiation.o $(SCRSPATH1)/recovering.o\
			$(SCRSPATH1)/snow.o $(SCRSPATH1)/tables.o\
			$(SCRSPATH1)/times.o $(SCRSPATH1)/turbulence.o\
			$(SCRSPATH1)/vegetation.o $(SCRSPATH1)/water.balance.o\
			$(LIBPATH1)/alloc.o $(LIBPATH1)/error.o\
			$(LIBPATH1)/list.o $(LIBPATH1)/t_io.o\
			$(LIBPATH1)/tensors3D.o  $(LIBPATH1)/utilities.o\
			$(LIBPATH1)/datamanipulation.o $(LIBPATH1)/random.o\
			$(LIBPATH1)/linearalgebra.o $(LIBPATH1)/write_dem.o\
			$(LIBPATH1)/probability.o $(LIBPATH1)/statistics.o\
			$(LIBPATH1)/string.o\
			$(LIBPATH2)/import_ascii.o $(LIBPATH2)/rw_maps.o\
			$(LIBPATH2)/write_ascii.o $(LIBPATH2)/tabs.o\
			$(LIBPATH2)/init.o\
			$(LIBPATH3)/networks.o $(LIBPATH3)/geomorphology.0875.o\
			$(LIBPATH3)/shadows.o  $(LIBPATH3)/dtm_resolution.o\
			$(LIBPATH3)/geomorphology.o\
			$(LIBPATH4)/util_math.o $(LIBPATH4)/sparse_matrix.o\
			$(LIBPATH5)/read_command_line.o\
			$(LIBPATH6)/meteoioplugin.o

HPATH1 	= $(LIBPATH1)
HPATH2	= $(LIBPATH2)
HPATH3  = $(LIBPATH3)
HPATH4  = $(LIBPATH4)
HPATH5  = $(LIBPATH5)
HPATH6  = $(LIBPATH6)
HPATH0  = $(SCRSPATH1)

CFLAGS	= -O3 -g 
#INCLUDE = -I$(HPATH1) -I$(HPATH2) -I$(HPATH3) -I$(HPATH4) -I$(HPATH0)

DEBUG   = -g -Wall
CC	= gcc $(DEBUG)
CPP	= g++

.cc.o: $*.cc $*.h
	$(CPP) $(CPPFLAGS) -c $< $(INCLUDE) -o $@

.c.o: $*.c $*.h
	$(CC) $(CFLAGS) -c $< $(INCLUDE) -o $@


all: geotop

geotop: $(OBJ)
	$(CC) -o $(BINS) $(OBJ) -lm -rdynamic -ldl -lstdc++ -lmeteoio
clean:
	rm -rf *.o *~ $(OBJ)
