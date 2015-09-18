#################################################################################
# Author: Matteo Dall'Amico
# Date:   13 oct 2013
# Comment: The binary will be created in the same folder where the makefile is.

HM          = .
PREV		= ..
LHM			= $(HM)
BINPATH 	= $(LHM)
NAME		= cloud_preproc
BINS		= $(BINPATH)/$(NAME)

SCRSPATH1	= $(HM)
SCRSPATH2	= $(PREV)/geotop
SCRSPATH3	= $(PREV)/libraries/fluidturtle
SCRSPATH4	= $(PREV)/libraries/ascii
SCRSPATH5	= $(PREV)/meteoio_plugin


OBJ		= 	$(SCRSPATH1)/cloud.o\
			$(SCRSPATH2)/radiation.o $(SCRSPATH2)/meteodata.o\
			$(SCRSPATH2)/clouds.o $(SCRSPATH2)/times.o\
			$(SCRSPATH2)/geotop_common.o $(SCRSPATH2)/parameters.o\
			$(SCRSPATH2)/meteo.o $(SCRSPATH2)/meteodistr.o\
			$(SCRSPATH2)/pedo.funct.o $(SCRSPATH2)/inputKeywords.o\
			$(SCRSPATH3)/tensors3D.o $(SCRSPATH3)/t_io.o\
			$(SCRSPATH3)/t_alloc.o\
			$(SCRSPATH4)/tabs.o $(SCRSPATH4)/rw_maps.o\
			$(SCRSPATH5)/meteoioplugin.o


BOOST_ROOT = /opt/local
METEOIO_ROOT = /usr/local
MAIN_DIR	= $(PREV)


CPPFLAGS = -g
INCLUDE = -I$(BOOST_ROOT)/include -I$(METEOIO_ROOT)/include -I$(SCRSPATH2) -I$(MAIN_DIR)

DEBUG   = -g -Wall
CPP	= g++

.cc.o: $*.cc $*.h
	$(CPP) $(CPPFLAGS) -c $< $(INCLUDE) -o $@

all: eurac

eurac: $(OBJ)
	$(CPP) -o $(BINS) $(OBJ) -L$(METEOIO_ROOT)/lib -L$(BOOST_ROOT)/lib -lm -lstdc++ -lmeteoio -lboost_program_options-mt -lboost_regex-mt -lboost_filesystem-mt -lboost_iostreams-mt -lboost_system-mt -lboost_unit_test_framework-mt -ldl

clean:
	rm -rf *.o *~ $(OBJ)
