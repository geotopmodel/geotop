
HM          = src
BINPATH 	= bin
NAME		= geotop-2.0.0
BINS		= $(BINPATH)/$(NAME)

SCRSPATH1	= $(HM)/geotop
LIBPATH1	= $(HM)/libraries/fluidturtle
LIBPATH2	= $(HM)/libraries/ascii
LIBPATH3	= $(HM)/libraries/geomorphology
LIBPATH4	= $(HM)/libraries/math

SRC	= $(SCRSPATH1)/blowingsnow.c $(SCRSPATH1)/clouds.c $(SCRSPATH1)/energy.balance.c $(SCRSPATH1)/geotop.c \
	$(SCRSPATH1)/input.c $(SCRSPATH1)/meteo.c $(SCRSPATH1)/meteodistr.c $(SCRSPATH1)/output.c \
	$(SCRSPATH1)/parameters.c $(SCRSPATH1)/PBSM.c $(SCRSPATH1)/pedo.funct.c $(SCRSPATH1)/radiation.c \
	$(SCRSPATH1)/snow.c $(SCRSPATH1)/tables.c $(SCRSPATH1)/times.c $(SCRSPATH1)/turbulence.c \
	$(SCRSPATH1)/vegetation.c $(SCRSPATH1)/water.balance.c \
	$(SCRSPATH1)/channels.c $(SCRSPATH1)/deallocate.c \
	$(SCRSPATH1)/indices.c $(SCRSPATH1)/meteodata.c  $(SCRSPATH1)/recovering.c \
	$(LIBPATH1)/alloc.c $(LIBPATH1)/datamanipulation.c $(LIBPATH1)/error.c  $(LIBPATH1)/linearalgebra.c \
	$(LIBPATH1)/list.c $(LIBPATH1)/probability.c  $(LIBPATH1)/random.c $(LIBPATH1)/statistics.c \
    $(LIBPATH1)/string.c  $(LIBPATH1)/t_io.c   $(LIBPATH1)/tensors3D.c  $(LIBPATH1)/utilities.c  $(LIBPATH1)/write_dem.c  \
    $(LIBPATH2)/import_ascii.c  $(LIBPATH2)/rw_maps.c  $(LIBPATH2)/write_ascii.c  $(LIBPATH2)/tabs.c $(LIBPATH2)/init.c\
	$(LIBPATH3)/networks.c  $(LIBPATH3)/geomorphology.0875.c $(LIBPATH3)/geomorphology.c  $(LIBPATH3)/shadows.c   $(LIBPATH3)/dtm_resolution.c  \
	$(LIBPATH4)/util_math.c $(LIBPATH4)/sparse_matrix.c 
            
HPATH1 	= $(LIBPATH1)
HPATH2	= $(LIBPATH2)
HPATH3  = $(LIBPATH3)
HPATH4  = $(LIBPATH4)
HPATH0  = $(SCRSPATH1)

FFLAGS	= -O3 -g
CFLAGS	= -O3 -g -I$(HPATH1) -I$(HPATH2) -I$(HPATH3) -I$(HPATH4) -I$(HPATH0) 

DEBUG = -g
CC	= gcc $(DEBUG)

geotop2_0_0:$(SRC)
			$(CC) $(CFLAGS) -o $(BINS) $(SRC) -lm

