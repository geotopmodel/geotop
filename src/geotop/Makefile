CC = gcc
CFLAGS = -I../libraries/fluidturtle -I../libraries/ascii -I../libraries/math -I../libraries/geomorphology
LDFLAGS = -lm

SRCS = blowingsnow.c channels.c clouds.c deallocate.c energy.balance.c \
	   geotop.c indices.c input.c meteo.c meteodata.c meteodistr.c \
	   output.c parameters.c PBSM.c pedo.funct.c radiation.c recovering.c \
	   snow.c tables.c times.c turbulence.c vegetation.c water.balance.c

OBJS	= ${SRCS:.c=.o}
LIBPATH = ../libraries
OBJLIBS = $(LIBPATH)/libascii.a $(LIBPATH)/libmath.a \
          $(LIBPATH)/libgeomorphology.a $(LIBPATH)/libfluidturtle.a

EXECUTABLE = geotop

all: $(EXECUTABLE)

.c.o:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

$(EXECUTABLE): $(OBJS) $(OBJLIBS)
	$(CC) $(OBJS)  $(LDFLAGS) $(OBJLIBS) -o $@ 

clean :
	-$(RM) -f $(OBJS) $(OBJLIB) $(EXECUTABLE)
