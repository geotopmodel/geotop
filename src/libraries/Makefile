CC = gcc

OBJLIBS = libascii.a libfluidturtle.a libgeomorphology.a libmath.a
DIRS = ascii fluidturtle geomorphology math

all: $(OBJLIBS)

libascii.a:
	cd ascii; $(MAKE) $(MFLAGS)

libfluidturtle.a:
	cd fluidturtle; $(MAKE) $(MFLAGS)

libgeomorphology.a:
	cd geomorphology; $(MAKE) $(MFLAGS)

libmath.a:
	cd math; $(MAKE) $(MFLAGS)

clean:
	$(RM) -f $(OBJLIBS)
	for d in $(DIRS); do (cd $$d; $(MAKE) clean); done	
