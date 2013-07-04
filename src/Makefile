EXECUTABLE = geotop/geotop

$(EXECUTABLE):
	(cd libraries; $(MAKE) $(MFLAGS))
	(cd geotop; $(MAKE) $(MFLAGS))

all: $(EXECUTABLE)

clean:
	(cd geotop; make clean)
	(cd libraries; make clean)
