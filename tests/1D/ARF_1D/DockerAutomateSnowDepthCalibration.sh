#!/bin/bash
for i in {1..1000}
do
	# define the prefix
	num=`printf %04d $i`

	# copy the appropriate geotop input file, which is created with R
	cp $(pwd)/SnowDepthCalibration/${num}_geotop.inpts $(pwd)/geotop.inpts

	# run geotop
	docker run --rm -v $(pwd):/work omslab/geotop
	
	# copy output file
	cp $(pwd)/output-tabs/point0001.txt $(pwd)/SnowDepthCalibration/${num}_point0001.txt
	
	# status update
	echo $num complete
done
