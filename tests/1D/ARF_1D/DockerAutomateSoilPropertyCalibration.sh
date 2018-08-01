#!/bin/bash

# declare fire array
declare -a fires=("Unburned" "Moderate" "Severe")

for i in {1..100}
	do for fire in "${fires[@]}"
		do
		# define the prefix
		num=`printf %04d $i`

		# copy the appropriate geotop input file, which is created with R
		cp $(pwd)/SoilPropertyCalibration/${num}_${fire}_SoilARF0001.txt $(pwd)/soil/SoilARF0001.txt

		# run geotop
		docker run --rm -v $(pwd):/work omslab/geotop

		# copy output file
		cp $(pwd)/output-tabs/point0001.txt $(pwd)/SoilPropertyCalibration/output_${num}_${fire}_point0001.txt
		cp $(pwd)/output-tabs/soiltemp0001.txt $(pwd)/SoilPropertyCalibration/output_${num}_${fire}_soiltemp0001.txt
		cp $(pwd)/output-tabs/thetaliq0001.txt $(pwd)/SoilPropertyCalibration/output_${num}_${fire}_thetaliq0001.txt
		cp $(pwd)/output-tabs/thetaice0001.txt $(pwd)/SoilPropertyCalibration/output_${num}_${fire}_thetaice0001.txt

		# status update
		echo "$num" "$fire" complete
		done
	done
