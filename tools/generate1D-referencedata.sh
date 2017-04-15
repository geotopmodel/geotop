#!/bin/bash
# This is part of the GEOtop 2.1 distribution 
# Bash Script to run GEOtop for 1D test case with versions: 
# SE27XX , METEOIO-ON (geotop_dev), METEIO-OFF (geotop-dev)  
# and prepare reference data in the corrected directories 
# See test case: https://github.com/geotopmodel/geotop/tree/geotop_dev/tests/1D
#
# @author: Emanuele Cordano (emanuele.cordano@gmail.com) / Stefano Cozzini (stefano.cozzini@exact-lab.it)
# @date : 2016-11-26
#
# Personal Settings

GEOTOP_HOME=/Users/cozzini/GEOTOP/github/geotop

#
export GEOTOP_BIN_SE27XX=${GEOTOP_HOME}/bin/geotop-2.0.0 
export GEOTOP_BIN_METEOIO_ON=${GEOTOP_HOME}/bin/geotop-2.1.0-METEOIO-ON
export GEOTOP_BIN_METEOIO_OFF=${GEOTOP_HOME}/bin/geotop-2.1.0-METEOIO-OFF
export TEST_DIR=${GEOTOP_HOME}/tests

## Global Settings

export TEST_1D_DIR=$TEST_DIR/1D 

echo $GEOTOP_BIN_SE27XX
echo $TEST_1D_DIR

cd $TEST_1D_DIR
for i in $( ls ); do
    if [ -d $i ]
    	then 
          echo item: $i
          echo GEOTOP_SIM_DIR: $TEST_1D_DIR/$i
    else 
           echo "$i is not a test directory"
	   continue 
    fi
 
# Run all within the directory:  
    cd $i
# first GEOTOP_SE27XX 
      echo "doing SE27XX for $i"
      $GEOTOP_BIN_SE27XX $TEST_1D_DIR/$i    1>$TEST_1D_DIR/$i/stdout.SE27XX 2>$TEST_1D_DIR/$i/stderr.SE27XX
      echo FINISHED: $GEOTOP_BIN_SE27XX $TEST_1D_DIR/$i
      rm -rf $TEST_1D_DIR/$i/output-tabs-SE27XX
      mv $TEST_1D_DIR/$i/output-tabs  $TEST_1D_DIR/$i/output-tabs-SE27XX
      rm  $TEST_1D_DIR/$i/output-tabs-SE27XX/placeholder
      mv $TEST_1D_DIR/$i/geotop.log  $TEST_1D_DIR/$i/geotop.log-SE27XX
      git checkout $TEST_1D_DIR/$i/output-tabs/placeholder
	

## Run GEOTOP_METEIO-ON 

      echo "doing METEO-IO on for $i"
      $GEOTOP_BIN_METEOIO_ON $TEST_1D_DIR/$i  1>$TEST_1D_DIR/$i/stdout.METEOIO_ON 2>$TEST_1D_DIR/$i/stderr.METEOIO_ON
      echo FINISHED: $GEOTOP_METEOIO_ON $TEST_1D_DIR/$i
      rm -rf $TEST_1D_DIR/$i/output-tabs-METEOIO-ON
      mv $TEST_1D_DIR/$i/output-tabs  $TEST_1D_DIR/$i/output-tabs-METEOIO-ON
      rm  $TEST_1D_DIR/$i/output-tabs-METEOIO-ON/placeholder
      mv $TEST_1D_DIR/$i/geotop2-1.log  $TEST_1D_DIR/$i/geotop2-1.log-METEOIO-ON
      git checkout $TEST_1D_DIR/$i/output-tabs/placeholder

 
## Run GEOTOP_METEIO-OFF

      echo "doing METEO-IO off  for $i"
      $GEOTOP_BIN_METEOIO_OFF $TEST_1D_DIR/$i  1>$TEST_1D_DIR/$i/stdout.METEOIO_OFF 2>$TEST_1D_DIR/$i/stderr.METEOIO_OFF
      echo FINISHED: $GEOTOP_BIN_METEOIO_OFF $TEST_1D_DIR/$i
      rm -rf $TEST_1D_DIR/$i/output-tabs-METEOIO-OFF
      mv $TEST_1D_DIR/$i/output-tabs  $TEST_1D_DIR/$i/output-tabs-METEOIO-OFF
      rm  $TEST_1D_DIR/$i/output-tabs-METEOIO-OFF/placeholder
      mv $TEST_1D_DIR/$i/geotop2-1.log  $TEST_1D_DIR/$i/geotop2-1.log-METEOIO-OFF
      git checkout $TEST_1D_DIR/$i/output-tabs/placeholder
 
  
done






