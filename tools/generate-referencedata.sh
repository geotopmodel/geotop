#!/bin/bash
# This is part of the GEOtop 2.1 distribution 
# Bash Script to run GEOtop for  test cases with versions: 
# SE27XX , METEOIO-ON (geotop_dev), METEIO-OFF (geotop-dev)  
# and prepare reference data in the corrected directories 
# See test case: https://github.com/geotopmodel/geotop/tree/geotop_dev/tests/1D
#
# @author: Emanuele Cordano (emanuele.cordano@gmail.com) / Stefano Cozzini (stefano.cozzini@exact-lab.it)
# @date : 2016-11-26
#

show_help ()
{ # some help here..  
  echo "This script execute one or all test cases (1D or 3D) "
  echo " option: "
  echo " -h : this help"
  echo " -a [3D/1D]: run all 1D or 3D tests"
  echo " -f test_dir: run specific test contained in specified directory"
 } # 

#function3D 
run_3D () 
 { # this function execute the three runs on directory $1  
   
  echo $1
# Run all within the directory:
  cd $1 
for SUFFIX in  SE27XX METEOIO-ON METEOIO-OFF ; do
  echo "doing $SUFFIX for $1 with the following exe: $GEOTOP_BIN_$SUFFIX "
  $GEOTOP_BIN_$SUFFIX $1   1>$1/stdout.$SUFFIX 2>$1/stderr.$SUFFIX
  if [ $? -eq 0 ]; then
     echo "FINISHED OK moving on in changing data.."

# remove all previuos reference data and store the just produced ones..
      if [ -d "output-tabs"  ]; then 
         rm -rf $1/output-tabs-$SUFFIX
         mv $1/output-tabs  $1/output-tabs-$SUFFIX
         rm $1/output-tabs-$SUFFIX/placeholder
         git checkout output-tabs/placeholder
      fi 
      if [ -d "output-maps"  ]; then
         rm -rf $1/output-maps-$SUFFIX
         mv $1/output-maps  $1/output-maps-$SUFFIX
         rm $1/output-maps-$SUFFIX/placeholder
         git checkout $1/output-maps/placeholder
      fi   

# rename log file       
     if [ $SUFFIX == "SE27XX"  ]
     then
       mv $1/geotop.log  $1/geotop.log-$SUFFIX
     else 
       echo "renaming log files for $SUFFIX"
       mv $1/geotop2-1.log $1/geotop2-1.log-$SUFFIX  
     fi 
# if execution wrong 
  else 
     echo "EXECUTION FAILED: please check test "
  fi   
## all done     
     
done
 } 

# personal setting: please define where geotop distribution is 

GEOTOP_HOME=/home/ubuntu/geotopmodel/geotop

test_to_be_run="not_set"
#
export GEOTOP_BIN_=${GEOTOP_HOME}/bin/geotop-
export TEST_DIR=${GEOTOP_HOME}/tests

## Global Settings


# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
output_file=""
verbose=0

while getopts "h?f:a:" opt; do
    case "$opt" in
        h|\?)
	        show_help
	        exit 0
        ;;
        a)  all_tests=$OPTARG
        ;;
        f)  test_to_be_run=$OPTARG
        ;;
        esac
done
shift $((OPTIND-1))

 [ "$1" = "--" ] && shift

 
## Global Settings

export TEST_XD_DIR=$TEST_DIR/$all_tests 
echo $GEOTOP_BIN_SE27XX

echo $TEST_XD_DIR
cd  $TEST_XD_DIR


if [ $test_to_be_run != "not_set"  ]
then  
 {  echo "ready to run test in '$test_to_be_run'"
    if [ -d $test_to_be_run ]
    then
        echo item: $test_to_be_run
        echo GEOTOP_SIM_DIR: $TEST_XD_DIR/$test_to_be_run
	run_3D  $TEST_XD_DIR/$test_to_be_run
    else
        echo "'$test_to_be_run' is not a test directory"
    fi 
    exit
  }
fi

for i in $( ls ); do
    if [ -d $TEST_XD_DIR/$i ]
    	then 
          echo item: $i
          echo GEOTOP_SIM_DIR: $TEST_XD_DIR/$i
	  run_3D  $TEST_XD_DIR/$i
    else 
           echo "$i is not a test directory"
	   continue 
    fi
done 

exit
