#!/bin/bash

DIR=/home/elisa/Scrivania/MHPC/geotop_3.0/tests/1D/Matsch_B2_Ref_007

rm -rf ${DIR}/geotop.inpts ${DIR}/meteo0001.txt ${DIR}/*SE27XX*
rm ${DIR}/output-tabs/*
cp -r * ${DIR}
