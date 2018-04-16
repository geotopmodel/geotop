#!/bin/bash

_cpp=geotop2-1.log
_cc=geotop.log-SE27XX

awk '/TimeStepEnergyAndWater/,/SubfolderRecoveryFiles/' ${_cpp} | sed -e 's/\[NOTICE\]: //;s/\[.*\]//;s/(default)//' > _cleaned_cpp 
awk '/TimeStepEnergyAndWater/,/SubfolderRecoveryFiles/' ${_cc} | sed -e 's/\[.*\]//;s/(default)//' > _cleaned_cc

keys=$(cat _cleaned_cc | cut -d '=' -f 1)

# keys=(TimeStepEnergyAndWater SubfolderRecoveryFiles)
for _k in ${keys[*]}
do
    # echo $_k
    awk -v k=${_k} '$1 == k {print $0}' _cleaned_cpp > _k_cpp
    awk -v k=${_k} '$1 == k {print $0}' _cleaned_cc > _k_c
    numdiff _k_cpp _k_c &> _out
    if [ $? -ne 0 ]
    then
	echo ${_k}
	cat _out | sed '/+++/d;/\*\*\*/d'
    fi
done
# kk=TimeStepEnergyAndWater
# kk=TimeStepEnergyAndWater
# awk -v x=$kk '$0 ~ x {print $0}' _cleaned_cpp
# awk '/TimeStepEnergyAndWater/ {print $NF}' _cleaned_cc
