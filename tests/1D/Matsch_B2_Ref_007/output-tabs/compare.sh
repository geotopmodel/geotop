#!/bin/bash

_files=$(find . -name '*.txt')

for _f in ${_files[*]}:
do
    _o=${_f//\.\//../output-tabs-METEOIO-OFF\/}
    echo comparing ${_o} and ${_f}
    numdiff  -a 1e-5 -r 1e-8 -s ' \t\n=,:;<>[](){}^' ${_f} ${_o}
done	  
