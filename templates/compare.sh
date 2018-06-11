#!/bin/bash

declare -i st
st=0
for _f in $(ls -1 output-tabs/*.txt output-maps/*.txt)
do
    _orig=${_f/\//-SE27XX\/}
    numdiff -a 1e-5 -r 1e-5 -s ' \t\n=,:;<>[](){}^' $_f $_orig &> _out.$$
    if (( $? != 0 ))
    then
	echo "Comparing $_f and $_orig"
	cat _out.$$
	st=1
    fi
	
done > failing_output

if (( $st == 1 ))
then
    echo ""
    echo ""
    echo "########################################################################"
    echo "Test in folder $PWD failed."
    echo "Here are the first 10 lines."
    head -n 10 failing_output
    echo ''
    echo ''
    echo "Full failing diff can be found here"
    echo ''
    echo "         $PWD/failing_output"
    echo ''
    echo "########################################################################"
    echo ""
    echo ""
fi

rm _out.$$

exit $st


