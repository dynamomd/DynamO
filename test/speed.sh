#!/bin/bash

NUMRUN=5
NCOLL=10000000
> testvals
if [ ! -e config.speed.xml.bz2 ]; then
    echo "Creating Configuration....."
    is_configmod -m0 -C 30 -o config.test.xml.bz2 > /dev/null
    echo "Equilibrating configuration...."
    is_mdrun -c 1000000 config.test.xml.bz2 > /dev/null
    rm config.test.xml.bz2
    mv config.0.end.xml.bz2 config.speed.xml.bz2
fi 

for i in $(seq 0 $NUMRUN); do
    echo -n "Running test $i..."
    val=$(is_mdrun config.speed.xml.bz2 -c $NCOLL | grep "Avg Coll" | gawk '{print $4}')    
    echo $val >> testvals
    echo $val $(cat testvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}')
done

rm testvals config.0.end.xml.bz2 output.0.xml.bz2
