#!/bin/bash

for file in $(ls -1d *.data.xml.bz2); 
  do
  E2=$(bzcat $file | xml sel -t -v '/OutputData/Energy/InternalEnergy/@SquareAvg')
  E=$(bzcat $file | xml sel -t -v '/OutputData/Energy/InternalEnergy/@Avg')
  T=$(bzcat $file | xml sel -t -v '/OutputData/Energy/T/@val')
  echo $E $E2 $T | gawk '{print $3,($2 - $1*$1) / ($3*$3)}'
done