#!/bin/bash
cat $1 | grep vertex | gawk '{print $2,$3,$4}' > tmp.vertex

xavg=$(cat tmp.vertex | gawk 'BEGIN {xval=0} {xval+=$1} END {print xval/NR}')
yavg=$(cat tmp.vertex | gawk 'BEGIN {xval=0} {xval+=$2} END {print xval/NR}')
zavg=$(cat tmp.vertex | gawk 'BEGIN {xval=0} {xval+=$3} END {print xval/NR}')

echo Averages are $xavg $yavg $zavg

cat tmp.vertex \
    | gawk '{print $1-'$xavg',$2-'$yavg',$3-'$zavg'}' \
    > tmp2.vertex

xmax=$(cat tmp2.vertex | gawk 'func fabs(val){return (val<0)?-val:val;} BEGIN {xval=0} {if (fabs($1)>xval) xval = fabs($1)} END {print xval}')
ymax=$(cat tmp2.vertex | gawk 'func fabs(val){return (val<0)?-val:val;} BEGIN {xval=0} {if (fabs($2)>xval) xval = fabs($2)} END {print xval}')
zmax=$(cat tmp2.vertex | gawk 'func fabs(val){return (val<0)?-val:val;} BEGIN {xval=0} {if (fabs($3)>xval) xval = fabs($3)} END {print xval}')

echo Maximums are $xmax $ymax $zmax

max=$xmax

maxord=$(echo $xmax $ymax $zmax | gawk '{val = $1; if ($2 > val) val = $2; if ($3 > val) val = $3} END {print val}')

echo Maximum ordinate val is $maxord

cat tmp2.vertex \
    | gawk '{print $1/'$maxord',$2/'$maxord',$3/'$maxord'}' \
    > fixed.vertex

rm tmp.vertex tmp2.vertex
