#!/bin/bash
cat $1 | grep vertex | gawk '{print $2+0,$3+0,$4+0}' > tmp.vertex

xmin=$(cat tmp.vertex | gawk 'BEGIN {xval=1e200} {if ($1 < xval) xval = $1 } END {print xval}')
ymin=$(cat tmp.vertex | gawk 'BEGIN {xval=1e200} {if ($2 < xval) xval = $2 } END {print xval}')
zmin=$(cat tmp.vertex | gawk 'BEGIN {xval=1e200} {if ($3 < xval) xval = $3 } END {print xval}')

echo Min vals are $xmin $ymin $zmin

cat tmp.vertex \
    | gawk '{print $1-'$xmin',$2-'$ymin',$3-'$zmin'}' \
    > tmp2.vertex

xmax=$(cat tmp2.vertex | gawk 'func fabs(val){return (val<0)?-val:val;} BEGIN {xval=0} {if (fabs($1)>xval) xval = fabs($1)} END {print xval}')
ymax=$(cat tmp2.vertex | gawk 'func fabs(val){return (val<0)?-val:val;} BEGIN {xval=0} {if (fabs($2)>xval) xval = fabs($2)} END {print xval}')
zmax=$(cat tmp2.vertex | gawk 'func fabs(val){return (val<0)?-val:val;} BEGIN {xval=0} {if (fabs($3)>xval) xval = fabs($3)} END {print xval}')

echo Max lengths are $xmax $ymax $zmax

max=$xmax

maxord=$(echo $xmax $ymax $zmax | gawk '{val = $1; if ($2 > val) val = $2; if ($3 > val) val = $3} END {print val}')

echo Maximum ordinate val is $maxord

cat tmp2.vertex \
    | gawk '{print $1/'$maxord' - 0.5, $2/'$maxord' - 0.5, $3/'$maxord' - 0.5}' \
    > fixed.vertex

#rm tmp.vertex tmp2.vertex
