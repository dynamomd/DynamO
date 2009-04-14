#/bin/bash

NCOLL=20000000

#for FCC in $(seq 4 10) $(seq 12 2 20) $(seq 24 4 40) $(seq 45 5 70) $(seq 80 10 100);do
#    for DENS in 0.2 0.5 0.8 ; do
#	dynamod -m0 -C $FCC -d $DENS > /dev/null
#	
#	CellCount=$(dynarun config.out.xml.bz2 -c $NCOLL | grep "Cells <x,y,z>" | gawk '{print $4}' | gawk  -F',' '{print $1}')
#	
#	t1=$(bzcat output.xml.bz2 | xmlstarlet sel -t -v '/OutputData/Misc/Timing/CollPerSec/@val')
#	dynarun config.out.xml.bz2 -c $NCOLL > /dev/null
#	t2=$(bzcat output.xml.bz2 | xmlstarlet sel -t -v '/OutputData/Misc/Timing/CollPerSec/@val')
#	dynarun config.out.xml.bz2 -c $NCOLL > /dev/null
#	t3=$(bzcat output.xml.bz2 | xmlstarlet sel -t -v '/OutputData/Misc/Timing/CollPerSec/@val')
#
#	N=$(bzcat output.xml.bz2 | grep -i ParticleCount | gawk -F'"' '{print $2}')
#	celldens=$(gawk 'BEGIN {print '$N'/(('$CellCount')^3)}')
#
#	echo "$N $celldens $CellCount $t1 $t2 $t3" >> profile.Rho$DENS.results
#	echo "$N $celldens $CellCount $t1 $t2 $t3"
#    done
#done

for FCC in  30 ;do #20 25
    for DENS in $(seq 0.05 0.05 1.2) ; do
	dynamod -m0 -C $FCC -d $DENS > /dev/null
	
	CellCount=$(dynarun config.out.xml.bz2 -c $NCOLL | grep "Cells <x,y,z>" | gawk '{print $4}' | gawk  -F',' '{print $1}')
	
	t1=$(bzcat output.xml.bz2 | xmlstarlet sel -t -v '/OutputData/Misc/Timing/CollPerSec/@val')
	dynarun config.out.xml.bz2 -c $NCOLL > /dev/null
	t2=$(bzcat output.xml.bz2 | xmlstarlet sel -t -v '/OutputData/Misc/Timing/CollPerSec/@val')
	dynarun config.out.xml.bz2 -c $NCOLL > /dev/null
	t3=$(bzcat output.xml.bz2 | xmlstarlet sel -t -v '/OutputData/Misc/Timing/CollPerSec/@val')

	N=$(bzcat output.xml.bz2 | grep -i ParticleCount | gawk -F'"' '{print $2}')
	celldens=$(gawk 'BEGIN {print '$N'/(('$CellCount')^3)}')

	echo "$N $celldens $CellCount $t1 $t2 $t3" >> profile.F$FCC.results
	echo "$N $celldens $CellCount $t1 $t2 $t3"
    done
done
