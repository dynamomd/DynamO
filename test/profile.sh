#>profile.results

for FCC in $(seq 51 60);do
    for DENS in 0.5; do
	./configmod -i2 -C $FCC -d $DENS > /dev/null
	
	CellCount=$(./mdrun config.out.xml.bz2 -c 1000000 | grep Cells | gawk '{print $4}' | gawk  -F',' '{print $1}')
	
	t=$(bzcat output.xml.bz2 | xml sel -t -v '/OutputData/Misc/Timing/CollPerSec/@val')
	
	N=$(bzcat output.xml.bz2 | grep -i ParticleCount | gawk -F'"' '{print $2}')
	celldens=$(gawk 'BEGIN {print '$N'/(('$CellCount')^3)}')

	echo "$N $celldens $t $CellCount" >> profile.results
	echo "$N $celldens $t $CellCount"
    done
done
