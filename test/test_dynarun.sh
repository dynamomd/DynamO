#!/bin/bash
#    DYNAMO:- Event driven molecular dynamics simulator 
#    http://www.marcusbannerman.co.uk/dynamo
#    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    version 3 as published by the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Dynarun="../bin/dynarun"
Dynamod="../bin/dynamod"

#Next is the name of XML starlet
Xml="xml"

if [ ! -x $Dynarun ]; then 
    echo "Could not find dynarun, have you built it?"
fi

if [ ! -x $Dynamod ]; then 
    echo "Could not find dynamod, have you built it?"
fi

which $Xml || Xml="xmlstarlet"

which $Xml || `echo "Could not find XMLStarlet"; exit`

which gawk || `echo "Could not find gawk"; exit`

#We create a local copy of the executables, so that recompilation won't break running tests
cp $Dynamod ./dynamod
cp $Dynarun ./dynarun

function HS_replex_test {
    for i in $(seq 0 2); do
	./dynamod -m 0 -C 7 -T $(echo "0.5*$i + 0.5" | bc -l) \
	    -o config.$i.start1.xml.bz2 > run.log

	bzcat config.$i.start1.xml.bz2 | \
	    $Xml ed -u '//Simulation/Scheduler/@Type' -v "$1" \
	    | bzip2 > config.$i.start.xml.bz2
	
	rm config.$i.start1.xml.bz2
    done

    #Equilibration
    time ./dynarun --engine 2 $2 -i 10 -f 100 \
	config.*.start.xml.bz2 >> run.log
    
    #Production
    time ./dynarun --engine 2 $2 -i 10 -f 200 config.*.end.xml.bz2 >> run.log

    if [ ! -e "output.0.xml.bz2" ]; then
	echo "$1 HS1 Replica Exchange -: FAILED Could not find output.0.xml.bz2"
	exit 1
    fi

    if [ ! -e "output.1.xml.bz2" ]; then
	echo "$1 HS1 Replica Exchange -: FAILED Could not find output.1.xml.bz2"
	exit 1
    fi

    if [ ! -e "output.2.xml.bz2" ]; then
	echo "$1 HS1 Replica Exchange -: FAILED Could not find output.2.xml.bz2"
	exit 1
    fi

    MFT1=$(bzcat output.0.xml.bz2 | $Xml sel -t -v "/OutputData/Misc/totMeanFreeTime/@val")

    MFT2=$(bzcat output.1.xml.bz2 | $Xml sel -t -v "/OutputData/Misc/totMeanFreeTime/@val")

    MFT3=$(bzcat output.2.xml.bz2 | $Xml sel -t -v "/OutputData/Misc/totMeanFreeTime/@val")

    MFT1=$(echo "$MFT1 * sqrt(0.5)" | bc -l)
    MFT3=$(echo "$MFT3 * sqrt(1.5)" | bc -l)

    AVG=$(echo "($MFT1 + $MFT2 + $MFT3) / 3.0" | bc -l)

    pass=$(echo $MFT1 $AVG | gawk '{a=int(100 * (($1 / $2)-1.0)); if (100 * (($1 / $2)-1.0)<0) a=-a; print a}' 2>&1)
    if [ $pass != 0 ]; then 
	echo "$1 HS1 Replica Exchange -: FAILED $MFT1 $MFT2 $MFT3 $AVG"
	exit 1
    else
	echo -n $(echo $MFT1 $AVG | gawk '{print $1 / $2}')" "
    fi

    pass=$(echo $MFT2 $AVG | gawk '{a=int(100 * (($1 / $2)-1.0)); if (100 * (($1 / $2)-1.0)<0) a=-a; print a}')
    if [ $pass != 0 ]; then 
	echo "$1 HS2 Replica Exchange -: FAILED $MFT1 $MFT2 $MFT3 $AVG"
	exit 1 	
    else
	echo -n $(echo $MFT2 $AVG | gawk '{print $1 / $AVG}')" "
    fi

    pass=$(echo $MFT3 $AVG | gawk '{a=int(100 * (($1 / $2)-1.0)); if (100 * (($1 / $2)-1.0)<0) a=-a; print a}')
    if [ $pass != 0 ]; then 
	echo "$1 HS3 Replica Exchange -: FAILED $MFT1 $MFT2 $MFT3 $AVG"
	exit 1 	
    else
	echo -n $(echo $MFT3 $AVG | gawk '{print $1 / $2}')" "
    fi

    echo "$1 HS Replica Exchange -: PASSED"
    rm config.*.end.xml.bz2 config.*.start.xml.bz2 \
	output.*.xml.bz2 *.dat replex.stats
}

function wallsw_bailout { 
    echo "$1 WallSW -: FAILED"
    exit 1 
}


function Ring_compressiontest { 
    ./dynamod -m 7 --f3 0 &> run.log
    ./dynamod -m 3 --s1 config.out.xml.bz2 --i1 1 -C 4 -d 0.01 -o tmp.xml.bz2 &> run.log
    bzcat tmp.xml.bz2 | xmlstarlet ed -u '//Interaction[@Type="SquareBond"]/@End' -v 2559 -u '//BC/@Type' -v "PBC" | bzip2 > config.out.xml.bz2
    ./dynarun --engine 3 --target-pack-frac 0.5 config.out.xml.bz2 &> run2.log

    #Check for errors
    > error.log
    ./dynamod --check config.out.xml.bz2 > run.log 2> error.log

    echo -n "Ring polymer compression -: "
    if [ -s error.log ]; then
	echo "FAILED: Check error.log for the overlaps!"
	exit 1;
    else
	echo "PASSED"
    fi

    #Cleanup 
    rm -Rf config.out.xml.bz2 output.xml.bz2 tmp.xml.bz2 run.log error.log
}

function SWpolymer_compressiontest {
    Nc=2
    C=8
    Mode=2
    T=1.3
    compress_density=0.1
    density=0.3
    N=$((Nc*((4/(2**Mode))*(C**3))))


    RANDOM_SEED="-s 1234"
    ./dynamod $RANDOM_SEED -m 7 --i1 $((Nc/2)) --f3 1.5 --b1 -o config.polymer.xml.bz2 &> run.log

    if [ $(echo "$density <= $compress_density" |bc) -eq 1 ]; then
    #If the density is lower than the compress_density, just directly make the big system
	./dynamod $RANDOM_SEED -m 3 -T $T -r $T --i1 $Mode --b1 --s1 config.polymer.xml.bz2 -d $density -C $C -o config.tmp2.xml.bz2 &>> run.log
	bzcat config.tmp2.xml.bz2 | xmlstarlet ed -u '//Interaction[1]/IDPairRange/@End' -v $((N-1)) -u '//BC/@Type' -v "PBC" | bzip2 > config.tmp.xml.bz2
    else
    #We have to compress as the density is too high
	./dynamod $RANDOM_SEED -m 3 --i1 $Mode --b1 --s1 config.polymer.xml.bz2 -d $compress_density -C $C -o config.tmp2.xml.bz2 &>> run.log
	bzcat config.tmp2.xml.bz2 | xmlstarlet ed -u '//Interaction[1]/IDPairRange/@End' -v $((N-1)) -u '//BC/@Type' -v "PBC" | bzip2 > config.packed.xml.bz2
	./dynarun $RANDOM_SEED config.packed.xml.bz2 --engine 3 --target-density $density -o config.compressed.xml.bz2 &>> run.log
	./dynamod $RANDOM_SEED config.compressed.xml.bz2 -r $T -T $T -Z -o config.rescaled.xml.bz2 &>> run.log
	cp config.rescaled.xml.bz2 config.tmp.xml.bz2
    fi

    #Finally, fix the definitions of the structures in the config file
    bzcat config.tmp.xml.bz2 | xmlstarlet ed -d '//Molecule' > config.tmp.xml
    rm config.tmp.xml.bz2
    for start in $(seq 0 $Nc $((N-1))); do
	cat config.tmp.xml | xmlstarlet ed \
            -s '//Structure' -t elem -n "Molecule" -v " " \
            -s '//Molecule[last()]' -t elem -n "IDRange" -v " " \
            -s '//Molecule[last()]/IDRange' -t attr -n "Type" -v "Ranged" \
            -s '//Molecule[last()]/IDRange' -t attr -n "Start" -v $start \
            -s '//Molecule[last()]/IDRange' -t attr -n "End" -v $((start+Nc-1)) \
            > config.tmp2.xml
	mv config.tmp2.xml config.tmp.xml
    done
    bzip2 config.tmp.xml
    mv config.tmp.xml.bz2 config.start.xml.bz2


    ./dynarun $RANDOM_SEED config.start.xml.bz2 -c 300000 &>> run.log
    MFT="0.0220444033732851"

    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{var=($1-'$MFT')/'$MFT'; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "SWpolymer -: FAILED, Measured MFT =" $(bzcat output.xml.bz2 \
		| $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val') \
		", expected MFT =" $MFT
	    exit 1
	else
	    echo "SWpolymer -: PASSED, Measured MFT =" $(bzcat output.xml.bz2 \
		| $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val') \
		", expected MFT =" $MFT
	fi
    else
	echo "Error, no output.0.xml.bz2 in SWpolymer test"
	exit 1
    fi
}

function wallsw {
#Just tests the square shoulder interaction between two walls
    cat wallsw.xml | \
	$Xml ed -u '//Simulation/Scheduler/@Type' -v "$1" \
	> tmp.xml
    #Any multiple of 5 will/should always give the original configuration
    #doing 9995 as this stops any 2 periodicity
    ./dynarun -c 9995 $2 tmp.xml &> run.log
    
    ./dynamod --round config.out.xml.bz2 > /dev/null
    ./dynamod config.out.xml.bz2 > /dev/null

    bzcat config.out.xml.bz2 | \
	$Xml sel -t -c '//ParticleData' > testresult.dat

    cat wallsw.xml | \
	$Xml sel -t -c '//ParticleData' > correct.dat
    
    diff testresult.dat correct.dat &> /dev/null \
	|| wallsw_bailout $1
    
    echo "$1 WallSW -: PASSED"

    #Cleanup
    rm -Rf config.out.xml.bz2 output.xml.bz2 tmp.xml.bz2 run.log \
	testresult.dat correct.dat
}

function umbrella_bailout { 
    echo "$1 Umbrella -: FAILED"
    exit 1 
}

function umbrella {
#Just tests the square shoulder interaction between two walls
    bzcat umbrella.xml | \
	$Xml ed -u '//Simulation/Scheduler/@Type' -v "$1" \
	| bzip2 > tmp.xml.bz2
    
    #Any multiple of 12 will/should always give the original configuration
    #doing only 12 to stop error creeping in
    ./dynarun -c 12 $2 tmp.xml.bz2 &> run.log

    ##This rounds the last digit off 
    ./dynamod --round config.out.xml.bz2 > /dev/null
    ./dynamod config.out.xml.bz2 > /dev/null

    bzcat config.out.xml.bz2 | \
	$Xml sel -t -c '//ParticleData' > testresult.dat

    bzcat umbrella.xml | \
	$Xml sel -t -c '//ParticleData' > correct.dat
    
    diff testresult.dat correct.dat &> /dev/null \
	|| umbrella_bailout $1
    
    echo "$1 Umbrella -: PASSED"

    #Cleanup
    rm -Rf config.out.xml.bz2 output.xml.bz2 tmp.xml.bz2 run.log \
	testresult.dat correct.dat
}

function BinarySphereTest {
    > run.log

    ./dynamod -s1 -m 8 --f3 0.05 -d 1.4 -C 10 --f1 0.5 &> run.log
    bzcat config.out.xml.bz2 \
	| $Xml ed -u "//Globals/Global[@Name='SchedulerNBList']/@Type" \
	-v "$1" | bzip2 > tmp.xml.bz2

    ./dynarun -c 1000000 tmp.xml.bz2 >> run.log 2>&1
    ./dynarun -c 1000000 config.out.xml.bz2 >> run.log 2>&1
    
    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{var=($1-0.0098213311089127)/0.0098213311089127; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "BinarySphereTest -: FAILED"
	    exit 1
	else
	    echo "BinarySphereTest -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in Binary Sphere Test"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 config.out.xml.bz2 output.xml.bz2 \
	tmp.xml.bz2 run.log
}

function ShearingTest {
    > run.log

    ./dynamod -s1 -m 4 --f1 0.9 &> run.log    
    ./dynarun -c 500000 config.out.xml.bz2 >> run.log 2>&1
    ./dynarun -c 1000000 config.out.xml.bz2 >> run.log 2>&1
    
    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{var=($1-0.113195634)/0.113195634; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "ShearingTest -: FAILED"
	    exit 1
	else
	    echo "ShearingTest -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in Shearing test"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 config.out.xml.bz2 output.xml.bz2 \
	tmp.xml.bz2 run.log
}

function IsolatedPolymerTest {
    > run.log

    ./dynamod -s1 -m 2 --i1 50  &> run.log
    ./dynarun -s2 -c 1000000 config.out.xml.bz2 >> run.log 2>&1
    ./dynarun -s3 -c 1000000 config.out.xml.bz2 >> run.log 2>&1
    
    MFT=0.0240657157464771
    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{var=($1-'$MFT')/'$MFT'; print ((var < 0.08) && (var > -0.08))}') != "1" ]; then
	    echo "IsolatedPolymerTest -: FAILED MFT_expected=" $MFT " MFT_measured=" $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val')
	    exit 1
	else
	    echo "IsolatedPolymerTest -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in IsolatedPolymer test"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 config.out.xml.bz2 output.xml.bz2 \
	tmp.xml.bz2 run.log
}

function HardLinesTest {
    > run.log

    dens=0.1
    ./dynamod -s1 -m 9 -C 1000 -d $dens  &> run.log
    ./dynarun -c 100000 config.out.xml.bz2 >> run.log 2>&1

    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{mft=1.0/(1.237662399*'$dens'); var=($1-mft)/mft; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "Hard Lines -: FAILED"
	    gawk 'BEGIN {mft=1.0/(1.237662399*'$dens'); print "MFT is supposed to be ",mft}'
	    bzcat output.xml.bz2 \
		| $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
		| gawk '{print "MFT is " $0}'
	    exit 1
	else
	    echo "Hard Lines -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in HardLinesTest"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 config.out.xml.bz2 output.xml.bz2 \
	tmp.xml.bz2 run.log
}

function GravityPlateTest {
    > run.log

    ./dynamod -s1 -m 22 -d 0.1  &> run.log
    ./dynarun -c 100000 config.out.xml.bz2 >> run.log 2>&1
    MFT=3.55501052762802

    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{mft='$MFT'; var=($1 - mft) / mft; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "Gravity Plate -: FAILED"
	    gawk 'BEGIN {mft='$MFT'; print "MFT is supposed to be ",mft}'
	    bzcat output.xml.bz2 \
		| $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
		| gawk '{print "MFT is " $0}'
	    exit 1
	else
	    echo "Gravity Plate -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in Gravity Plate"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 config.out.xml.bz2 output.xml.bz2 \
	tmp.xml.bz2 run.log
}

function twoDsteppedPotentialTest {
    > run.log
    ./dynamod -m 16 -z 1 --i1=2 --zero-vel 2 --rectangular-box \
	--s1 "1.0,0.1:0.9,0.2:0.8,0.3:0.7,0.4:0.6,0.5:0.5,0.6:0.4,0.7:0.3,0.8:0.2,0.9:0.1,1.0" \
	-x128 -y128 -d 1.0 -s 1  >> run.log 2>&1

    #Equilibration
    ./dynarun config.out.xml.bz2 -c 1000000 >> run.log 2>&1

    #Collection of data
    ./dynarun config.out.xml.bz2 -c 1000000 >> run.log 2>&1

    MFT="0.0419518"

    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{mft='$MFT'; var=($1-mft)/mft; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "2D-Stepped-Potential-Test -: FAILED"
	    exit 1
	else
	    echo "2D-Stepped-Potential-Test -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in 2D-Stepped-Potential-Test"
	exit 1
    fi

}

function StaticSpheresTest {
    > run.log

    ./dynarun -c 500000 static-spheres.xml >> run.log 2>&1
    
    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{mft=7.81945252098576; var=($1-mft)/mft; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "StaticSpheresTest -: FAILED"
	    exit 1
	else
	    echo "StaticSpheresTest -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in StaticSpheresTest"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 config.out.xml.bz2 output.xml.bz2 \
	tmp.xml.bz2 run.log
}

function SwingSpheresTest {
    > run.log

    ./dynarun -c 500000 swing.xml >> run.log 2>&1
    
    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{mft=0.00191272168715021; var=($1-mft)/mft; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "SwingSphereTest -: FAILED"
	    exit 1
	else
	    echo "SwingSphereTest -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in StaticSpheresTest"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 config.out.xml.bz2 output.xml.bz2 \
	tmp.xml.bz2 run.log
}

function BinaryThermalisedGranulate {
    > run.log

    ./dynarun -c 100000 -s 1 thermalisedMix.xml >> run.log 2>&1

    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val' \
	    | gawk '{mft=0.354633475131049; var=($1-mft)/mft; print ((var < 0.02) && (var > -0.02))}') != "1" ]; then
	    echo "BinaryThermalisedGranulate -: FAILED"
	    exit 1
	else
	    echo "BinaryThermalisedGranulate -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in BinaryThermalisedGranulate"
	exit 1
    fi
    
#Cleanup
    rm -Rf output.xml.bz2 config.out.xml.bz2 run.log
}

echo "INTERACTIONS+Dynamod Systems"
echo "Testing binary hard spheres, NeighbourLists and BoundedPQ's"
BinarySphereTest "Cells"
echo "Testing Lines, NeighbourLists and BoundedPQ's"
HardLinesTest
echo "Testing static spheres in gravity, NeighbourLists and BoundedPQ's"
StaticSpheresTest
#######THIS TEST IS GOOD, BUT A RECENT PATCH CHANGED THE cubic root finder, altering the result##### PLEASE RECALIBRATE
#echo "Testing static and bonded spheres in gravity, NeighbourLists and BoundedPQ's"
#SwingSpheresTest
echo "Testing *2D* stepped potential spheres, NeighbourLists and BoundedPQ's"
twoDsteppedPotentialTest

echo ""
echo "GLOBALS"
echo "Testing shearing boundary conditions with inelastic particles"
ShearingTest
echo "Testing infinite systems with neighbour lists and a 50mer polymer"
IsolatedPolymerTest
echo "Testing infinite systems with neighbour lists and gravity!"
GravityPlateTest
#echo "Testing binary spheres and the ListAndCell neighbourlist"
#BinarySphereTest "ListAndCell"

echo ""
echo "SYSTEM EVENTS"
#echo "Testing the square umbrella potential, NeighbourLists and BoundedPQ's"
#umbrella "NeighbourList"

echo ""
echo "LOCAL EVENTS"
echo "Testing local events (walls) and square wells"
wallsw "NeighbourList"
echo "Testing thermalised and normal walls in gravity with binary granulate implemented using properties"
BinaryThermalisedGranulate

echo ""
echo "ENGINE TESTING"
echo "COMPRESSION"
echo "Testing local events (walls) and square wells with a " \
    "Null compression"
wallsw "NeighbourList" "--engine 3 --growth-rate 0.0"
echo "Testing compression of polymers (very sensitive to errors in the algorithm)"
Ring_compressiontest
echo "Packing of a square well polymer into a larger system, and compressing the system"
SWpolymer_compressiontest

echo "REPLICA EXCHANGE"
echo "Testing replica exchange of hard spheres"
HS_replex_test "NeighbourList"
echo "Testing replica exchange of hard spheres with 3 threads"
HS_replex_test "NeighbourList" "-N3"
