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

if [ ! -x $Dynarun ]; then 
    echo "Could not find dynamod, have you built it?"
fi

if [ $(whereis -b $Xml | wc -w) -lt 2 ]; then 
    Xml="xmlstarlet"
    if [ $(whereis -b $Xml | wc -w) -lt 2 ]; then 
	echo "Could not find XMLStarlet"
	exit 
    fi
fi


function HS_replex_test {
    for i in $(seq 0 2); do
	$Dynamod -m 0 -C 3 -T $(echo "0.5*$i + 0.5" | bc -l) \
	    -o config.$i.start1.xml.bz2 > /dev/null

	bzcat config.$i.start1.xml.bz2 | \
	    $Xml ed -u '/DYNAMOconfig/Simulation/Scheduler/@Type' -v "$1" \
	    | bzip2 > config.$i.start.xml.bz2
	
	rm config.$i.start1.xml.bz2
    done

    #Equilibration
    $Dynarun --engine 2 $2 -c 10000000000 -i 10 -f 600 \
	config.*.start.xml.bz2 > /dev/null

    #Production
    $Dynarun --engine 2 $2 -c 10000000000 -i 10 -f 1200 \
	-L CollisionMatrix \
	-L KEnergy config.*.end.xml.bz2 > /dev/null

    MFT1=$(bzcat output.0.xml.bz2 | $Xml sel -t -v "/OutputData/CollCounters/Totals/TotCount[@Event='CORE']/@MFT")

    MFT2=$(bzcat output.1.xml.bz2 | $Xml sel -t -v "/OutputData/CollCounters/Totals/TotCount[@Event='CORE']/@MFT")

    MFT3=$(bzcat output.2.xml.bz2 | $Xml sel -t -v "/OutputData/CollCounters/Totals/TotCount[@Event='CORE']/@MFT")

    MFT1=$(echo "$MFT1 * sqrt(0.5)" | bc -l)
    MFT3=$(echo "$MFT3 * sqrt(1.5)" | bc -l)

    AVG=$(echo "($MFT1 + $MFT2 + $MFT3) / 3.0" | bc -l)

    pass=$(echo $MFT1 $AVG | gawk '{print int(100 * (($1 / $2)-1.0))}')
    if [ $pass != 0 ]; then 
	echo "$1 HS1 Replica Exchange -: FAILED $pass"
	exit 1 	
    else
	echo -n $(echo $MFT1 $AVG | gawk '{print $1 / $2}')" "
    fi

    pass=$(echo $MFT2 $AVG | gawk '{print int(100 * (($1 / $2)-1.0))}')
    if [ $pass != 0 ]; then 
	echo "$1 HS2 Replica Exchange -: FAILED $pass"
	exit 1 	
    else
	echo -n $(echo $MFT2 $AVG | gawk '{print $1 / $AVG}')" "
    fi

    pass=$(echo $MFT3 $AVG | gawk '{print int(100 * (($1 / $2)-1.0))}')
    if [ $pass != 0 ]; then 
	echo "$1 HS3 Replica Exchange -: FAILED $pass"
	exit 1 	
    else
	echo -n $(echo $MFT3 $AVG | gawk '{print $1 / $2}')" "
    fi

    echo "$1 HS Replica Exchange -: PASSED"
    rm config.*.end.xml.bz2 config.*.start.xml.bz2 output.*.xml.bz2 *.dat replex.stats
}

function HS_compressiontest { 
    #Compresses some N=256, rho=0.5 hard spheres in a fcc, which
    #should compress into a maximum packed FCC crystal quite rapidly
    $Dynamod -m0 -d0.5 -C 4 &> run.log

    #Set the scheduler
    bzcat config.out.xml.bz2 | \
	$Xml ed -u '/DYNAMOconfig/Simulation/Scheduler/@Type' -v "$1" \
	| bzip2 > tmp.xml.bz2

    #Run the simulation
    $Dynarun -c 100000 --engine 3 --growth-rate 100.0 $2 tmp.xml.bz2 &> run2.log
    cat run2.log >> run.log
    rm run2.log
    $Dynarun -c 1 config.out.xml.bz2 &> run2.log
    cat run2.log >> run.log
    rm run2.log

    #The 0.xxx is the fractional error, this will fail if the density
    #is higher than the max pack frac, and fail if its fractionally
    #more than 0.xxxx under max pack frac
    pass=$(echo "0.0001 > (-("$(bzcat output.xml.bz2 | $Xml sel -t -v \
	'/OutputData/Misc/Density/@val')"/sqrt(2) - 1))" | bc -l)
    
    if [ $pass != 1 ]; then
	echo "$1 HS compression -: FAILED"
	exit 1 
    else
	echo "$1 compression -: PASSED"
    fi
    
    
    #Cleanup 
    rm -Rf config.out.xml.bz2 \
	output.xml.bz2 tmp.xml.bz2 run.log \
	testresult.dat correct.dat
    }

function wallsw_bailout { 
    echo "$1 WallSW -: FAILED"
    exit 1 
}

function wallsw {
#Just tests the square shoulder interaction between two walls
    bzcat wallsw.xml.bz2 | \
	$Xml ed -u '/DYNAMOconfig/Simulation/Scheduler/@Type' -v "$1" \
	| bzip2 > tmp.xml.bz2
    
    #Any multiple of 5 will/should always give the original configuration
    #doing 9995 as this stops any 2 periodicity
    $Dynarun -c 9995 $2 tmp.xml.bz2 &> run.log
    
    bzcat config.out.xml.bz2 | \
	$Xml sel -t -c '/DYNAMOconfig/ParticleData' > testresult.dat

    bzcat wallsw.xml.bz2 | \
	$Xml sel -t -c '/DYNAMOconfig/ParticleData' > correct.dat

    diff testresult.dat correct.dat &> /dev/null \
	|| wallsw_bailout $1
    
    echo "$1 WallSW -: PASSED"

    #Cleanup
    rm -Rf config.out.xml.bz2 output.xml.bz2 tmp.xml.bz2 run.log \
	testresult.dat correct.dat
}

function cannon {
    #collision cannon test
    bzcat cannon.xml.bz2 | \
	$Xml ed -u '/DYNAMOconfig/Simulation/Scheduler/@Type' -v "$1" \
	| bzip2 > tmp2.xml.bz2

    bzcat tmp2.xml.bz2 | \
	$Xml ed -u '/DYNAMOconfig/Simulation/Scheduler/Sorter/@Type' -v "$2" \
	| bzip2 > tmp.xml.bz2
    
    rm tmp2.xml.bz2
    
    $Dynarun -c 1000 tmp.xml.bz2 &> run.log
    
    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 | \
	    $Xml sel -t -v '/OutputData/Misc/totMeanFreeTime/@val') != "3" ]; then
	    echo "$1 Cannon -: FAILED"
	    exit 1
	else
	    echo "$1 Cannon -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in $1 cannon test"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 output.xml.bz2 tmp.xml.bz2 run.log
}

function linescannon {
    #collision cannon test
    bzcat config.lines-cannon.xml.bz2 | \
	$Xml ed -u '/DYNAMOconfig/Simulation/Scheduler/@Type' -v "$1" \
	| bzip2 > tmp2.xml.bz2

    bzcat config.lines-cannon.xml.bz2 | \
	$Xml ed -u '/DYNAMOconfig/Simulation/Scheduler/Sorter/@Type' -v "$2" \
	| bzip2 > tmp.xml.bz2
    
    rm tmp2.xml.bz2
    
    $Dynarun -c 10 tmp.xml.bz2 &> run.log
    
    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/Misc/SimLength/@Time' \
	    | gawk '{printf "%.0f",$1}') != "40" ]; then
	    echo "$1 Lines Cannon -: FAILED"
	    exit 1
	else
	    echo "$1 Lines Cannon -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in $1 cannon test"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 output.xml.bz2 tmp.xml.bz2 run.log
}

function ThermostatTest {
    #Testing the Andersen thermostat holds the right temperature
    > run.log

    $Dynamod -m 0 -r 0.1 -T 1.0 -o tmp.xml.bz2 &> run.log    
    $Dynarun -c 100000 tmp.xml.bz2 >> run.log 2>&1
    $Dynarun -c 100000 config.out.xml.bz2 >> run.log 2>&1
    $Dynarun -c 500000 config.out.xml.bz2 -L KEnergy >> run.log 2>&1
    
    if [ -e output.xml.bz2 ]; then
	if [ $(bzcat output.xml.bz2 \
	    | $Xml sel -t -v '/OutputData/KEnergy/T/@val' \
	    | gawk '{printf "%.1f",$1}') != "1.0" ]; then
	    echo "Thermostat -: FAILED"
	    exit 1
	else
	    echo "Thermostat -: PASSED"
	fi
    else
	echo "Error, no output.0.xml.bz2 in $1 cannon test"
	exit 1
    fi
    
#Cleanup
    rm -Rf config.end.xml.bz2 config.out.xml.bz2 output.xml.bz2 \
	tmp.xml.bz2 run.log
}

echo "SCHEDULER AND SORTER TESTING"
echo "Testing basic system, zero + infinite time events, hard spheres, PBC, Dumb Scheduler, CBT"
cannon "Dumb" "CBT"
echo "Testing basic system, zero + infinite time events, hard spheres, PBC, Dumb Scheduler, boundedPQ"
cannon "Dumb" "BoundedPQ"
echo "Testing basic system, zero + infinite time events, hard sphere, PBC, Neighbour lists + scheduler, globals, CBT"
cannon "NeighbourList" "CBT"
echo "Testing basic system, zero + infinite time events, hard sphere, PBC, Neighbour lists + scheduler, globals, boundedPQ"
cannon "NeighbourList" "BoundedPQ"

echo ""
echo "SYSTEM EVENTS"
echo "Testing the Andersen Thermostat, NeighbourLists and BoundedPQ's"
ThermostatTest

echo ""
echo "LOCAL EVENTS"
echo "Testing local events (walls) and square wells"
wallsw "NeighbourList"

echo ""
echo "INTERACTIONS"
echo "Testing Lines, NeighbourLists and BoundedPQ's"
linescannon "NeighbourList" "BoundedPQ"


echo ""
echo "ENGINE TESTING"
echo "Testing local events (walls) and square wells with a " \
    "Null compression"
wallsw "NeighbourList" "--engine 3 --growth-rate 0.0"
echo "Testing compression of hard spheres"
HS_compressiontest "NeighbourList"

#echo "Testing replica exchange of hard spheres"
#HS_replex_test "NeighbourList"
#
#echo ""
#echo "THREADING TESTING"
#echo "Testing replica exchange with 2 threads"
#HS_replex_test "NeighbourList" "-N2"
