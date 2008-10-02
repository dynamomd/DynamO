#!/bin/bash

#How many runs to do
PresentWD=$(pwd)
MDRUN_DIR=$(pwd)"/dynamo"
MDRUN="$MDRUN_DIR/mdrun"

NRuns=8
NColls=1000000
RANDOMSEED=1
DYNAMICS=1
#1 = SLLOD, 2 = Elastic, 3 = Growth, 4= LE + Newtonian
SCHEDULER=1 
#1 = Cellular, 2 = single

cd $MDRUN_DIR
svn update
make -j9 mdrun
cd $PresentWD

for FCC in $(seq 10 30); do
#FCC=7
#BCC=5
#SC=5

XCELLS=$FCC
YCELLS=$FCC
ZCELLS=$FCC

#for DENS in $(seq 0.5 0.1 1.01); do
#    for ELAS in $(seq 0.5 0.1 0.9); do
DENS=0.5
ELAS=0.5
	
NAME="FCC"$FCC".x"$XCELLS".d"$DENS".e"$ELAS
OUTPUT_DIR="$NAME"

qsub -N $NAME -q Xeon -l nodes=1:ppn=1 -o $OUTPUT_DIR/Std.out -e $OUTPUT_DIR/Std.err -v DYNAMICS=$DYNAMICS,SCHEDULER=$SCHEDULER,SEED=$RANDOMSEED,XCELLS=$XCELLS,YCELLS=$YCELLS,ZCELLS=$ZCELLS,NRuns=$NRuns,DENS=$DENS,ELAS=$ELAS,NColls=$NColls,MDRUN=$MDRUN,OUTPUT_DIR=$OUTPUT_DIR,FCC=$FCC,BCC=$BCC,SC=$SC,NAME=$NAME dynamo/scripts/qsub_script.sh

#done
#done
done