#!/bin/bash

PresentWD=$(pwd)
MDRUN_DIR=$(pwd)"/dynamo"
MDRUN=$MDRUN_DIR"/mdrun"

#How many runs to do
NRuns=8
NColls=1000000
SEED=1

cd $MDRUN_DIR
svn update
HORACE=YES make -j9 mdrun
cd $PresentWD

#Simulation settings, all passed automagically to bsub script
FCC=18
XCELLS=7
YCELLS=$FCC
ZCELLS=$FCC
#BCC=5
#SC=5

#for ELAS in $(seq 0.5 0.1 1.01); do
ELAS=0.9	

for DENS in $(seq 0.5 0.1 1.01); do
    NAME="FCC"$FCC".x"$XCELLS".d"$DENS".e"$ELAS
    
    OUTPUT_DIR="$NAME"
    
    mkdir $OUTPUT_DIR
    
    echo -e "PresentWD=$PresentWD\nMDRUN_DIR=$MDRUN_DIR\nMDRUN=$MDRUN\nNRuns=$NRuns\nNColls=$NColls\nSEED=$SEED\nFCC=$FCC\nXCELLS=$XCELLS\nYCELLS=$YCELLS\nZCELLS=$ZCELLS\nBCC=$BCC\nSC=$SC\nELAS=$ELAS\nDENS=$DENS\nNAME=$NAME\nOUTPUT_DIR=$OUTPUT_DIR\n" > script.tmp.sh
    
    cat dynamo/scripts/qsub_script.sh | tail -n +2 >> script.tmp.sh
    
    bsub -J $NAME -o $OUTPUT_DIR/Std.out -e $OUTPUT_DIR/Std.err -W 24:00 -n 1 -q single -B < script.tmp.sh    

done

rm script.tmp.sh
