#!/bin/bash
MDRUN="/home/marcus/git_tree/code/dynamo/mdrun"
NERuns=1
NRuns=40
MDRUN_EOPTS='--replex -i 0.1 -f 100 -c 1000000000 -N 4 -E -L UEnergy -L Misc'
MDRUN_OPTS='--replex -i 0.1 -f 100 -c 1000000000 -N 4 -E -L Misc -L UEnergy -L CollisionMatrix -L KEnergy -L IntEnergyHist'
PE="-pe shm 4" #Using the parallel memory environment
function worker {
    if [ ! -d $OUTPUT_DIR ]; then
	echo "Could not find output dir $OUTPUT_DIR"
	exit
    fi
    	
    #-hold_jid $NAME # to introduce backup scripts
    qsub $PE -o $OUTPUT_DIR/terminal.out -S /bin/bash -j y -cwd -N $NAME -v MDRUN=$MDRUN,NRuns=$NRuns,NERuns=$NERuns,MDRUN_EOPTS="$MDRUN_EOPTS",MDRUN_OPTS="$MDRUN_OPTS",NEColls=$NEColls,OUTPUT_DIR="$OUTPUT_DIR" ~/git_tree/code/dynamo/scripts/qsub_script.sh 
}

for eta in $(seq -f "%.2f" 0.05 0.05 0.7); do
    OUTPUT_DIR="$PWD/eta$eta"
    NAME="Homo.eta$eta"
    worker
done
