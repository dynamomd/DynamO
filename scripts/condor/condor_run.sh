#!/bin/bash
DYNARUN="/home/mjki2mb2/dynamo/bin/dynarun"
DYNAMOD="/home/mjki2mb2/dynamo/bin/dynamod"
NERuns=1
NRuns=5
MDRUN_EOPTS='-c 10000'
MDRUN_OPTS='-c 100000'

#Clear the DAGMan job file

function worker {
    if [ ! -d $OUTPUT_DIR ]; then
	echo "Could not find the output directory, $OUTPUT_DIR"
	exit
    fi

    if [ ! -e $OUTPUT_DIR/Equilibration/Run0/config.out.xml ]; then
	echo "Could not find the initial configuration"\
	    "$OUTPUT_DIR/Equilibration/Run0/config.out.xml"
	exit
    fi
    mkdir -p $OUTPUT_DIR/Production

    for i in $(seq 1 $NERuns); do
	mkdir -p $OUTPUT_DIR/Equilibration/Run$i
    done

    for i in $(seq 1 $NRuns); do
	mkdir -p $OUTPUT_DIR/Production/Run$i
    done
    
    for i in $(seq 1 $NERuns); do
	cat condor.sub | sed 's,#INITDIR#,../Equilibration/Run'$i',' > $OUTPUT_DIR/condor/E$i.sub
	ARG="$MDRUN_EOPTS --uncompressed $OUTPUT_DIR/Equilibration/Run$((i-1))/config.out.xml -o $OUTPUT_DIR/Equilibration/Run$i/config.out.xml"

	echo "JOB E$i E$i.sub" >> $OUTPUT_DIR/condor/job.dag
	echo "VARS E$i DYNARUN=\"$DYNARUN\" ARG=\"$ARG\" OUTPUTDIR=\"../Equilibration/Run$i\"" \
	    >> $OUTPUT_DIR/condor/job.dag
    done

    echo "SCRIPT POST E$NERuns $PWD/endofequil.sh $NERuns $OUTPUT_DIR" >> $OUTPUT_DIR/condor/job.dag
    
    for i in $(seq 1 $NRuns); do
	cat condor.sub | sed 's,#INITDIR#,../Production/Run'$i',' > $OUTPUT_DIR/condor/P$i.sub
	ARG="$MDRUN_OPTS --uncompressed $OUTPUT_DIR/Production/Run$((i-1))/config.out.xml -o $OUTPUT_DIR/Production/Run$i/config.out.xml"

	echo "JOB P$i P$i.sub" >> $OUTPUT_DIR/condor/job.dag
	echo "VARS P$i DYNARUN=\"$DYNARUN\" ARG=\"$ARG\" OUTPUTDIR=\"$OUTPUT_DIR/Production/Run$i\"" \
	    >> $OUTPUT_DIR/condor/job.dag
    done

    for i in $(seq 1 $((NERuns - 1))); do
	echo "PARENT E$i CHILD E$((i+1))" >> $OUTPUT_DIR/condor/job.dag
    done

    echo "PARENT E$NERuns CHILD P1" >> $OUTPUT_DIR/condor/job.dag

    for i in $(seq 1 $((NRuns - 1))); do
	echo "PARENT P$i CHILD P$((i+1))" >> $OUTPUT_DIR/condor/job.dag
    done
    
    echo "DOT dag.dot UPDATE" >> $OUTPUT_DIR/condor/job.dag

    cd $OUTPUT_DIR/condor/
    condor_submit_dag job.dag
}


NAME="condorTest"
OUTPUT_DIR="$PWD/$NAME"
mkdir -p $OUTPUT_DIR/condor
> $OUTPUT_DIR/condor/job.dag
mkdir -p $OUTPUT_DIR/Equilibration/Run0
mkdir -p $OUTPUT_DIR/Production/Run0
$DYNAMOD -m 0 -o $OUTPUT_DIR/Equilibration/Run0/config.out.xml --uncompressed
worker
