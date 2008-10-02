#!/bin/bash
echo Testing for mdrun executable at $MDRUN
if [ ! -x "$MDRUN" ]; then
    echo "mdrun executable not found!"
    exit
fi

if [ ! -d $OUTPUT_DIR ]; then
    mkdir $OUTPUT_DIR
fi

cd $OUTPUT_DIR

if [ ! -d ./equilibration ]; then
    mkdir equilibration
fi

cd equilibration

if [ ! -d ./tmp ]; then
    mkdir ./tmp
fi
cd tmp

##Equilibration loop
for RUN in $(seq 1 $NERuns); do
    if [ -d "../Run$RUN" ]; then
	echo "Output directory for equilibration run $RUN already exists, skipping running the system"
    else
        ##If this is the first equilibration, the files are stored in the initial configuration folder
	if [ $RUN == 1 ]; then
	    if [ ! -d "../../initial" ]; then
		echo "Initial configuration directory doesn't exist! Required to start the simulations"
		exit
	    fi
	    
	    FILES="../../initial/config.out.xml.bz2"
	    if [ ! -e $FILES ]; then
		FILES=$(ls ../../initial/config.*.end.xml.bz2)
	    fi
	else
	    FILES="../Run"$(($RUN - 1))"/config.out.xml.bz2"
	    if [ ! -e $FILES ]; then
		FILES=$(ls ../Run$(($RUN - 1))/config.*.end.xml.bz2)
	    fi
	fi
	
	echo Config files are $FILES 
	
	echo "Running Equilibration trajectory "$RUN
	
	echo "Running the command $MDRUN $MDRUN_EOPTS $FILES || exit \$?"
	$MDRUN $MDRUN_EOPTS $FILES || exit $?

	#Move the output and final config files
	mkdir ../Run$RUN
	mv * ../Run$RUN
    fi
done

cd ..

if [ -d ./tmp ]; then
    rmdir tmp
fi

cd ..

if [ ! -d ./production ]; then
    mkdir production
fi

cd production

if [ ! -d ./tmp ]; then
    mkdir ./tmp
fi

cd tmp
for RUN in $(seq 1 $NRuns); do
    if [ -d "../Run$RUN" ]; then
	echo "Output directory for production run $RUN already exists, skipping running the system"
    else
	FILES=""
        ##If this is the first equilibration, the files are stored in the initial configuration folder
	if [ $RUN == 1 ]; then
	    if [ ! -d "../../equilibration/Run$NERuns" ]; then
		if [ ! -d "../../initial" ]; then
		    echo "Initial configuration directory doesn't exist and neither does the output of the equilibration!"
		    exit
		else
		    echo "Found the initial configuration folder, guessing that equilibration has been skipped"

		    FILES="../../initial/config.out.xml.bz2"
		    if [ ! -e $FILES ]; then
			FILES=$(ls ../../initial/config.*.end.xml.bz2)
		    fi
		fi
	    else
		FILES="../../equilibration/Run$NERuns/config.out.xml.bz2"
		if [ ! -e $FILES ]; then
		    FILES=$(ls ../../equilibration/Run$NERuns/config.*.end.xml.bz2)
		fi
	    fi
	else
	    FILES="../Run$(($RUN - 1))/config.out.xml.bz2"
	    if [ ! -e $FILES ]; then
		FILES=$(ls ../Run$(($RUN - 1))/config.*.end.xml.bz2)
	    fi
	fi
	
	echo "Running trajectory "$RUN

	echo "Running the command $MDRUN $MDRUN_OPTS $FILES || exit \$?"
	$MDRUN $MDRUN_OPTS $FILES || exit $?
	
	#Move the output and final config files
	mkdir ../Run$RUN
	mv * ../Run$RUN
    fi
done

cd ..
if [ -d tmp ]; then
    rmdir tmp
fi

cd ..
