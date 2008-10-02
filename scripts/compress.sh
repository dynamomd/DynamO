#!/bin/bash
cp config.start.xml.bz2 config.out.xml.bz2
for eta in $(seq 0.05 0.05 0.7); do
    echo $(date) Started pack frac $eta
    ~/git_tree/code/dynamo/mdrun config.out.xml.bz2 -D2 --growth-rate 0.01 --target-pack-frac $eta -c 1000000000 #> /dev/null
    cp config.out.xml.bz2 config.eta$eta.xml.bz2
done

#cp config.start.xml.bz2 config.out.xml.bz2
#for eta in $(seq 0.05 0.05 0.9); do
#    echo $(date) Started pack frac $eta
#    ~/git_tree/code/dynamo/configmod config.eta$eta.xml.bz2 -o config.eta$eta.xml.bz2 -Z -r 1
#    echo "Running equilibration now"
#    ~/git_tree/code/dynamo/mdrun config.eta$eta.xml.bz2 -o  config.eta$eta.xml.bz2 -c 5000000
#    cp config.out.xml.bz2 config.eta$eta.xml.bz2
#done
