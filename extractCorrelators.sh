#!/bin/bash
file="output.xml.bz2"
file="peek.data.xml.bz2"


bzcat $file \
    | xml sel -t -v "/OutputData/EinsteinCorrelator[@name='ThermalDiffusionE']"\
    | gawk '{if (NF) print $1,($2+$3+$4)/3.0}' > thermalDiff.dat

bzcat $file \
    | xml sel -t \
    -v "/OutputData/EinsteinCorrelator[@name='ThermalConductivityE']"\
    | gawk '{if (NF) print $1,($2+$3+$4)/3.0}' > thermalCond.dat

bzcat $file \
    | xml sel -t -v "/OutputData/Correlator[@name='MutualDiffusionGK']"\
    | gawk '{if (NF) print $1,($2+$3+$4)/3.0}' > mutualDiff.GK.dat

bzcat $file \
    | xml sel -t -v "/OutputData/EinsteinCorrelator[@name='MutualDiffusionE']"\
    | gawk '{if (NF) print $1,($2+$3+$4)/3.0}' > mutualDiff.E.dat

bzcat $file \
    | xml sel -t \
    -v "/OutputData/EinsteinCorrelator[@name='ViscosityE']"\
    | gawk '{if (NF) print $1,($3+$4+$5+$7+$8+$9)/6.0}' > shearVisc.dat

bzcat $file \
    | xml sel -t \
    -v "/OutputData/EinsteinCorrelator[@name='ViscosityE']"\
    | gawk '{if (NF) print $1,($2+$6+$10)/3.0}' > bulkVisc.dat

