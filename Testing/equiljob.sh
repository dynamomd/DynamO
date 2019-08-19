#!/bin/bash
#SBATCH --job-name DynamO_test
#SBATCH --partition=laird,cebc,sixhour
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL

module load compiler/gcc/8.3

dynamod -m0  --i1=2 -z 1 --zero-vel 2 --rectangular-box -d 0.3 -C 50 -o config.out.xml

./add_wall.py -f config.out.xml -r 15 -N 500 -L 100 -o config.fixed.xml

dynarun --out-data-file output.test.xml -o config.test.xml -p 1000 -c 10000 config.fixed.xml -L RadialDistribution:BinWidth=0.1,Length=300,rdfpairs="Wall Bulk"
