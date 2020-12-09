#!/bin/bash
#SBATCH --job-name DynamO_compile
#SBATCH --partition=laird,sixhour,cebc
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --mem=5gb
#SBATCH --time=6:00:00
##SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load compiler/gcc/8.3
module load cmake
module load anaconda/4.7

#rm -r build
#mkdir build
cd build
export BOOST_ROOT=$WORK/boost_1_64_0/
export BOOST_LIBRARYDIR=$WORK/boost_1_64_0/stage/lib
#cmake ../ -DCMAKE_BUILD_TYPE="release"
make
