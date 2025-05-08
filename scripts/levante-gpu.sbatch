#!/bin/bash
  
#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=%x.%j.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

# input file
IN_FILE=$(pwd)/tasks/dbg.nc
# output file
OUT_FILE=output.nc
# reference file associated with input file
REF_FILE=$(pwd)/reference_results/seq_dbg_double.nc

# setup env
. scripts/levante-setup.sh nvidia cpu

# build the code
PERCISION=single
. scripts/build_cpu_$PERCISION.sh

module load cdo

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 

. scripts/check-correctness.sh

# run the executable

BIN=$(pwd)/$BUILD/bin/graupel

$BIN $IN_FILE $OUT_FILE
check $OUT_FILE $REF_FILE

$BIN $(pwd)/tasks/11k.nc $OUT_FILE
check $OUT_FILE $(pwd)/reference_results/seq_11k_$PERCISION.nc

$BIN $(pwd)/tasks/20k.nc $OUT_FILE
check $OUT_FILE $(pwd)/reference_results/seq_20k_$PERCISION.nc

echo "Job done!"
