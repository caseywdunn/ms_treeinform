#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 8
#SBATCH -c 16
#SBATCH --mem=60G
#SBATCH -C intel
#SBATCH --array=1-10
# array=10 is no genes reassigned

source activate agalma

DATADIR=/gpfs/data/cbc/aguang/treeinform/drosophila
export AGALMA_DB=$DATADIR/agalma_trinity.sqlite
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

set -e

cd $DATADIR/scratch
threshold=`echo "5000*10^(-${SLURM_ARRAY_TASK_ID})" | bc -l`
#-p option is for the previous genetree run to link treeinform to.
#In our case it is 10
agalma treeinform -p 10 --id threshold-${SLURM_ARRAY_TASK_ID} -t $threshold #--restart --stage 2
agalma homologize --id threshold-${SLURM_ARRAY_TASK_ID}
agalma multalign --id threshold-${SLURM_ARRAY_TASK_ID}
agalma genetree --id threshold-${SLURM_ARRAY_TASK_ID}
