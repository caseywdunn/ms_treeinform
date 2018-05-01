#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 8
#SBATCH -c 16
#SBATCH --mem=60G
#SBATCH -C intel
#SBATCH --array=1-9

source activate agalma

additional=( 0.01 0.02 0.035 0.15 0.25 0.35 1.5 2.5 3.5 )

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity
export AGALMA_DB=$DATADIR/agalma_trinity.sqlite
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

set -e

cd $DATADIR/scratch/
threshold=${additional[$SLURM_ARRAY_TASK_ID-1]}
#-p option is for the previous genetree run to link treeinform to.
#In our case it is 27 but if you haven't run into errors it should be 12.
agalma treeinform -p 27 --id threshold-${threshold} -t $threshold
agalma homologize --id threshold-${threshold}
agalma multalign --id threshold-${threshold}
agalma genetree --id threshold-${threshold}