#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 8
#SBATCH -c 16
#SBATCH --mem=60G
#SBATCH -C intel

source activate agalma

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity
export AGALMA_DB=$DATADIR/agalma_trinity.sqlite
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

set -e

SLURM_ARRAY_TASK_ID=4
cd $DATADIR/scratch/treeinform-$((SLURM_ARRAY_TASK_ID+27))/multalign-49
threshold=`echo "5000*10^(-${SLURM_ARRAY_TASK_ID})" | bc -l`
#-p option is for the previous genetree run to link treeinform to.
#In our case it is 27 but if you haven't run into errors it should be 12.
#agalma treeinform -p 27 --id threshold-${SLURM_ARRAY_TASK_ID} -t $threshold #--restart --stage 2
#agalma homologize --id threshold-${SLURM_ARRAY_TASK_ID}
agalma multalign --id threshold-${SLURM_ARRAY_TASK_ID} --restart --stage 2
agalma genetree --id threshold-${SLURM_ARRAY_TASK_ID}