#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --array=1-10
#SBATCH --mem=16G

module load agalma/1.0.0

# make a new database & new folders for each threshold
mkdir -p threshold-${SLURM_ARRAY_TASK_ID}
cd threshold-${SLURM_ARRAY_TASK_ID}
cp ~/repos/ms_treeinform/data/agalma-7taxa-20170207.sqlite ./threshold-${SLURM_ARRAY_TASK_ID}.sqlite

export AGALMA_DB=threshold-${SLURM_ARRAY_TASK_ID}.sqlite

threshold=`echo "5000*10^(-${SLURM_ARRAY_TASK_ID})" | bc -l`
agalma treeinform -p 12 --id threshold-${SLURM_ARRAY_TASK_ID} -t $threshold
agalma homologize --id threshold-${SLURM_ARRAY_TASK_ID}
agalma multalign --id threshold-${SLURM_ARRAY_TASK_ID}
agalma genetree --id threshold-${SLURM_ARRAY_TASK_ID}