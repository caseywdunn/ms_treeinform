#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --array=1-9
#SBATCH --mem=16G

module load agalma/1.0.0
arr=(0.015 0.025 0.035 0.15 0.25 0.35 1.5 2.5 3.5)

# make a new database & new folders for each threshold
mkdir -p threshold-${arr[${SLURM_ARRAY_TASK_ID}]}
cd threshold-${arr[${SLURM_ARRAY_TASK_ID}]}
cp ~/repos/ms_treeinform/data/agalma-7taxa-20170207.sqlite ./threshold-${arr[${SLURM_ARRAY_TASK_ID}]}.sqlite

export AGALMA_DB=threshold-${arr[${SLURM_ARRAY_TASK_ID}]}.sqlite

#threshold=`echo "5000*10^(-${arr[${SLURM_ARRAY_TASK_ID}]})" | bc -l`
threshold=${arr[${SLURM_ARRAY_TASK_ID}]}
agalma treeinform -p 12 --id threshold-${arr[${SLURM_ARRAY_TASK_ID}]} -t $threshold
agalma homologize --id threshold-${arr[${SLURM_ARRAY_TASK_ID}]}
agalma multalign --id threshold-${arr[${SLURM_ARRAY_TASK_ID}]}
agalma genetree --id threshold-${arr[${SLURM_ARRAY_TASK_ID}]}