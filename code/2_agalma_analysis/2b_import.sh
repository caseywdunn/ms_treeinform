#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -n 8
#SBATCH --mem=20G
#SBATCH --array=1-2

source activate agalma

DATADIR=/gpfs/data/datasci/aguang/treeinform/agalma/trinity
export AGALMA_DB=$DATADIR/agalma_trinity.sqlite
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

IDS=( JGI_NEMVEC NCBI_HYDMAG )
ID=${IDS[$SLURM_ARRAY_TASK_ID-1]}
echo $ID

cd $DATADIR/scratch
agalma import --id $ID
agalma translate --id $ID
