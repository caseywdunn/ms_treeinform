#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -c 20
#SBATCH --mem=120G
#SBATCH -C intel
#SBATCH --array=1-5

source activate agalma

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions
export AGALMA_DB=$DATADIR/trinity/agalma_trinity.sqlite
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

IDS=(
SRX288276
SRX288285
SRX288430
SRX288431
SRX288432
)

ID=${IDS[$SLURM_ARRAY_TASK_ID-1]}
echo $ID

cd $DATADIR/trinity/scratch

agalma transcriptome --id $ID

cd $DATADIR/trinity/reports
agalma report --id $ID