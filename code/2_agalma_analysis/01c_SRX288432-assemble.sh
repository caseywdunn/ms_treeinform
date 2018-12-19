#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -c 20
#SBATCH --mem=120G
#SBATCH -C intel

source activate agalma

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions
export AGALMA_DB=$DATADIR/trinity/agalma_trinity.sqlite
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"
ID=SRX288432

echo $ID

cd $DATADIR/trinity/scratch/transcriptome-22

agalma transcriptome --restart --stage 11 --id $ID

cd $DATADIR/trinity/reports
agalma report --id $ID