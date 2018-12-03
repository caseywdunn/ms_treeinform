#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 8
#SBATCH -c 16
#SBATCH --mem=60G
#SBATCH -C intel

source activate agalma

DATADIR=/gpfs/scratch/aguang/ms_treeinform/data/drosophila_analysis/trinity/
export AGALMA_DB=$DATADIR/agalma_trinity.sqlite
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

set -e

cd $DATADIR/scratch
agalma homologize --id TrinityPhylogeny
agalma multalign --id TrinityPhylogeny
agalma genetree --id TrinityPhylogeny