#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 8
#SBATCH -c 16
#SBATCH --mem=60G
#SBATCH -C intel

source activate agalma

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/corset
export AGALMA_DB=/users/aguang/data/aguang/corset/agalma_corset.sqlite
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

set -e

mkdir -p $DATADIR/scratch
cd $DATADIR/scratch
agalma homologize --id CorsetPhylogeny
agalma multalign --id CorsetPhylogeny
agalma genetree --id CorsetPhylogeny