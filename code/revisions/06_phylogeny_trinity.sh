#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH -n 8
#SBATCH --mem=16G

module load conda
source activate agalma

DATADIR=/gpfs/scratch/aguang/treeinform/revisions/trinity
export AGALMA_DB=$DATADIR/agalma_trinity.sqlite
agalma -t 8 -m 16G

cd $DATADIR/trinity/scratch
agalma homologize --id TrinityPhylogeny
agalma multalign --id TrinityPhylogeny
agalma genetree --id TrinityPhylogeny