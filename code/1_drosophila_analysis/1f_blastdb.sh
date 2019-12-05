#!/bin/bash
#SBATCH -t 8:00:00

export PATH=/gpfs/runtime/cbc_conda/bin/:$PATH
source activate_cbc_conda
conda activate agalma

DATADIR=/gpfs/data/cbc/aguang/treeinform/drosophila/data
export BIOLITE_RESOURCES="database=$DATADIR/agalma_trinity.sqlite"
export AGALMA_DB=$DATADIR/agalma_trinity.sqlite

makeblastdb -in $DATADIR/dana-all-transcript-r1.06.fasta -out $DATADIR/SRX247001 -dbtype nucl
makeblastdb -in $DATADIR/dpse-all-transcript-r3.04.fasta -out $DATADIR/SRX054483 -dbtype nucl
makeblastdb -in $DATADIR/dsim-all-transcript-r2.02.fasta -out $DATADIR/SRX054470 -dbtype nucl
makeblastdb -in $DATADIR/dvir-all-transcript-r1.07.fasta -out $DATADIR/SRX247003 -dbtype nucl
