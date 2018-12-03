#!/bin/bash
#SBATCH -t 8:00:00

# this script replicates setting up a conda environment of agalma1.1.0 with trinity/2.5.1 loaded,
# and downloads the dataset for the analyses.

source activate agalma

DATADIR=/gpfs/data/datasci/aguang/treeinform/drosophila/
export BIOLITE_RESOURCES="database=$DATADIR/agalma_trinity.sqlite"
export AGALMA_DB=$DATADIR/agalma_trinity.sqlite

# make a new database & new folders for each threshold
mkdir -p $DATADIR/
mkdir -p $DATADIR/data
mkdir -p $DATADIR/scratch
mkdir -p $DATADIR/reports

cd $DATADIR/data

EMAIL=august_guang@brown.edu
bl-sra import -c -e $EMAIL SRX246999 &
bl-sra import -c -e $EMAIL SRX247001 &
bl-sra import -c -e $EMAIL SRX247003 &
bl-sra import -c -e $EMAIL SRX054483 &
bl-sra import -c -e $EMAIL SRX054470 &
bl-sra import -c -e $EMAIL SRX054462 &
bl-sra import -c -e $EMAIL SRX054487
wait
