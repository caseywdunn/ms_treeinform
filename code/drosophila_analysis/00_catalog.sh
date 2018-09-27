#!/bin/bash
#SBATCH -t 8:00:00

# this script replicates setting up a conda environment of agalma1.1.0 with trinity/2.5.1 loaded,
# and downloads the dataset for the analyses. At this time (February 19, 2018) this is the
# latest version of trinity.

conda create -n agalma -c dunnlab agalma

source activate agalma

DATADIR=/gpfs/scratch/aguang/ms_treeinform/data/drosophila_analysis
export BIOLITE_RESOURCES="database=$DATADIR/trinity/agalma_trinity.sqlite"
export AGALMA_DB=$DATADIR/trinity/agalma_trinity.sqlite

# make a new database & new folders for each threshold
mkdir -p $DATADIR/trinity
mkdir -p $DATADIR/trinity/data
mkdir -p $DATADIR/trinity/scratch
mkdir -p $DATADIR/trinity/reports

cd $DATADIR/trinity/data

EMAIL=august_guang@brown.edu
bl-sra import -c -e $EMAIL SRX246999 &
bl-sra import -c -e $EMAIL SRX247001 &
bl-sra import -c -e $EMAIL SRX247003 &
bl-sra import -c -e $EMAIL SRX054483 &
bl-sra import -c -e $EMAIL SRX054470 &
bl-sra import -c -e $EMAIL SRX054462 &
bl-sra import -c -e $EMAIL SRX054487
wait