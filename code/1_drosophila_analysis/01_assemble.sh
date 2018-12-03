#!/bin/bash
#SBATCH -t 4-00:12:00
#SBATCH -n 8
#SBATCH --mem=48G

source activate agalma

DATADIR=/gpfs/scratch/aguang/ms_treeinform/data/drosophila_analysis
export BIOLITE_RESOURCES="database=$DATADIR/trinity/agalma_trinity.sqlite"
export AGALMA_DB=$DATADIR/trinity/agalma_trinity.sqlite

cd $DATADIR/trinity/data
agalma catalog insert --id SRX246999 --paths SRR768436_1.fastq SRR768436_2.fastq --species "Drosophila yakuba" --ncbi_id 7245
agalma catalog insert --id SRX247001 --paths SRR768438_1.fastq SRR768438_2.fastq --species "Drosophila ananassae" --ncbi_id 7217
agalma catalog insert --id SRX247003 --paths SRR768440_1.fastq SRR768440_2.fastq --species "Drosophila virilis" --ncbi_id 7244
agalma catalog insert --id SRX054483 --paths SRR166831_1.fastq SRR166831_2.fastq --species "Drosophila pseudoobscura" --ncbi_id 7237
agalma catalog insert --id SRX054470 --paths SRR166818_1.fastq SRR166818_2.fastq --species "Drosophila simulans" --ncbi_id 7240
agalma catalog insert --id SRX054462 --paths SRR166810_1.fastq SRR166810_2.fastq --species "Drosophila melanogaster" --ncbi_id 7227
agalma catalog insert --id SRX054487 --paths SRR166835_1.fastq SRR166835_2.fastq  --species "Drosophila mojavensis" --ncbi_id 7230

cd $DATADIR/trinity/scratch
agalma -t 8 -m 48G transcriptome --id SRX246999
agalma -t 8 -m 48G transcriptome --id SRX247001
agalma -t 8 -m 48G transcriptome --id SRX247003
agalma -t 8 -m 48G transcriptome --id SRX054483
agalma -t 8 -m 48G transcriptome --id SRX054470
agalma -t 8 -m 48G transcriptome --id SRX054462
agalma -t 8 -m 48G transcriptome --id SRX054487

cd $DATADIR/trinity/reports
agalma report --id SRX246999
agalma report --id SRX247001
agalma report --id SRX247003
agalma report --id SRX054483
agalma report --id SRX054470
agalma report --id SRX054462
agalma report --id SRX054487