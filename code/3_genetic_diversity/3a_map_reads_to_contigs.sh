#!/bin/bash
#SBATCH -t 8:00:00

# Assumes you have reads and contigs from running code in drosophila_assembly
# and agalma_assembly already or through an alternative method.

source activate agalma

DATAIDR=/gpfs/data/aguang/treeinform/data/genetic_diversity
export BIOLITE_RESOURCES="database=$DATADIR/trinity/agalma_trinity.sqlite"
export AGALMA_DB=$DATADIR/trinity/agalma_trinity.sqlite

# map reads to contigs with BWA
bwa index ref.fa
bwa mem ref.fa read1.fq read2.fq > aln-pe.sam

# discard contigs with coverage less than average 2.5x/individual
# defined as total length of mapped reads / contig length
samtools view -bS aln-pe.sam > aln-pe.bam
samtools depth aln-pe.bam > coverage.txt
# filter by coverage