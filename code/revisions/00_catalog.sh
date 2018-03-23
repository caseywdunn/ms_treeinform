#!/bin/bash
#SBATCH -t 1:00:00

# this script replicates setting up a conda environment of agalma1.1.0 with trinity/2.5.1 loaded,
# and downloads the dataset for the analyses. At this time (February 19, 2018) this is the
# latest version of trinity.

#conda create -n agalma -c dunnlab agalma

source activate agalma

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions
export BIOLITE_RESOURCES="database=$DATADIR/trinity/agalma_trinity.sqlite"
export AGALMA_DB=$DATADIR/trinity/agalma_trinity.sqlite

# make a new database & new folders for each threshold
mkdir -p $DATADIR/trinity/data
mkdir -p $DATADIR/trinity/scratch
mkdir -p $DATADIR/trinity/reports

cd $DATADIR/trinity/data

EMAIL=august_guang@brown.edu
bl-sra import -c -e $EMAIL SRX288432 &
bl-sra import -c -e $EMAIL SRX288431 &
bl-sra import -c -e $EMAIL SRX288430 &
bl-sra import -c -e $EMAIL SRX288285 &
bl-sra import -c -e $EMAIL SRX288276 &
wait

wget ftp://ftp.jgi-psf.org/pub/JGI_data/Nematostella_vectensis/v1.0/annotation/transcripts.Nemve1FilteredModels1.fasta.gz
gunzip transcripts.Nemve1FilteredModels1.fasta.gz
mv transcripts.Nemve1FilteredModels1.fasta Nematostella_vectensis_rna.fa
agalma catalog insert --id "JGI_NEMVEC" --paths "Nematostella_vectensis_rna.fa" --species "Nematostella vectensis" --ncbi_id "45351" --itis_id "52498" --library_type "genome" --note "Gene predictions from genome sequencing"

wget http://ftp.ncbi.nih.gov/genomes/Hydra_magnipapillata/RNA/rna.fa.gz
gunzip rna.fa.gz
mv rna.fa Hydra_magnipapillata_rna.fa 
agalma catalog insert --id "NCBI_HYDMAG" --paths "Hydra_magnipapillata_rna.fa" --species "Hydra magnipapillata" --ncbi_id "6085" --itis_id "50845" --library_type "genome" --note "Gene predictions from genome sequencing"
