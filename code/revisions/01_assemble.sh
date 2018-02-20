#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=24G

module load conda
source activate agalma

DATADIR=/gpfs/scratch/bguang/treeinform/revisions

# make a new database & new folders for each threshold
mkdir -p $DATADIR/trinity/data
mkdir -p $DATADIR/trinity/scratch
mkdir -p $DATADIR/trinity/reports

export AGALMA_DB=$DATADIR/trinity/agalma_trinity.sqlite

cd $DATADIR/trinity/data
#agalma testdata
agalma catalog insert --id SRX288276 --paths SRX288276_1.fq SRX288276_2.fq --species "Abylopsis tetragona" --ncbi_id 316209
agalma catalog insert --id SRX288285 --paths SRX288285_1.fq SRX288285_2.fq --species "Agalma elegans" --ncbi_id 316166
agalma catalog insert --id SRX288432 --paths SRX288432_1.fq SRX288432_2.fq --species "Craseoa lathetica" --ncbi_id 316205
agalma catalog insert --id SRX288431 --paths SRX288431_1.fq SRX288431_2.fq --species "Physalia physalis" --ncbi_id 168775
agalma catalog insert --id SRX288430 --paths SRX288430_1.fq SRX288430_2.fq --species "Nanomia bijuga" --ncbi_id 168759
agalma catalog insert --id JGI_NEMVEC --paths JGI_NEMVEC.fa --species "Nematostella vectensis" --ncbi_id 45351
agalma catalog insert --id NCBI_HYDMAG --paths NCBI_HYDMAG.pfa --species "Hydra magnipapillata" --ncbi_id 6085

cd $DATADIR/trinity/scratch
agalma transcriptome --id SRX288276
agalma transcriptome --id SRX288285
agalma transcriptome --id SRX288430
agalma transcriptome --id SRX288431
agalma transcriptome --id SRX288432

agalma import --id JGI_NEMVEC
agalma translate --id JGI_NEMVEC
agalma annotate --id JGI_NEMVEC

agalma import --id NCBI_HYDMAG --seq_type aa
agalma annotate --id NCBI_HYDMAG

cd $DATADIR/trinity/reports
agalma report --id SRX288276
agalma report --id SRX288285
agalma report --id SRX288430
agalma report --id SRX288431
agalma report --id SRX288432