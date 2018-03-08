#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=24G

module load conda
conda create -n corset
conda config --add bioconda
conda install -n corset samtools
conda install -n corset corset
conda install -n corset bowtie2

source activate corset

$DATADIR=/gpfs/scratch/bguang/treeinform/revisions/trinity/scratch
mkdir -p /users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
$CORSETDIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
cd $CORSETDIR

#bowtie
for i in {1..5}
do
    bowtie2-build $DATADIR/transcriptome-${i}/trinity_out_dir/Trinity.fasta $CORSETDIR/transcriptome-${i}
    FILES=`ls $DATADIR/SRX*_P1.fastq | sed 's/_P1.fastq//g'` # need to work on this, also it is fastq data
    for F in $FILES ; do
        R1=${F}_P1.fastq
        R2=${F}_P2.fastq
        bowtie2 --all -S Trinity -1 $R1 -2 $R2 > ${F}.sam  
        samtools view -S -b ${F}.sam > ${F}.bam
    done
done

#corset
corset -g 1,1,1,2,2,2 -n A1,A2,A3,B1,B2,B3 *.bam