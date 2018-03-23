#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=24G

source activate corset

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity
mkdir -p /users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
CORSETDIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
cd $CORSETDIR

# transcriptome SRX IDs
fastq=( SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )

bowtie
for i in {0..4}
do
    bowtie2-build $DATADIR/scratch/transcriptome-$((i+1))/trinity_out_dir/Trinity.fasta $CORSETDIR/transcriptome-$((i+1))
    FILES=`ls $DATADIR/data/${fastq[$i]}_1.fq | sed 's/_1.fq//g'`
    for F in $FILES ; do
        R1=${F}_1.fq
        R2=${F}_2.fq
        bowtie2 --all -x $CORSETDIR/transcriptome-$((i+1)) -1 $R1 -2 $R2 > ${F}.sam  
        samtools view -S -b ${F}.sam > ${F}.bam
	corset $DATADIR/data/${F}.bam -p ${fastq[$i]}
    done
done