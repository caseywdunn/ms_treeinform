#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=24G

module load conda
source activate corset

DATADIR=/gpfs/scratch/aguang/treeinform/revisions/trinity
mkdir -p /users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
CORSETDIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
cd $CORSETDIR

# transcriptome SRX IDs
fastq=( SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )

#bowtie
for i in {1..5}
do
    bowtie2-build $DATADIR/scratch/transcriptome-${i}/trinity_out_dir/Trinity.fasta $CORSETDIR/transcriptome-${i}
    FILES=`ls $DATADIR/data/${fastq[$i]}_1.fq | sed 's/_1.fq//g'` # need to work on this, also it is fastq data
    for F in $FILES ; do
        R1=${F}_1.fq
        R2=${F}_2.fq
        bowtie2 --all -S $CORSETDIR/transcriptome-${i} -1 $R1 -2 $R2 > ${F}.sam  
        samtools view -S -b ${F}.sam > ${F}.bam
    done
done

#corset
corset $CORSETDIR/*.bam