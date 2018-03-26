#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -c 16
#SBATCH --mem=60G
#SBATCH --array=2-5
# don't forget to change this back to 1-5

source activate corset

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity
mkdir -p /users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
CORSETDIR=/users/aguang/data/aguang/corset
cd $CORSETDIR

# transcriptome SRX IDs
srx=( SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )
fastqs=( SRR871525 SRR871526 SRR871527 SRR871528 SRR871529 )

# transcriptome ids
ids=( 20 21 23 24 22 )

ID=${ids[$SLURM_ARRAY_TASK_ID-1]}
fastq=${fastqs[$SLURM_ARRAY_TASK_ID-1]}
srx_id=${srx[$SLURM_ARRAY_TASK_ID-1]}

#bowtie2-build $DATADIR/scratch/transcriptome-$ID/trinity_out_dir/Trinity.fasta $CORSETDIR/transcriptome-$ID
R1=$DATADIR/data/${fastq}_1.fastq
R2=$DATADIR/data/${fastq}_2.fastq
bowtie2 -p 16 --all -x $CORSETDIR/transcriptome-$ID -1 $R1 -2 $R2 > ${srx_id}.sam  
samtools view -S -b ${srx_id}.sam > ${srx_id}.bam
rm ${srx_id}.sam # if we don't do this we will be way over quota
corset ${srx_id}.bam -p $srx_id